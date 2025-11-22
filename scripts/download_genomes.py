#!/usr/bin/env python3
"""
Download and merge viral genomes using NCBI datasets CLI.

Downloads genomes for taxonomy IDs listed in input files and merges them into
single FASTA files. Processes downloads in batches to handle API limits.
"""

import argparse
import shutil
import subprocess
import sys
import zipfile
from pathlib import Path


def download_and_merge_genomes(taxid_file, output_fasta, batch_size=100, quiet=False):
    """Download genomes using --inputfile and merge into single FASTA"""
    taxid_file = Path(taxid_file)

    if not taxid_file.exists():
        if not quiet:
            print(f"Warning: {taxid_file} does not exist, skipping...")
        return False

    if not quiet:
        print(f"Processing {taxid_file.name}...")

    # Read taxids
    with open(taxid_file, "r") as f:
        taxids = [line.strip() for line in f if line.strip()]

    if not taxids:
        if not quiet:
            print(f"  No taxids to download")
        return False

    if not quiet:
        print(f"  Found {len(taxids)} taxid(s) to download")

    # Split into batches
    batches = [taxids[i : i + batch_size] for i in range(0, len(taxids), batch_size)]

    if not quiet:
        print(f"  Processing {len(batches)} batch(es)...")

    temp_base = Path(f"temp_{taxid_file.stem}")
    temp_base.mkdir(exist_ok=True)

    try:
        # Download all batches
        for batch_num, batch in enumerate(batches, 1):
            if not quiet:
                print(f"  Batch {batch_num}/{len(batches)} ({len(batch)} taxids)...")

            batch_file = temp_base / f"batch_{batch_num}.txt"
            with open(batch_file, "w") as f:
                for taxid in batch:
                    f.write(f"{taxid}\n")

            zip_file = temp_base / f"batch_{batch_num}.zip"
            cmd = [
                "datasets",
                "download",
                "virus",
                "genome",
                "taxon",
                "--inputfile",
                str(batch_file),
                "--filename",
                str(zip_file),
            ]

            result = subprocess.run(cmd, capture_output=True, text=True)

            if result.returncode == 0 and zip_file.exists():
                batch_dir = temp_base / f"batch_{batch_num}"
                with zipfile.ZipFile(zip_file, "r") as zip_ref:
                    zip_ref.extractall(batch_dir)
                zip_file.unlink()
            else:
                if not quiet:
                    print(f"    Warning: Failed to download batch {batch_num}")
                    if result.stderr:
                        print(f"    Error: {result.stderr}")

            batch_file.unlink()

        # Merge all FASTA files from all batches
        fasta_files = sorted(temp_base.rglob("*.fna"))

        if not fasta_files:
            if not quiet:
                print(f"  Warning: No FASTA files found to merge")
            return False

        with open(output_fasta, "w") as outfile:
            for fasta_file in fasta_files:
                with open(fasta_file, "r") as infile:
                    outfile.write(infile.read())

        if not quiet:
            print(f"  Saved to {output_fasta}")

        return True

    except Exception as e:
        print(f"  Error processing {taxid_file.name}: {e}", file=sys.stderr)
        return False

    finally:
        if temp_base.exists():
            shutil.rmtree(temp_base)


def check_datasets_installed():
    """Check if NCBI datasets CLI is installed"""
    try:
        result = subprocess.run(
            ["datasets", "--version"], capture_output=True, text=True, check=False
        )
        return result.returncode == 0
    except FileNotFoundError:
        return False


def process_directory(input_dir, output_dir, batch_size, quiet):
    """Process all taxid files in a directory"""
    input_path = Path(input_dir)

    if not input_path.exists():
        print(f"Error: Input directory '{input_dir}' does not exist", file=sys.stderr)
        sys.exit(1)

    if not input_path.is_dir():
        print(f"Error: '{input_dir}' is not a directory", file=sys.stderr)
        sys.exit(1)

    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)

    taxid_files = sorted(input_path.glob("*.txt"))

    if not taxid_files:
        print(f"Warning: No .txt files found in {input_dir}")
        return

    if not quiet:
        print(f"Found {len(taxid_files)} taxid file(s) to process\n")

    success_count = 0
    for taxid_file in taxid_files:
        output_fasta = output_path / f"{taxid_file.stem}.fasta"
        if download_and_merge_genomes(taxid_file, output_fasta, batch_size, quiet):
            success_count += 1
        if not quiet:
            print()

    if not quiet:
        print(
            f"Processing complete. {success_count}/{len(taxid_files)} file(s) processed successfully."
        )


def process_single_file(input_file, output_file, batch_size, quiet):
    """Process a single taxid file"""
    if not Path(input_file).exists():
        print(f"Error: Input file '{input_file}' does not exist", file=sys.stderr)
        sys.exit(1)

    # Create output directory if needed
    output_path = Path(output_file)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    success = download_and_merge_genomes(input_file, output_file, batch_size, quiet)

    if not quiet:
        if success:
            print("\nProcessing complete.")
        else:
            print("\nProcessing completed with warnings.")


def main():
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )

    parser.add_argument(
        "-i", "--input", help="Input taxid file or directory containing taxid files"
    )

    parser.add_argument(
        "-o",
        "--output",
        help="Output FASTA file (for single file) or directory (for directory mode)",
    )

    parser.add_argument(
        "--input-dir",
        default="taxid_files",
        help="Input directory containing taxid files (default: taxid_files)",
    )

    parser.add_argument(
        "--output-dir",
        default="genomes_output",
        help="Output directory for genome FASTA files (default: genomes_output)",
    )

    parser.add_argument(
        "-b",
        "--batch-size",
        type=int,
        default=100,
        help="Number of taxids to download per batch (default: 100)",
    )

    parser.add_argument(
        "-q", "--quiet", action="store_true", help="Suppress progress messages"
    )

    args = parser.parse_args()

    # Check if datasets CLI is installed
    if not check_datasets_installed():
        print(
            "Error: NCBI datasets CLI is not installed or not in PATH", file=sys.stderr
        )
        print(
            "Install from: https://www.ncbi.nlm.nih.gov/datasets/docs/v2/download-and-install/",
            file=sys.stderr,
        )
        sys.exit(1)

    # Determine mode: single file or directory
    if args.input and args.output:
        # Single file mode
        process_single_file(args.input, args.output, args.batch_size, args.quiet)
    elif args.input:
        # Directory mode with custom input
        input_path = Path(args.input)
        if input_path.is_file():
            print("Error: --input is a file but no --output specified", file=sys.stderr)
            print("Use: -i <file> -o <output.fasta>", file=sys.stderr)
            sys.exit(1)
        output_dir = args.output if args.output else args.output_dir
        process_directory(args.input, output_dir, args.batch_size, args.quiet)
    else:
        # Default directory mode
        process_directory(args.input_dir, args.output_dir, args.batch_size, args.quiet)


if __name__ == "__main__":
    main()
