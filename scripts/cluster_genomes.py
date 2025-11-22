#!/usr/bin/env python3
"""
Cluster genome sequences using MMseqs2 easy-cluster.

Reduces sequence redundancy by clustering at specified identity and coverage
thresholds, outputting representative sequences for each cluster.
"""

import argparse
import multiprocessing
import shutil
import subprocess
import sys
from pathlib import Path


def check_mmseqs_installed():
    """Check if MMseqs2 is installed"""
    try:
        result = subprocess.run(
            ["mmseqs", "version"], capture_output=True, text=True, check=False
        )
        return result.returncode == 0
    except FileNotFoundError:
        return False


def count_sequences(fasta_file):
    """Count sequences in a FASTA file"""
    try:
        with open(fasta_file, "r") as f:
            return sum(1 for line in f if line.startswith(">"))
    except Exception:
        return 0


def cluster_genomes(
    input_fasta,
    output_fasta,
    identity=0.8,
    coverage=0.8,
    threads=None,
    cov_mode=0,
    quiet=False,
):
    """Cluster sequences using MMseqs2 easy-cluster and extract representatives"""
    input_fasta = Path(input_fasta)

    if not input_fasta.exists():
        if not quiet:
            print(f"Warning: {input_fasta} does not exist, skipping...")
        return False

    # Check if file is empty or has no sequences
    original_count = count_sequences(input_fasta)
    if original_count == 0:
        if not quiet:
            print(f"Warning: {input_fasta.name} contains no sequences, skipping...")
        return False

    if not quiet:
        print(
            f"Clustering {input_fasta.name} at {identity*100}% identity, {coverage*100}% coverage..."
        )
        print(f"  Input sequences: {original_count}")

    # Create temporary directory for MMseqs2 output
    temp_dir = Path(f"temp_mmseqs_{input_fasta.stem}")
    temp_dir.mkdir(exist_ok=True)

    try:
        # Run MMseqs2 easy-cluster
        cluster_prefix = temp_dir / "cluster"
        cmd = [
            "mmseqs",
            "easy-cluster",
            str(input_fasta),
            str(cluster_prefix),
            str(temp_dir / "tmp"),
            "--min-seq-id",
            str(identity),
            "-c",
            str(coverage),
            "--cov-mode",
            str(cov_mode),
        ]

        # Add threads if specified
        if threads:
            cmd.extend(["--threads", str(threads)])

        result = subprocess.run(cmd, capture_output=True, text=True)

        if result.returncode != 0:
            if not quiet:
                print(f"  Error running MMseqs2: {result.stderr}")
            return False

        # Get representative sequences file
        rep_seq_file = temp_dir / "cluster_rep_seq.fasta"

        if not rep_seq_file.exists():
            if not quiet:
                print(f"  Error: Representative sequences file not found")
            return False

        # Create output directory if needed
        output_path = Path(output_fasta)
        output_path.parent.mkdir(parents=True, exist_ok=True)

        # Copy representative sequences to output
        shutil.copy(rep_seq_file, output_fasta)

        # Count clustered sequences
        clustered_count = count_sequences(output_fasta)

        if not quiet:
            reduction = (
                ((original_count - clustered_count) / original_count * 100)
                if original_count > 0
                else 0
            )
            print(
                f"  Clustered to {clustered_count} representatives ({reduction:.1f}% reduction)"
            )
            print(f"  Output: {output_fasta}")

        return True

    except Exception as e:
        if not quiet:
            print(f"  Error processing {input_fasta.name}: {e}", file=sys.stderr)
        return False

    finally:
        # Clean up temporary directory
        if temp_dir.exists():
            shutil.rmtree(temp_dir)


def process_directory(
    input_dir, output_dir, identity, coverage, threads, cov_mode, quiet
):
    """Process all FASTA files in a directory"""
    input_path = Path(input_dir)

    if not input_path.exists():
        print(f"Error: Input directory '{input_dir}' does not exist", file=sys.stderr)
        sys.exit(1)

    if not input_path.is_dir():
        print(f"Error: '{input_dir}' is not a directory", file=sys.stderr)
        sys.exit(1)

    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)

    # Find all FASTA files
    fasta_files = sorted(
        list(input_path.glob("*.fasta"))
        + list(input_path.glob("*.fa"))
        + list(input_path.glob("*.fna"))
    )

    if not fasta_files:
        print(f"Warning: No FASTA files found in {input_dir}")
        return

    if not quiet:
        print(f"Found {len(fasta_files)} FASTA file(s) to process")
        print(
            f"Parameters: identity={identity*100}%, coverage={coverage*100}%, threads={threads or 'auto'}\n"
        )

    success_count = 0
    for fasta_file in fasta_files:
        output_fasta = output_path / fasta_file.name
        if cluster_genomes(
            fasta_file, output_fasta, identity, coverage, threads, cov_mode, quiet
        ):
            success_count += 1
        if not quiet:
            print()

    if not quiet:
        print(
            f"Processing complete. {success_count}/{len(fasta_files)} file(s) clustered successfully."
        )


def process_single_file(
    input_file, output_file, identity, coverage, threads, cov_mode, quiet
):
    """Process a single FASTA file"""
    if not Path(input_file).exists():
        print(f"Error: Input file '{input_file}' does not exist", file=sys.stderr)
        sys.exit(1)

    success = cluster_genomes(
        input_file, output_file, identity, coverage, threads, cov_mode, quiet
    )

    if not quiet:
        if success:
            print("\nClustering complete.")
        else:
            print("\nClustering failed.")

    sys.exit(0 if success else 1)


def main():
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )

    parser.add_argument(
        "-i", "--input", help="Input FASTA file or directory containing FASTA files"
    )

    parser.add_argument(
        "-o",
        "--output",
        help="Output FASTA file (for single file) or directory (for directory mode)",
    )

    parser.add_argument(
        "--input-dir",
        default="genomes_output",
        help="Input directory containing FASTA files (default: genomes_output)",
    )

    parser.add_argument(
        "--output-dir",
        default="genomes_clustered",
        help="Output directory for clustered FASTA files (default: genomes_clustered)",
    )

    parser.add_argument(
        "--identity",
        type=float,
        default=0.8,
        help="Minimum sequence identity threshold (0.0-1.0, default: 0.8)",
    )

    parser.add_argument(
        "--coverage",
        type=float,
        default=0.8,
        help="Minimum coverage threshold (0.0-1.0, default: 0.8)",
    )

    parser.add_argument(
        "--cov-mode",
        type=int,
        choices=[0, 1, 2, 3, 4, 5],
        default=0,
        help="Coverage mode: 0=bidirectional, 1=target, 2=query, 3=target/query, 4=query/target, 5=length (default: 0)",
    )

    parser.add_argument(
        "-t",
        "--threads",
        type=int,
        help=f"Number of threads to use (default: all available = {multiprocessing.cpu_count()})",
    )

    parser.add_argument(
        "-q", "--quiet", action="store_true", help="Suppress progress messages"
    )

    args = parser.parse_args()

    # Validate identity and coverage
    if not 0.0 <= args.identity <= 1.0:
        print("Error: --identity must be between 0.0 and 1.0", file=sys.stderr)
        sys.exit(1)

    if not 0.0 <= args.coverage <= 1.0:
        print("Error: --coverage must be between 0.0 and 1.0", file=sys.stderr)
        sys.exit(1)

    # Check if MMseqs2 is installed
    if not check_mmseqs_installed():
        print("Error: MMseqs2 is not installed or not in PATH", file=sys.stderr)
        print("Install from: https://github.com/soedinglab/MMseqs2", file=sys.stderr)
        sys.exit(1)

    # Determine mode: single file or directory
    if args.input and args.output:
        # Single file mode
        process_single_file(
            args.input,
            args.output,
            args.identity,
            args.coverage,
            args.threads,
            args.cov_mode,
            args.quiet,
        )
    elif args.input:
        # Directory mode with custom input
        input_path = Path(args.input)
        if input_path.is_file():
            print("Error: --input is a file but no --output specified", file=sys.stderr)
            print("Use: -i <file> -o <output.fasta>", file=sys.stderr)
            sys.exit(1)
        output_dir = args.output if args.output else args.output_dir
        process_directory(
            args.input,
            output_dir,
            args.identity,
            args.coverage,
            args.threads,
            args.cov_mode,
            args.quiet,
        )
    else:
        # Default directory mode
        process_directory(
            args.input_dir,
            args.output_dir,
            args.identity,
            args.coverage,
            args.threads,
            args.cov_mode,
            args.quiet,
        )


if __name__ == "__main__":
    main()
