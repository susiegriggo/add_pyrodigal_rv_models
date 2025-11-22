#!/usr/bin/env python3
"""
Download RefSeq RNA virus genomes and annotations using NCBI Datasets.
"""

import argparse
import os
import random
import shutil
import subprocess


def check_datasets_installed():
    """Check if NCBI datasets CLI is installed."""
    print("Checking for NCBI datasets CLI...")
    if shutil.which("datasets") is None:
        print("❌ Error: NCBI datasets CLI not found!")
        print("Install: conda install -c conda-forge ncbi-datasets-cli")
        exit(1)
    print("✓ NCBI datasets CLI found\n")


def parse_taxdump_fast(taxdump_dir="./", num_needed=10):
    """Parse taxdump files to find RNA virus species (optimized)."""
    nodes_file = os.path.join(taxdump_dir, "nodes.dmp")
    names_file = os.path.join(taxdump_dir, "names.dmp")

    print(f"Looking for taxdump files in {taxdump_dir}...")
    if not os.path.exists(nodes_file) or not os.path.exists(names_file):
        print(f"❌ Error: taxdump files not found in {taxdump_dir}")
        print(
            "Download with: wget https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz"
        )
        exit(1)
    print("✓ Found nodes.dmp and names.dmp\n")

    print("Parsing taxonomy tree (nodes.dmp)...")
    taxid_to_parent = {}
    taxid_to_rank = {}

    with open(nodes_file, "r") as f:
        for line in f:
            parts = line.strip().split("|")
            taxid = parts[0].strip()
            parent = parts[1].strip()
            rank = parts[2].strip()
            taxid_to_parent[taxid] = parent
            taxid_to_rank[taxid] = rank

    print(f"  Loaded {len(taxid_to_parent)} taxa")

    # Find RNA virus descendants
    print("Finding RNA virus species (taxid 2559587)...")
    rna_virus_taxid = "2559587"

    # Build path from each taxid to root for efficient lookup
    def is_rna_virus(taxid):
        current = taxid
        while current != "1" and current in taxid_to_parent:
            if current == rna_virus_taxid:
                return True
            current = taxid_to_parent[current]
        return False

    # Find species-level RNA virus taxa
    species_taxids = []
    for taxid, rank in taxid_to_rank.items():
        if rank == "species" and is_rna_virus(taxid):
            species_taxids.append(taxid)

    print(f"  Found {len(species_taxids)} RNA virus species")

    # Randomly sample the species we need
    num_to_sample = min(num_needed, len(species_taxids))
    selected_taxids = random.sample(species_taxids, num_to_sample)
    selected_set = set(selected_taxids)

    print(f"  Randomly selected {num_to_sample} species")

    # Get names only for selected species
    print("Loading names for selected species only...")
    taxid_to_name = {}
    with open(names_file, "r") as f:
        for line in f:
            parts = line.strip().split("|")
            taxid = parts[0].strip()
            name = parts[1].strip()
            name_class = parts[3].strip()

            if name_class == "scientific name" and taxid in selected_set:
                taxid_to_name[taxid] = name
                # Early exit once we have all names
                if len(taxid_to_name) == len(selected_set):
                    break

    print(f"✓ Loaded names for {len(taxid_to_name)} species\n")

    return [(taxid, taxid_to_name.get(taxid, "Unknown")) for taxid in selected_taxids]


def download_genomes(taxids, output_dir):
    """Download genomes using taxid input file."""
    # Create input file with taxids
    input_file = os.path.join(output_dir, "taxids.txt")
    print(f"Creating taxid input file: {input_file}")
    with open(input_file, "w") as f:
        for taxid in taxids:
            f.write(f"{taxid}\n")
    print(f"✓ Wrote {len(taxids)} taxids to file\n")

    print(f"Downloading {len(taxids)} genomes from NCBI...")
    print("This may take a few minutes depending on genome sizes...\n")

    cmd = [
        "datasets",
        "download",
        "virus",
        "genome",
        "taxon",
        "--inputfile",
        input_file,
        "--include",
        "genome,annotation",
        "--refseq",
        "--filename",
        os.path.join(output_dir, "dataset.zip"),
    ]

    result = subprocess.run(cmd, capture_output=True, text=True)

    if result.returncode != 0:
        print(f"❌ Download error: {result.stderr}")
        return

    print("✓ Download complete!\n")

    print("Extracting files...")
    subprocess.run(
        ["unzip", "-q", "-o", os.path.join(output_dir, "dataset.zip"), "-d", output_dir]
    )
    print(f"✓ Extraction complete!\n")

    # Clean up input file
    print("Cleaning up temporary files...")
    os.remove(input_file)

    # Show what was downloaded
    data_dir = os.path.join(output_dir, "ncbi_dataset", "data")
    if os.path.exists(data_dir):
        genome_count = len(
            [
                d
                for d in os.listdir(data_dir)
                if os.path.isdir(os.path.join(data_dir, d))
            ]
        )
        print(f"✓ Downloaded {genome_count} genome(s)")
        print(f"\nFiles located in: {data_dir}/")
        print("\nEach genome directory contains:")
        print("  - *.fna (genome sequence)")
        print("  - genomic.gff (annotations)")
    else:
        print("⚠ Warning: Expected data directory not found")


def main():
    parser = argparse.ArgumentParser(description="Download RefSeq RNA virus genomes")
    parser.add_argument(
        "-n",
        "--num-genomes",
        type=int,
        default=10,
        help="Number of genomes (default: 10)",
    )
    parser.add_argument(
        "-o",
        "--output-dir",
        default="rna_virus_genomes",
        help="Output directory (default: rna_virus_genomes)",
    )
    parser.add_argument(
        "-t", "--taxdump-dir", default="./", help="Taxdump directory (default: ./)"
    )

    args = parser.parse_args()

    print("=" * 60)
    print("RNA Virus Genome Downloader")
    print("=" * 60 + "\n")

    check_datasets_installed()

    print(f"Creating output directory: {args.output_dir}")
    os.makedirs(args.output_dir, exist_ok=True)
    print(f"✓ Output directory ready\n")

    # Get species (now much faster)
    selected = parse_taxdump_fast(args.taxdump_dir, args.num_genomes)

    print(f"\n{'='*60}")
    print(f"Selected {len(selected)} species:")
    print(f"{'='*60}")
    for i, (taxid, name) in enumerate(selected, 1):
        print(f"{i:2d}. [{taxid}] {name}")
    print(f"{'='*60}\n")

    # Download
    taxids = [taxid for taxid, _ in selected]
    download_genomes(taxids, args.output_dir)

    print("\n" + "=" * 60)
    print("✓ All done!")
    print("=" * 60)


if __name__ == "__main__":
    main()
