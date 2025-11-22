#!/usr/bin/env python3
"""
Create a TSV file with taxonomy information from taxid files.

Reads taxonomy IDs from files and extracts detailed information from NCBI
taxdump files (nodes.dmp, names.dmp) including family, rank, scientific name,
genetic code, and lineage.
"""

import argparse
import csv
import sys
from pathlib import Path


def load_nodes(nodes_file):
    """Load taxonomy nodes from nodes.dmp"""
    tax_to_parent = {}
    tax_to_rank = {}
    tax_to_gencode = {}

    try:
        with open(nodes_file, "r") as f:
            for line in f:
                parts = [p.strip() for p in line.split("|")]
                tax_id = int(parts[0])
                parent_tax_id = int(parts[1])
                rank = parts[2]
                gencode = int(parts[6])

                tax_to_parent[tax_id] = parent_tax_id
                tax_to_rank[tax_id] = rank
                tax_to_gencode[tax_id] = gencode
    except FileNotFoundError:
        print(f"Error: nodes.dmp file '{nodes_file}' not found", file=sys.stderr)
        sys.exit(1)
    except Exception as e:
        print(f"Error loading nodes.dmp: {e}", file=sys.stderr)
        sys.exit(1)

    return tax_to_parent, tax_to_rank, tax_to_gencode


def load_names(names_file):
    """Load scientific names from names.dmp"""
    tax_to_name = {}

    try:
        with open(names_file, "r") as f:
            for line in f:
                parts = [p.strip() for p in line.split("|")]
                tax_id = int(parts[0])
                name = parts[1]
                name_class = parts[3]

                # Only keep scientific names
                if name_class == "scientific name":
                    tax_to_name[tax_id] = name
    except FileNotFoundError:
        print(f"Error: names.dmp file '{names_file}' not found", file=sys.stderr)
        sys.exit(1)
    except Exception as e:
        print(f"Error loading names.dmp: {e}", file=sys.stderr)
        sys.exit(1)

    return tax_to_name


def get_lineage(tax_id, tax_to_parent, tax_to_name):
    """Get full lineage for a taxid"""
    lineage_parts = []
    current = tax_id

    # Traverse up to root
    visited = set()
    while current != 1 and current in tax_to_parent:
        if current in visited:  # Prevent infinite loops
            break
        visited.add(current)

        name = tax_to_name.get(current, f"unknown_{current}")
        lineage_parts.append(name)

        current = tax_to_parent[current]

    # Reverse to get root -> leaf order and join with semicolons
    return ";".join(reversed(lineage_parts))


def get_family_from_lineage(tax_id, tax_to_parent, tax_to_rank, tax_to_name):
    """Get the family name from lineage by traversing up"""
    current = tax_id
    visited = set()

    while current != 1 and current in tax_to_parent:
        if current in visited:
            break
        visited.add(current)

        rank = tax_to_rank.get(current, "")
        if rank == "family":
            return tax_to_name.get(current, f"unknown_family_{current}")

        current = tax_to_parent[current]

    return "no_family"


def extract_family_from_filename(filename):
    """Extract family name from filename (e.g., 'Familyname_11.txt' -> 'Familyname')"""
    stem = Path(filename).stem
    # Remove genetic code suffix if present (pattern: _\d+)
    parts = stem.rsplit("_", 1)
    if len(parts) == 2 and parts[1].isdigit():
        return parts[0]
    return stem


def process_taxid_files(
    input_dir, output_tsv, nodes_file, names_file, use_filename_family, quiet=False
):
    """Process all taxid files and create TSV"""

    input_path = Path(input_dir)

    if not input_path.exists():
        print(f"Error: Input directory '{input_dir}' does not exist", file=sys.stderr)
        sys.exit(1)

    if not input_path.is_dir():
        print(f"Error: '{input_dir}' is not a directory", file=sys.stderr)
        sys.exit(1)

    # Load taxonomy data
    if not quiet:
        print("Loading NCBI taxonomy data...")

    tax_to_parent, tax_to_rank, tax_to_gencode = load_nodes(nodes_file)
    tax_to_name = load_names(names_file)

    if not quiet:
        print(f"Loaded {len(tax_to_parent)} taxonomy nodes")
        print(f"Loaded {len(tax_to_name)} scientific names\n")

    # Find all taxid files
    taxid_files = sorted(input_path.glob("*.txt"))

    if not taxid_files:
        print(f"Warning: No .txt files found in {input_dir}")
        return

    if not quiet:
        print(f"Found {len(taxid_files)} taxid file(s) to process\n")

    # Collect all data
    rows = []
    total_taxids = 0
    skipped_taxids = 0

    for taxid_file in taxid_files:
        if not quiet:
            print(f"Processing {taxid_file.name}...")

        # Get family name from filename if requested
        if use_filename_family:
            family_from_file = extract_family_from_filename(taxid_file)
        else:
            family_from_file = None

        # Read taxids from file
        try:
            with open(taxid_file, "r") as f:
                taxids = [line.strip() for line in f if line.strip()]
        except Exception as e:
            print(f"  Error reading file: {e}", file=sys.stderr)
            continue

        if not quiet:
            print(f"  Found {len(taxids)} taxid(s)")

        # Process each taxid
        for taxid_str in taxids:
            try:
                taxid = int(taxid_str)
            except ValueError:
                if not quiet:
                    print(f"  Warning: Invalid taxid '{taxid_str}', skipping")
                skipped_taxids += 1
                continue

            total_taxids += 1

            # Check if taxid exists in taxonomy
            if taxid not in tax_to_rank:
                if not quiet:
                    print(f"  Warning: taxid {taxid} not found in taxonomy, skipping")
                skipped_taxids += 1
                continue

            # Get taxonomy information
            rank = tax_to_rank.get(taxid, "unknown")
            scientific_name = tax_to_name.get(taxid, f"unknown_{taxid}")
            genetic_code = tax_to_gencode.get(taxid, 1)  # Default to standard code
            lineage = get_lineage(taxid, tax_to_parent, tax_to_name)

            # Determine family
            if use_filename_family:
                family = family_from_file
            else:
                family = get_family_from_lineage(
                    taxid, tax_to_parent, tax_to_rank, tax_to_name
                )

            # Add row
            rows.append(
                {
                    "family": family,
                    "taxid": taxid,
                    "rank": rank,
                    "scientific_name": scientific_name,
                    "genetic_code": genetic_code,
                    "lineage": lineage,
                }
            )

        if not quiet:
            print()

    if not rows:
        print("Warning: No valid taxids found to process")
        return

    # Write to TSV
    output_path = Path(output_tsv)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    try:
        with open(output_path, "w", newline="") as f:
            writer = csv.DictWriter(
                f,
                fieldnames=[
                    "family",
                    "taxid",
                    "rank",
                    "scientific_name",
                    "genetic_code",
                    "lineage",
                ],
                delimiter="\t",
            )
            writer.writeheader()
            writer.writerows(rows)

        if not quiet:
            print(f"Successfully wrote {len(rows)} row(s) to {output_path}")
            print(f"Total taxids processed: {total_taxids}")
            if skipped_taxids > 0:
                print(f"Skipped taxids: {skipped_taxids}")

    except Exception as e:
        print(f"Error writing output file: {e}", file=sys.stderr)
        sys.exit(1)


def main():
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )

    parser.add_argument(
        "-i",
        "--input-dir",
        default="taxid_files",
        help="Input directory containing taxid files (default: taxid_files)",
    )

    parser.add_argument(
        "-o",
        "--output",
        default="extra_models.tsv",
        help="Output TSV file (default: extra_models.tsv)",
    )

    parser.add_argument(
        "--nodes",
        default="nodes.dmp",
        help="Path to NCBI nodes.dmp file (default: nodes.dmp)",
    )

    parser.add_argument(
        "--names",
        default="names.dmp",
        help="Path to NCBI names.dmp file (default: names.dmp)",
    )

    parser.add_argument(
        "--use-filename-family",
        action="store_true",
        help="Extract family name from filename instead of lineage",
    )

    parser.add_argument(
        "-q", "--quiet", action="store_true", help="Suppress progress messages"
    )

    args = parser.parse_args()

    process_taxid_files(
        args.input_dir,
        args.output,
        args.nodes,
        args.names,
        args.use_filename_family,
        args.quiet,
    )


if __name__ == "__main__":
    main()
