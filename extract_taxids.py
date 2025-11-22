#!/usr/bin/env python3
"""
Extract taxonomy IDs grouped by family and genetic code.

For families with multiple genetic codes, filters descendants by their genetic code.
For families with a single genetic code, uses the family taxid directly.
"""

import argparse
import sys
from pathlib import Path

import polars as pl


def load_taxonomy(nodes_file):
    """Load taxonomy tree from nodes.dmp"""
    tax_to_parent = {}
    tax_to_gencode = {}

    try:
        with open(nodes_file, "r") as f:
            for line in f:
                parts = [p.strip() for p in line.split("|")]
                tax_id = int(parts[0])
                parent_tax_id = int(parts[1])
                gencode = int(parts[6])

                tax_to_parent[tax_id] = parent_tax_id
                tax_to_gencode[tax_id] = gencode
    except FileNotFoundError:
        print(f"Error: Taxonomy file '{nodes_file}' not found", file=sys.stderr)
        sys.exit(1)
    except Exception as e:
        print(f"Error loading taxonomy: {e}", file=sys.stderr)
        sys.exit(1)

    return tax_to_parent, tax_to_gencode


def get_all_descendants(tax_id, tax_to_parent):
    """Get all descendant taxids of a given taxid"""
    # Build parent-to-children mapping
    children = {}
    for child, parent in tax_to_parent.items():
        if parent not in children:
            children[parent] = []
        children[parent].append(child)

    # Recursively get all descendants
    descendants = set()

    def get_descendants_recursive(tid):
        descendants.add(tid)
        if tid in children:
            for child in children[tid]:
                get_descendants_recursive(child)

    get_descendants_recursive(tax_id)
    return descendants


def process_taxonomy(input_file, nodes_file, output_dir, quiet=False):
    """Main processing function"""
    # Load taxonomy
    if not quiet:
        print("Loading taxonomy...")
    tax_to_parent, tax_to_gencode = load_taxonomy(nodes_file)
    if not quiet:
        print(f"Loaded {len(tax_to_parent)} taxonomy nodes\n")

    # Read the TSV
    try:
        df = pl.read_csv(input_file, separator="\t")
    except FileNotFoundError:
        print(f"Error: Input file '{input_file}' not found", file=sys.stderr)
        sys.exit(1)
    except Exception as e:
        print(f"Error reading input file: {e}", file=sys.stderr)
        sys.exit(1)

    # Create output directory
    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)

    # Find families with multiple genetic codes
    families_with_multiple_gencodes = set(
        df.group_by("family")
        .agg(pl.col("genetic_code").n_unique().alias("gencode_count"))
        .filter(pl.col("gencode_count") > 1)["family"]
        .to_list()
    )

    if not quiet:
        print(
            f"Families with multiple genetic codes: {len(families_with_multiple_gencodes)}\n"
        )

    files_written = 0
    # Group by family and genetic code
    for (family_name, genetic_code), group_df in df.group_by(
        ["family", "genetic_code"]
    ):
        output_file = output_path / f"{family_name}_{genetic_code}.txt"

        # Get family taxid (should only be one per family-genetic_code combo)
        family_taxid = group_df.filter(pl.col("rank") == "family")["taxid"].to_list()

        if family_taxid and family_name in families_with_multiple_gencodes:
            # Family with multiple genetic codes: get all descendants and filter by genetic code
            all_descendants = get_all_descendants(family_taxid[0], tax_to_parent)
            taxids = [
                tid
                for tid in all_descendants
                if tax_to_gencode.get(tid, 1) == genetic_code
            ]
        elif family_taxid:
            # Single genetic code family: just use the family taxid
            taxids = family_taxid
        else:
            # Non-family entries: use the taxids directly
            taxids = group_df["taxid"].to_list()

        # Write to file
        taxids = sorted(set(taxids))
        with open(output_file, "w") as f:
            for taxid in taxids:
                f.write(f"{taxid}\n")

        if not quiet:
            print(f"Saved {len(taxids)} taxids to {output_file}")
        files_written += 1

    if not quiet:
        print(f"\nProcessing complete. {files_written} files written to {output_dir}")


def main():
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )

    parser.add_argument(
        "-i",
        "--input",
        default="riboviria_genetic_codes.tsv",
        help="Input TSV file with taxonomy and genetic code data (default: riboviria_genetic_codes.tsv)",
    )

    parser.add_argument(
        "-n",
        "--nodes",
        default="nodes.dmp",
        help="NCBI taxonomy nodes.dmp file (default: nodes.dmp)",
    )

    parser.add_argument(
        "-o",
        "--output",
        default="taxid_files",
        help="Output directory for taxid files (default: taxid_files)",
    )

    parser.add_argument(
        "-q", "--quiet", action="store_true", help="Suppress progress messages"
    )

    args = parser.parse_args()

    process_taxonomy(
        input_file=args.input,
        nodes_file=args.nodes,
        output_dir=args.output,
        quiet=args.quiet,
    )


if __name__ == "__main__":
    main()
