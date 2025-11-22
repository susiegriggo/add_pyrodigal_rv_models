#!/usr/bin/env python3
"""
Train Pyrodigal gene prediction models for viral genomes.

Trains ORF finding models on clustered genome sequences using specified genetic
codes and updates or creates a meta.json file with model metadata.
"""

import argparse
import json
import sys
from pathlib import Path

try:
    import polars as pl
except ImportError:
    print(
        "Error: polars is not installed. Install with: pip install polars",
        file=sys.stderr,
    )
    sys.exit(1)

try:
    import pyrodigal
except ImportError:
    print(
        "Error: pyrodigal is not installed. Install with: pip install pyrodigal",
        file=sys.stderr,
    )
    sys.exit(1)

try:
    import Bio.SeqIO
except ImportError:
    print(
        "Error: biopython is not installed. Install with: pip install biopython",
        file=sys.stderr,
    )
    sys.exit(1)


def compact_json_dump(data, file_handle, indent=2):
    """
    Dump JSON with compact array formatting (keeps matrices on single/few lines)
    """

    def format_value(obj, current_indent=0):
        indent_str = " " * current_indent
        next_indent = " " * (current_indent + indent)

        if isinstance(obj, dict):
            if not obj:
                return "{}"
            items = []
            for key, value in obj.items():
                formatted_value = format_value(value, current_indent + indent)
                items.append(f'{next_indent}"{key}": {formatted_value}')
            return "{\n" + ",\n".join(items) + "\n" + indent_str + "}"

        elif isinstance(obj, list):
            if not obj:
                return "[]"
            # Check if all elements are simple types (int, float, bool, None, short strings)
            is_simple = all(
                isinstance(item, (int, float, bool, type(None)))
                or (isinstance(item, str) and len(item) < 50)
                for item in obj
            )

            if is_simple and len(obj) <= 20:
                # Compact format for simple, short arrays
                return "[" + ", ".join(json.dumps(item) for item in obj) + "]"
            elif is_simple:
                # For longer simple arrays, break into lines of ~10 elements
                chunks = [obj[i : i + 10] for i in range(0, len(obj), 10)]
                lines = [
                    ", ".join(json.dumps(item) for item in chunk) for chunk in chunks
                ]
                return (
                    "[\n"
                    + next_indent
                    + (",\n" + next_indent).join(lines)
                    + "\n"
                    + indent_str
                    + "]"
                )
            else:
                # Regular formatting for complex arrays
                items = [format_value(item, current_indent + indent) for item in obj]
                return (
                    "[\n"
                    + next_indent
                    + (",\n" + next_indent).join(items)
                    + "\n"
                    + indent_str
                    + "]"
                )

        else:
            return json.dumps(obj)

    formatted = format_value(data)
    file_handle.write(formatted)


def count_sequences_and_bases(fasta_file):
    """Count sequences and total bases in a FASTA file"""
    try:
        records = list(Bio.SeqIO.parse(fasta_file, "fasta"))
        num_seqs = len(records)
        total_bases = sum(len(rec.seq) for rec in records)
        return num_seqs, total_bases, records
    except Exception as e:
        return 0, 0, []


def train_pyrodigal_model(
    fasta_file, genetic_code, output_json, min_length=5000, quiet=False
):
    """Train a Pyrodigal model on a FASTA file with specified genetic code"""
    fasta_file = Path(fasta_file)

    if not fasta_file.exists():
        if not quiet:
            print(f"Warning: {fasta_file} does not exist, skipping...")
        return None

    if not quiet:
        print(f"Training model for {fasta_file.name} (genetic code {genetic_code})...")

    # Read sequences from FASTA using Biopython
    num_seqs, total_bases, records = count_sequences_and_bases(fasta_file)

    if not records:
        if not quiet:
            print(f"  Warning: No sequences found in {fasta_file}")
        return None

    if not quiet:
        print(f"  Found {num_seqs} sequence(s), {total_bases:,} bp total")

    if total_bases < min_length:
        if not quiet:
            print(f"  Total length < {min_length:,} bp, skipping training")
        return None

    # Train ORF finder
    try:
        orf_finder = pyrodigal.GeneFinder()
        training_info = orf_finder.train(
            *(bytes(seq.seq) for seq in records), translation_table=genetic_code
        )
        if not quiet:
            print(f"  Successfully trained model")
    except Exception as e:
        if not quiet:
            print(f"  Error: Could not train model: {e}")
        return None

    # Create output directory if needed
    output_path = Path(output_json)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    # Save training info to JSON
    try:
        with open(output_json, "w") as f:
            compact_json_dump(training_info.to_dict(), f, indent=2)
        if not quiet:
            print(f"  Saved training model to {output_json}")
    except Exception as e:
        if not quiet:
            print(f"  Error saving model: {e}")
        return None

    # Return the training info dict and genetic code
    return training_info.to_dict(), genetic_code


def load_existing_meta(meta_json_path, quiet=False):
    """Load existing meta.json file"""
    existing_models = []
    max_index = 0

    if meta_json_path.exists():
        try:
            with open(meta_json_path, "r") as f:
                content = f.read().strip()
                if content:
                    existing_models = json.loads(content)
                    if not quiet:
                        print(
                            f"Loaded {len(existing_models)} existing model(s) from meta.json"
                        )
                else:
                    if not quiet:
                        print("meta.json is empty, starting fresh")
        except json.JSONDecodeError as e:
            if not quiet:
                print(f"Warning: meta.json is not valid JSON ({e}), starting fresh")
    else:
        if not quiet:
            print("meta.json does not exist, will create new file")

    # Get the highest existing index
    for model in existing_models:
        desc = model.get("description", "")
        if "|" in desc:
            try:
                index = int(desc.split("|")[0])
                max_index = max(max_index, index)
            except:
                pass

    return existing_models, max_index


def process_directory(
    tsv_file, genomes_dir, models_dir, meta_json_path, min_length, append_meta, quiet
):
    """Process all families from TSV file"""

    # Check if TSV file exists
    if not Path(tsv_file).exists():
        print(f"Error: TSV file '{tsv_file}' does not exist", file=sys.stderr)
        sys.exit(1)

    # Check if genomes directory exists
    if not Path(genomes_dir).exists():
        print(
            f"Error: Genomes directory '{genomes_dir}' does not exist", file=sys.stderr
        )
        sys.exit(1)

    # Read the TSV to get genetic codes for each family
    try:
        df = pl.read_csv(tsv_file, separator="\t")
    except Exception as e:
        print(f"Error reading TSV file: {e}", file=sys.stderr)
        sys.exit(1)

    # Get unique family-genetic_code combinations
    family_gc = df.select(["family", "genetic_code"]).unique()

    # Also include models for genera
    genus_gc = (
        df.filter(pl.col("rank") == "genus")
        .select(["scientific_name", "genetic_code"])
        .unique()
        .rename({"scientific_name": "family"})
    )

    family_gc = pl.concat([family_gc, genus_gc])

    if not quiet:
        print(f"Found {len(family_gc)} unique family-genetic code combination(s)\n")

    # Create output directory for models
    models_path = Path(models_dir)
    models_path.mkdir(parents=True, exist_ok=True)

    genomes_path = Path(genomes_dir)
    meta_path = Path(meta_json_path)

    # Load existing meta.json
    if append_meta:
        existing_models, max_index = load_existing_meta(meta_path, quiet)
    else:
        existing_models = []
        max_index = 0
        if not quiet:
            print("Creating new meta.json (not appending to existing)")

    current_index = max_index + 1
    if not quiet and append_meta and max_index > 0:
        print(f"Starting new model index at: {current_index}\n")

    # Collect all new models
    new_models = []
    success_count = 0

    # Train model for each family
    for row in family_gc.iter_rows(named=True):
        family = row["family"]
        genetic_code = row["genetic_code"]

        fasta_file = genomes_path / f"{family}_{genetic_code}.fasta"
        output_json = models_path / f"{family}_{genetic_code}_model.json"

        result = train_pyrodigal_model(
            fasta_file, genetic_code, output_json, min_length, quiet
        )

        if result is not None:
            training_info_dict, gc_table = result

            # Get GC content and convert to percentage (rounded to 1 decimal)
            gc_content = training_info_dict.get("gc", 0)
            gc_percent = round(gc_content * 100, 1)

            # Get uses_sd from training_info
            uses_sd = training_info_dict.get("uses_sd", False)

            # Create description in format: index|model_name|V|gc|translation_table|uses_sd
            description = f"{current_index}|{family}_{genetic_code}_model|V|{gc_percent}|{gc_table}|{uses_sd}"

            # Add to new models list
            model_entry = {
                "description": description,
                "training_info": training_info_dict,
            }
            new_models.append(model_entry)

            if not quiet:
                print(f"  Model description: {description}")

            current_index += 1
            success_count += 1

        if not quiet:
            print()

    if not quiet:
        print(f"Successfully trained {success_count}/{len(family_gc)} model(s)")

    # Combine existing and new models
    if append_meta:
        all_models = existing_models + new_models
        if not quiet:
            print(
                f"Appending {len(new_models)} new model(s) to existing {len(existing_models)}"
            )
    else:
        all_models = new_models

    # Create meta.json parent directory if needed
    meta_path.parent.mkdir(parents=True, exist_ok=True)

    # Write to meta.json
    try:
        with open(meta_path, "w") as f:
            compact_json_dump(all_models, f, indent=2)
        if not quiet:
            print(
                f"Wrote meta.json with {len(all_models)} total model(s) to {meta_path}"
            )
    except Exception as e:
        print(f"Error writing meta.json: {e}", file=sys.stderr)
        sys.exit(1)


def process_single_file(
    fasta_file,
    genetic_code,
    output_json,
    model_name,
    meta_json_path,
    min_length,
    append_meta,
    quiet,
):
    """Process a single FASTA file"""

    if not Path(fasta_file).exists():
        print(f"Error: Input file '{fasta_file}' does not exist", file=sys.stderr)
        sys.exit(1)

    result = train_pyrodigal_model(
        fasta_file, genetic_code, output_json, min_length, quiet
    )

    if result is None:
        if not quiet:
            print("\nTraining failed.")
        sys.exit(1)

    training_info_dict, gc_table = result

    # Load existing meta.json if appending
    meta_path = Path(meta_json_path)
    if append_meta:
        existing_models, max_index = load_existing_meta(meta_path, quiet)
    else:
        existing_models = []
        max_index = 0

    current_index = max_index + 1

    # Get GC content
    gc_content = training_info_dict.get("gc", 0)
    gc_percent = round(gc_content * 100, 1)
    uses_sd = training_info_dict.get("uses_sd", False)

    # Use provided model name or derive from file
    if not model_name:
        model_name = Path(fasta_file).stem + "_model"

    description = f"{current_index}|{model_name}|V|{gc_percent}|{gc_table}|{uses_sd}"

    model_entry = {"description": description, "training_info": training_info_dict}

    if not quiet:
        print(f"\nModel description: {description}")

    # Update meta.json if specified
    if meta_json_path:
        if append_meta:
            all_models = existing_models + [model_entry]
        else:
            all_models = [model_entry]

        meta_path.parent.mkdir(parents=True, exist_ok=True)

        try:
            with open(meta_path, "w") as f:
                compact_json_dump(all_models, f, indent=2)
            if not quiet:
                print(f"Updated meta.json: {meta_path}")
        except Exception as e:
            print(f"Error writing meta.json: {e}", file=sys.stderr)
            sys.exit(1)

    if not quiet:
        print("\nTraining complete.")


def main():
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )

    # Mode selection
    mode_group = parser.add_mutually_exclusive_group()
    mode_group.add_argument(
        "--single", action="store_true", help="Process a single FASTA file"
    )

    # Single file mode arguments
    parser.add_argument("-i", "--input", help="Input FASTA file (single file mode)")

    parser.add_argument(
        "-g",
        "--genetic-code",
        type=int,
        help="Genetic code table number (single file mode, required with --single)",
    )

    parser.add_argument(
        "-o", "--output", help="Output JSON model file (single file mode)"
    )

    parser.add_argument(
        "-n",
        "--model-name",
        help="Model name for meta.json (single file mode, default: derived from filename)",
    )

    # Directory mode arguments
    parser.add_argument(
        "--tsv",
        default="riboviria_genetic_codes.tsv",
        help="Input TSV file with family-genetic code mappings (default: riboviria_genetic_codes.tsv)",
    )

    parser.add_argument(
        "--genomes-dir",
        default="genomes_clustered",
        help="Directory containing genome FASTA files (default: genomes_clustered)",
    )

    parser.add_argument(
        "--models-dir",
        default="pyrodigal_models",
        help="Output directory for model JSON files (default: pyrodigal_models)",
    )

    # Common arguments
    parser.add_argument(
        "--meta-json", help="Path to meta.json file (required for meta.json updates)"
    )

    parser.add_argument(
        "--no-append",
        action="store_true",
        help="Create new meta.json instead of appending to existing",
    )

    parser.add_argument(
        "--min-length",
        type=int,
        default=5000,
        help="Minimum total sequence length required for training (default: 5000)",
    )

    parser.add_argument(
        "-q", "--quiet", action="store_true", help="Suppress progress messages"
    )

    args = parser.parse_args()

    # Validate arguments
    if args.single:
        if not args.input:
            print("Error: --input is required in single file mode", file=sys.stderr)
            sys.exit(1)
        if not args.genetic_code:
            print(
                "Error: --genetic-code is required in single file mode", file=sys.stderr
            )
            sys.exit(1)
        if not args.output:
            print("Error: --output is required in single file mode", file=sys.stderr)
            sys.exit(1)

        # Process single file
        process_single_file(
            args.input,
            args.genetic_code,
            args.output,
            args.model_name,
            args.meta_json,
            args.min_length,
            not args.no_append,
            args.quiet,
        )
    else:
        # Directory mode
        if not args.meta_json:
            print("Error: --meta-json is required in directory mode", file=sys.stderr)
            sys.exit(1)

        process_directory(
            args.tsv,
            args.genomes_dir,
            args.models_dir,
            args.meta_json,
            args.min_length,
            not args.no_append,
            args.quiet,
        )


if __name__ == "__main__":
    main()
