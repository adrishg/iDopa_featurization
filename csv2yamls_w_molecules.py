import csv
import yaml
import os
import argparse
import json

def load_molecules(molecules_arg: str) -> dict:
    """Load molecules dictionary from a JSON string or a JSON file path."""
    if os.path.isfile(molecules_arg):
        with open(molecules_arg, 'r', encoding='utf-8') as f:
            return json.load(f)
    try:
        return json.loads(molecules_arg)
    except json.JSONDecodeError:
        raise ValueError("--molecules must be valid JSON or a path to a JSON file.")


def generate_yaml_files(csv_file_path: str, output_dir: str, molecules_dict: dict) -> None:
    """
    Reads a CSV file containing protein tags and sequences, and generates YAML files for each entry,
    pairing the protein with ligands provided in the molecules dictionary.

    The generated YAML matches your desired schema and includes a properties block with affinity:

    properties:
        - affinity:
            binder: B
    """

    # Ensure the output directory exists
    os.makedirs(output_dir, exist_ok=True)
    print(f"Output directory '{output_dir}' ensured.")

    try:
        with open(csv_file_path, mode='r', newline='', encoding='utf-8') as csvfile:
            reader = csv.DictReader(csvfile)
            print(f"Processing CSV file: {csv_file_path}")

            for row_num, row in enumerate(reader):
                tag = row.get('Tag')
                sequence = row.get('sequence')

                if not tag or not sequence:
                    print(f"Skipping row {row_num + 2} due to missing 'Tag' or 'sequence'.")
                    continue

                # Sanitize the tag for use in filenames: replace spaces with underscores
                sanitized_tag = tag.replace(' ', '_')

                for ligand_code, ligand_smiles in molecules_dict.items():
                    # Build YAML structure
                    yaml_data = {
                        "sequences": [
                            {
                                "protein": {
                                    "id": "A",
                                    "sequence": sequence
                                }
                            },
                            {
                                "ligand": {
                                    "id": "B",
                                    "smiles": ligand_smiles
                                }
                            }
                        ],
                        "properties": [
                            {
                                "affinity": {
                                    "binder": "B"
                                }
                            }
                        ]
                    }

                    # Output path
                    output_filename = f"{sanitized_tag}_{ligand_code}.yaml"
                    output_file_path = os.path.join(output_dir, output_filename)

                    # Write the YAML data to the file
                    try:
                        with open(output_file_path, 'w', encoding='utf-8') as outfile:
                            yaml.safe_dump(yaml_data, outfile, indent=2, sort_keys=False)
                        print(f"Generated: {output_file_path}")
                    except IOError as e:
                        print(f"Error writing file {output_file_path}: {e}")

    except FileNotFoundError:
        print(f"Error: CSV file not found at {csv_file_path}")
    except Exception as e:
        print(f"An unexpected error occurred: {e}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description=(
            "Generate YAML files from a CSV of protein sequences, using a molecules dictionary "
            "(JSON string or path), and append a properties.affinity block with binder B."
        )
    )
    parser.add_argument("csv_file", type=str, help="Path to the input CSV file.")
    parser.add_argument(
        "--output-dir",
        type=str,
        default="output_yamls",
        help="Directory to save the generated YAML files. Creates it if it doesn't exist. (default: output_yamls)"
    )
    parser.add_argument(
        "--molecules",
        type=str,
        required=True,
        help=(
            "JSON string or path to JSON file defining molecule dictionary, e.g. "
            "'{\"DOP\": \"C1=...\", \"5HT\": \"C1=...\"}' or molecules.json"
        )
    )

    args = parser.parse_args()

    try:
        molecules_dict = load_molecules(args.molecules)
    except ValueError as ve:
        print(f"Error: {ve}")
        raise SystemExit(1)

    generate_yaml_files(args.csv_file, args.output_dir, molecules_dict)
