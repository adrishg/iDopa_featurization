# -*- coding: utf-8 -*-
import os
import json
import csv
import argparse

def extract_affinity_values(input_dir, output_csv):
    rows = []
    for folder_name in os.listdir(input_dir):
        folder_path = os.path.join(input_dir, folder_name)
        if not os.path.isdir(folder_path):
            continue

        # Search for JSON file with pattern: affinity_*nameOfFolder*.json
        expected_prefix = f"affinity_{folder_name}"
        for file in os.listdir(folder_path):
            if file.startswith(expected_prefix) and file.endswith(".json"):
                json_path = os.path.join(folder_path, file)
                try:
                    with open(json_path, 'r') as f:
                        data = json.load(f)
                        row = {
                            "Tag": folder_name,
                            "affinity_pred_value": data.get("affinity_pred_value", ""),
                            "affinity_probability_binary": data.get("affinity_probability_binary", "")
                        }
                        rows.append(row)
                except Exception as e:
                    print(f"[WARNING] Failed to read {json_path}: {e}")
                break  # stop after finding first match

    # Write to CSV
    with open(output_csv, 'w', newline='') as csvfile:
        fieldnames = ["Tag", "affinity_pred_value", "affinity_probability_binary"]
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            writer.writerow(row)

    print(f"[INFO] CSV written to {output_csv}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Extract affinity values from JSON files in subfolders.")
    parser.add_argument('--input-dir', required=True, help="Path to the parent folder containing subfolders.")
    parser.add_argument('--output-csv', required=True, help="Output CSV file path.")

    args = parser.parse_args()
    extract_affinity_values(args.input_dir, args.output_csv)

