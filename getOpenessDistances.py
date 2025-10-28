# -*- coding: utf-8 -*-
import os
import csv
import numpy as np

def extract_ca_coordinates(pdb_file, res1, res2, chain_id):
    ca_coords = {}
    with open(pdb_file, 'r') as file:
        for line in file:
            if line.startswith("ATOM") and line[12:16].strip() == "CA" and line[21] == chain_id:
                res_num = int(line[22:26].strip())
                if res_num in (res1, res2):
                    x = float(line[30:38])
                    y = float(line[38:46])
                    z = float(line[46:54])
                    ca_coords[res_num] = np.array([x, y, z])
    if res1 in ca_coords and res2 in ca_coords:
        return np.linalg.norm(ca_coords[res1] - ca_coords[res2])
    return None

def analyze_openess(parent_folder, res1, res2, chain_id, output_csv="openess_summary.csv"):
    results = []

    for tag in os.listdir(parent_folder):
        subfolder = os.path.join(parent_folder, tag)
        if not os.path.isdir(subfolder):
            continue

        distances = []
        for file in os.listdir(subfolder):
            if file.endswith(".pdb"):
                pdb_path = os.path.join(subfolder, file)
                distance = extract_ca_coordinates(pdb_path, res1, res2, chain_id)
                if distance is not None:
                    distances.append(distance)

        if distances:
            distances = np.array(distances)
            openess_avg = distances.mean()
            openess_min = distances.min()
            openess_max = distances.max()
            openess_range = openess_max - openess_min
            results.append([tag, openess_avg, openess_min, openess_max, openess_range])

    # Write CSV
    with open(output_csv, 'w', newline='') as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(["Tag", "openess_avg", "openess_min", "openess_max", "openess_range"])
        writer.writerows(results)

    print(f"Done. Output written to {output_csv}")

# Example usage:
# Replace with actual residue numbers and chain
if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description="Calculate openess metrics from PDB folders.")
    parser.add_argument("--parent-folder", required=True, help="Parent folder containing subfolders for each Tag")
    parser.add_argument("--res1", type=int, required=True, help="First residue number")
    parser.add_argument("--res2", type=int, required=True, help="Second residue number")
    parser.add_argument("--chain", type=str, required=True, help="Chain ID")
    parser.add_argument("--output-csv", default="openess_summary.csv", help="Output CSV filename")
    args = parser.parse_args()

    analyze_openess(args.parent_folder, args.res1, args.res2, args.chain, args.output_csv)

