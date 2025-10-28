import os
import csv
import numpy as np
import argparse

def extract_ca_coordinates(pdb_file, res1, res2, chain_id):
    """
    Extracts the coordinates of the C-alpha atoms for two specified residues
    in a given PDB file and calculates the Euclidean distance between them.

    Args:
        pdb_file (str): Path to the PDB file.
        res1 (int): The residue number of the first residue.
        res2 (int): The residue number of the second residue.
        chain_id (str): The chain ID.

    Returns:
        float or None: The distance between the C-alpha atoms in Angstroms,
                       or None if coordinates are not found.
    """
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

def analyze_openess(parent_folder, res1, res2, chain_id, open_threshold, output_csv="openess_summary.csv"):
    """
    Analyzes 'openess' metrics, including the proportion of open vs. closed
    models, for PDB files within a nested folder structure.

    Args:
        parent_folder (str): The parent directory containing subfolders (tags).
        res1 (int): The first residue number.
        res2 (int): The second residue number.
        chain_id (str): The chain ID.
        open_threshold (float): The distance threshold to define an "open" model.
        output_csv (str): The name of the output CSV file.
    """
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
            
            # Openness metrics
            openess_avg = distances.mean()
            openess_min = distances.min()
            openess_max = distances.max()
            openess_range = openess_max - openess_min
            
            # Proportional analysis
            total_models = len(distances)
            open_models = np.sum(distances > open_threshold)
            closed_models = total_models - open_models
            
            proportion_open = open_models / total_models if total_models > 0 else 0
            proportion_closed = closed_models / total_models if total_models > 0 else 0
            
            results.append([
                tag, 
                openess_avg, 
                openess_min, 
                openess_max, 
                openess_range,
                proportion_open,
                proportion_closed
            ])

    # Write CSV
    with open(output_csv, 'w', newline='') as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow([
            "Tag", 
            "openess_avg", 
            "openess_min", 
            "openess_max", 
            "openess_range",
            "proportion_open",
            "proportion_closed"
        ])
        writer.writerows(results)

    print(f"Done. Output written to {output_csv}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Calculate openess metrics from PDB folders.")
    parser.add_argument("--parent-folder", required=True, help="Parent folder containing subfolders for each Tag")
    parser.add_argument("--res1", type=int, required=True, help="First residue number")
    parser.add_argument("--res2", type=int, required=True, help="Second residue number")
    parser.add_argument("--chain", type=str, required=True, help="Chain ID")
    parser.add_argument("--open-threshold", type=float, default=17.5, help="Distance threshold (in Angstroms) to define an 'open' model. Default is 17.5.")
    parser.add_argument("--output-csv", default="openess_summary.csv", help="Output CSV filename")
    args = parser.parse_args()

    analyze_openess(args.parent_folder, args.res1, args.res2, args.chain, args.open_threshold, args.output_csv)
