# -*- coding: utf-8 -*-
import os
import json
import numpy as np
import matplotlib.pyplot as plt
from scipy.spatial.distance import pdist, squareform
import csv
import argparse
import sys
sys.stdout.reconfigure(encoding='utf-8')

def distance_map(pdb_file):
    coords = []
    plddt_scores = []
    pdb_code = os.path.basename(pdb_file).split('.')[0]

    with open(pdb_file, 'r') as f:
        for line in f:
            if line.startswith('ATOM') and line[12:16].strip() == 'CA' and line[21] == 'A':
                x = float(line[30:38])
                y = float(line[38:46])
                z = float(line[46:54])
                bfactor = float(line[60:66])
                coords.append([x, y, z])
                plddt_scores.append(bfactor / 100.0)

    if not coords:
        raise ValueError(f"No CA atoms found in {pdb_file}")

    coords = np.array(coords)
    plddt_scores = np.array(plddt_scores)
    dist_matrix = squareform(pdist(coords, 'euclidean'))
    plddt_weight_matrix = np.outer(plddt_scores, plddt_scores)
    return dist_matrix, plddt_weight_matrix, pdb_code

def load_pae_matrix(npz_file, key):
    npz = np.load(npz_file)
    if key not in npz:
        raise ValueError(f"Key '{key}' not found in {npz_file}")
    matrix = npz[key]
    scaling_factor = 5.0
    return np.exp(-matrix / scaling_factor), matrix

def process_all_folders(parent_folder, output_folder):
    os.makedirs(output_folder, exist_ok=True)
    subfolders = [f.path for f in os.scandir(parent_folder) if f.is_dir()]
    composite_variances = []

    for folder in subfolders:
        folder_name = os.path.basename(folder)
        pdb_files = [os.path.join(folder, f) for f in os.listdir(folder) if f.endswith('.pdb')]
        matrices_unweighted, matrices_plddt_weighted = [], []
        matrices_pae_weighted, matrices_pde_weighted = [], []
        pae_raw_values, pde_raw_values = [], []
        complex_pde_values = []
        pae_min_values, pae_max_values = [], []
        pde_min_values, pde_max_values = [], []

        for pdb_file in pdb_files:
            try:
                out_dir = os.path.join(output_folder, folder_name)
                os.makedirs(out_dir, exist_ok=True)
                dist_matrix, plddt_weight_matrix, pdb_code = distance_map(pdb_file)
                basename = os.path.basename(pdb_file).replace('.pdb', '')
                json_file, pae_file, pde_file = None, None, None

                for f in os.listdir(folder):
                    if f.startswith('pae_' + basename) and f.endswith('.npz'):
                        pae_file = os.path.join(folder, f)
                    elif f.startswith('pde_' + basename) and f.endswith('.npz'):
                        pde_file = os.path.join(folder, f)
                    elif f.startswith('confidence_' + basename) and f.endswith('.json'):
                        json_file = os.path.join(folder, f)

                if pae_file:
                    try:
                        pae_weight, pae_raw = load_pae_matrix(pae_file, 'pae')
                        pae_raw_values.append(pae_raw)
                        if pae_weight.shape == dist_matrix.shape:
                            matrices_pae_weighted.append(dist_matrix * pae_weight)
                            pae_min_values.append(np.min(pae_raw))
                            pae_max_values.append(np.max(pae_raw))
                        else:
                            print(f"Skipping PAE weighting for {basename}: shape mismatch {pae_weight.shape} vs {dist_matrix.shape}")
                    except Exception as e:
                        print(f"Warning (PAE): {e}")

                if pde_file:
                    try:
                        pde_weight, pde_raw = load_pae_matrix(pde_file, 'pde')
                        pde_raw_values.append(pde_raw)
                        if pde_weight.shape == dist_matrix.shape:
                            matrices_pde_weighted.append(dist_matrix * pde_weight)
                            pde_min_values.append(np.min(pde_raw))
                            pde_max_values.append(np.max(pde_raw))
                        else:
                            print(f"Skipping PDE weighting for {basename}: shape mismatch {pde_weight.shape} vs {dist_matrix.shape}")
                    except Exception as e:
                        print(f"Warning (PDE): {e}")

                if json_file:
                    try:
                        with open(json_file, 'r') as jf:
                            conf_data = json.load(jf)
                        if 'complex_pde' in conf_data:
                            complex_pde_values.append(conf_data['complex_pde'])
                    except Exception as e:
                        print(f"Warning reading complex_pde from {json_file}: {e}")

                weighted_plddt = dist_matrix * plddt_weight_matrix
                matrices_unweighted.append(dist_matrix)
                matrices_plddt_weighted.append(weighted_plddt)

            except Exception as e:
                print(f"Failed: {pdb_file} - {e}")

        if not matrices_unweighted:
            print(f"Skipping {folder_name}, no valid PDBs.")
            continue

        shapes = {m.shape for m in matrices_unweighted}
        if len(shapes) > 1:
            print(f"Mixed dimensions in {folder_name}, skipping variance calculation.")
            continue

        def compute_variance(matrices):
            stacked = np.stack(matrices)
            return np.var(stacked, axis=0), np.mean(np.var(stacked, axis=0))

        var_unweighted, compvar_unweighted = compute_variance(matrices_unweighted)
        var_plddt, compvar_plddt = compute_variance(matrices_plddt_weighted)
        var_pae, compvar_pae = compute_variance(matrices_pae_weighted) if matrices_pae_weighted else (None, 'NA')
        var_pde, compvar_pde = compute_variance(matrices_pde_weighted) if matrices_pde_weighted else (None, 'NA')

        if complex_pde_values:
            mean_complex_pde = np.mean(complex_pde_values)
            var_complex_pde = np.var(complex_pde_values) * 1000
            min_complex_pde = np.min(complex_pde_values)
            max_complex_pde = np.max(complex_pde_values)
        else:
            mean_complex_pde = var_complex_pde = min_complex_pde = max_complex_pde = 'NA'

        pae_avg = f"{np.mean([v for mat in pae_raw_values for v in mat.flatten()]):.3f}" if pae_raw_values else 'NA'
        pae_min = f"{np.min(pae_min_values):.3f}" if pae_min_values else 'NA'
        pae_max = f"{np.max(pae_max_values):.3f}" if pae_max_values else 'NA'

        pde_avg = f"{np.mean([v for mat in pde_raw_values for v in mat.flatten()]):.3f}" if pde_raw_values else 'NA'
        pde_min = f"{np.min(pde_min_values):.3f}" if pde_min_values else 'NA'
        pde_max = f"{np.max(pde_max_values):.3f}" if pde_max_values else 'NA'

        composite_variances.append([
            folder_name,
            f"{compvar_unweighted:.3f}",
            f"{compvar_plddt:.3f}",
            f"{compvar_pae:.3f}" if compvar_pae != 'NA' else 'NA',
            f"{compvar_pde:.3f}" if compvar_pde != 'NA' else 'NA',
            f"{mean_complex_pde:.3f}" if mean_complex_pde != 'NA' else 'NA',
            f"{var_complex_pde:.3f}" if var_complex_pde != 'NA' else 'NA',
            f"{min_complex_pde:.3f}" if min_complex_pde != 'NA' else 'NA',
            f"{max_complex_pde:.3f}" if max_complex_pde != 'NA' else 'NA',
            pae_min, pae_max, pae_avg,
            pde_min, pde_max, pde_avg
        ])

        print(f"Done with {folder_name}")

    csv_path = os.path.join(output_folder, "composite_variances.csv")
    with open(csv_path, 'w', newline='') as csvfile:
        csv_writer = csv.writer(csvfile)
        csv_writer.writerow([
            'Tag', 'variance_avg', 'variance_pLDDT_w', 'variance_PAE_w', 'variance_PDE_w',
            'complex_PDE_avg', 'complex_PDE_var', 'complex_PDE_min', 'complex_PDE_max',
            'PAE_min', 'PAE_max', 'PAE_avg', 'PDE_min', 'PDE_max', 'PDE_avg'])
        csv_writer.writerows(composite_variances)
    print(f"\nComposite variances saved to: {csv_path}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Compute composite variance and complex_pde statistics from AlphaFold models.")
    parser.add_argument("--input_dir", required=True, help="Path to input parent folder (folder of folders)")
    parser.add_argument("--output_dir", required=True, help="Path to output folder")
    args = parser.parse_args()
    process_all_folders(args.input_dir, args.output_dir)

