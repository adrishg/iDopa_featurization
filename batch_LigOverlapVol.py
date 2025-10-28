# -*- coding: utf-8 -*-
import os
import sys
import argparse
from Bio.PDB import PDBParser, Superimposer
import numpy as np
from collections import defaultdict
import csv

BINDING_POCKET_RESIDUES = {
    "A": [12, 65, 67, 80, 355]
}
LIGAND_RESIDUE_NAME = "LIG"

def default_vdw_radius_factory():
    return 1.5

VAN_DER_WAALS_RADII = defaultdict(default_vdw_radius_factory)
VAN_DER_WAALS_RADII.update({"C": 1.70, "O": 1.52, "N": 1.55, "S": 1.80, "H": 1.20})

def get_atoms_from_selection(structure, selection_dict):
    selected_atoms = []
    for model in structure:
        for chain_id, res_nums in selection_dict.items():
            if chain_id in model:
                chain = model[chain_id]
                for res_num in res_nums:
                    residue_id = (' ', res_num, ' ')
                    if residue_id in chain:
                        selected_atoms.extend(chain[residue_id].get_atoms())
    print(f"      Found {len(selected_atoms)} binding pocket atoms for selection.")
    return selected_atoms

def get_ligand_atoms(structure, ligand_name):
    ligand_atoms = []
    for model in structure:
        for chain in model:
            for residue in chain:
                if residue.get_resname().strip() == ligand_name.strip() and residue.id[0] != ' ':
                    ligand_atoms.extend(residue.get_atoms())
    print(f"      Found {len(ligand_atoms)} ligand atoms for '{ligand_name}'.")
    return ligand_atoms

def calculate_ligand_volume_monte_carlo(atoms, vdw_radii_dict, n_points=100000):
    if not atoms:
        return 0.0
    print(f"        Starting Monte Carlo volume estimation with {n_points} points for {len(atoms)} atoms...")
    coords = np.array([atom.get_coord() for atom in atoms])
    radii = np.array([
        vdw_radii_dict.get(atom.element.strip() if atom.element else atom.get_name()[0], 1.5)
        for atom in atoms
    ])
    min_coords = np.min(coords - radii[:, np.newaxis], axis=0)
    max_coords = np.max(coords + radii[:, np.newaxis], axis=0)
    box_min = min_coords - 0.5
    box_max = max_coords + 0.5
    box_dimensions = box_max - box_min
    box_volume = np.prod(box_dimensions)
    random_points = box_min + box_dimensions * np.random.rand(n_points, 3)
    diffs = random_points[:, np.newaxis, :] - coords[np.newaxis, :, :]
    dists_squared = np.sum(diffs ** 2, axis=2)
    within_any_sphere = np.any(dists_squared <= radii[np.newaxis, :] ** 2, axis=1)
    points_in_molecule = np.sum(within_any_sphere)
    estimated_volume = (points_in_molecule / n_points) * box_volume
    print(f"        Monte Carlo estimation complete. Estimated volume: {estimated_volume:.2f} A^3")
    return estimated_volume

def get_average_plddt(atoms):
    if not atoms:
        return 0.0
    return np.mean([atom.get_bfactor() / 100.0 for atom in atoms])

def process_pdb_files_in_subfolder(subfolder_path, binding_pocket_residues, ligand_name, vdw_radii):
    print(f"\n--- Processing Subfolder: {os.path.basename(subfolder_path)} ---")
    parser = PDBParser()
    pdb_files = [f for f in os.listdir(subfolder_path) if f.endswith(".pdb")]
    subfolder_name = os.path.basename(subfolder_path)

    all_ligand_atoms_aligned = []
    individual_results = []
    reference_structure = None
    ref_atoms_for_superimposition = []

    for pdb_file in pdb_files:
        pdb_path = os.path.join(subfolder_path, pdb_file)
        print(f"  Checking {pdb_file} for reference structure...")
        structure = parser.get_structure(pdb_file, pdb_path)
        binding_pocket_atoms = get_atoms_from_selection(structure, binding_pocket_residues)
        if binding_pocket_atoms:
            reference_structure = structure
            ref_atoms_for_superimposition = binding_pocket_atoms
            print(f"  Reference structure set: {pdb_file}")
            break

    if not reference_structure:
        print(f"  No suitable reference structure found in {subfolder_name}. Skipping.")
        return [], 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, subfolder_name

    for pdb_file in pdb_files:
        pdb_path = os.path.join(subfolder_path, pdb_file)
        print(f"  Processing {pdb_file}...")
        structure = parser.get_structure(pdb_file, pdb_path)
        current_binding_pocket_atoms = get_atoms_from_selection(structure, binding_pocket_residues)
        ligand_atoms = get_ligand_atoms(structure, ligand_name)
        if not ligand_atoms:
            print(f"    No ligand found. Skipping.")
            continue

        if structure.id != reference_structure.id:
            if len(ref_atoms_for_superimposition) == len(current_binding_pocket_atoms):
                print(f"    Superimposing onto reference...")
                sup = Superimposer()
                ref_atoms_sorted = sorted(ref_atoms_for_superimposition, key=lambda a: a.get_name())
                current_atoms_sorted = sorted(current_binding_pocket_atoms, key=lambda a: a.get_name())
                sup.set_atoms(ref_atoms_sorted, current_atoms_sorted)
                sup.apply(structure.get_atoms())
                print(f"    RMSD: {sup.rms:.3f} Å")
                ligand_atoms = get_ligand_atoms(structure, ligand_name)

        print(f"    Calculating volume...")
        unweighted_vol = calculate_ligand_volume_monte_carlo(ligand_atoms, vdw_radii)
        avg_plddt = get_average_plddt(ligand_atoms)
        weighted_vol = unweighted_vol * avg_plddt
        print(f"    Volume: {unweighted_vol:.2f} Å^3 | pLDDT avg: {avg_plddt:.2f} | Weighted+: {weighted_vol:.2f}")

        individual_results.append({
            'PDB_File': pdb_file,
            'Individual_Unweighted_Ligand_Volume_A^3': unweighted_vol,
            'Individual_Ligand_Avg_pLDDT': avg_plddt,
            'Individual_Weighted_Ligand_Volume_A^3': weighted_vol
        })
        all_ligand_atoms_aligned.extend(ligand_atoms)

    if not all_ligand_atoms_aligned:
        print(f"  No aligned ligand atoms for combined volume calculation in {subfolder_name}.")
        return individual_results, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, subfolder_name

    print(f"  Calculating combined volume for {subfolder_name}...")
    combined_unweighted_volume = calculate_ligand_volume_monte_carlo(all_ligand_atoms_aligned, vdw_radii, 500000)
    plddt_vals = [atom.get_bfactor() / 100.0 for atom in all_ligand_atoms_aligned]
    combined_avg_plddt = np.mean(plddt_vals)
    combined_min_plddt = np.min(plddt_vals)
    combined_max_plddt = np.max(plddt_vals)
    combined_weighted_vol_pos = combined_unweighted_volume * combined_avg_plddt
    combined_weighted_vol_neg = combined_unweighted_volume * (1.0 - combined_avg_plddt)
    print(f"    Combined Volume: {combined_unweighted_volume:.2f} Å^3")
    print(f"    Weighted+: {combined_weighted_vol_pos:.2f} | Weighted-: {combined_weighted_vol_neg:.2f}")
    print(f"    pLDDT avg: {combined_avg_plddt:.2f} | min: {combined_min_plddt:.2f} | max: {combined_max_plddt:.2f}")

    return individual_results, combined_weighted_vol_pos, combined_weighted_vol_neg, combined_unweighted_volume, combined_avg_plddt, combined_min_plddt, combined_max_plddt, subfolder_name

def write_overall_summary_to_csv(summary_list, output_filepath):
    fieldnames = [
        'Tag', 'Folder_Path', 'overlap_volume', 'overlap_w_pos_volume',
        'overlap_w_neg_volume', 'ligand_pLDDT_avg', 'ligand_pLDDT_min', 'ligand_pLDDT_max'
    ]
    with open(output_filepath, 'w', newline='') as csvfile:
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
        writer.writeheader()
        for row in summary_list:
            writer.writerow(row)

def write_individual_results_to_csv(results_list, output_filepath):
    fieldnames = [
        'PDB_File', 'Individual_Unweighted_Ligand_Volume_A^3',
        'Individual_Ligand_Avg_pLDDT', 'Individual_Weighted_Ligand_Volume_A^3'
    ]
    with open(output_filepath, 'w', newline='') as csvfile:
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
        writer.writeheader()
        for row in results_list:
            writer.writerow(row)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Ligand Volume Analysis")
    parser.add_argument("--input_dir", required=True, help="Path to parent folder containing subfolders with PDBs")
    parser.add_argument("--output_dir", required=True, help="Directory where summary CSV and images will be saved")
    args = parser.parse_args()

    parent_folder_path = args.input_dir
    output_dir = args.output_dir
    os.makedirs(output_dir, exist_ok=True)

    overall_summary_csv_path = os.path.join(output_dir, "overall_folder_summary.csv")
    all_subfolder_summary_results = []

    for item_name in os.listdir(parent_folder_path):
        current_subfolder_path = os.path.join(parent_folder_path, item_name)
        if not os.path.isdir(current_subfolder_path):
            continue
        result = process_pdb_files_in_subfolder(
            current_subfolder_path,
            BINDING_POCKET_RESIDUES,
            LIGAND_RESIDUE_NAME,
            VAN_DER_WAALS_RADII
        )
        (individual_results, pos_vol, neg_vol, raw_vol,
         avg_plddt, min_plddt, max_plddt, subfolder_name) = result

        individual_csv_path = os.path.join(output_dir, f"{subfolder_name}_individual_ligand_analysis.csv")
        write_individual_results_to_csv(individual_results, individual_csv_path)

        summary = {
            'Tag': subfolder_name,
            'Folder_Path': current_subfolder_path,
            'overlap_volume': raw_vol,
            'overlap_w_pos_volume': pos_vol,
            'overlap_w_neg_volume': neg_vol,
            'ligand_pLDDT_avg': avg_plddt,
            'ligand_pLDDT_min': min_plddt,
            'ligand_pLDDT_max': max_plddt
        }
        all_subfolder_summary_results.append(summary)

    write_overall_summary_to_csv(all_subfolder_summary_results, overall_summary_csv_path)
    print("\nAnalysis complete. Results written to:")
    print(f"  Summary CSV: {overall_summary_csv_path}")
    for summary in all_subfolder_summary_results:
        print(f"  Tag: {summary['Tag']} → {summary['overlap_volume']:.2f} Å^3 total")

