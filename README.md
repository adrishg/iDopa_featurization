# iDopa_featurization

Scripts for automated generation, modeling, and structural featurization of iDopaSnFR biosensor variants using [Boltz2](https://github.com/jwohlwend/boltz), an all-atom diffusion model for joint protein–ligand structure prediction.

Starting from a CSV of variant sequences, the pipeline produces a unified feature table describing each biosensor variant and its predicted interactions with small-molecule ligands.

---

## Repository structure

```
pipeline.sh                  — Full pipeline (YAML generation → Boltz2 → analysis), SLURM-ready
analysis.sh                  — Feature extraction + CSV merging for one prediction set

csv2yamls_w_molecules.py     — Generate Boltz2 YAML inputs from a CSV of sequences
batch_LigOverlapVol.py       — Compute ligand overlap volume across ensemble models
batch_distanceMaps_variance.py — Compute distance-map variance and confidence metrics
getAffinities.py             — Extract affinity predictions from Boltz2 output files
getOpenessDistances.py       — Measure PBP domain openness (Cα–Cα distance)
getOpenessDistancesProp.py   — Openness with open/closed proportion statistics
merge_csv_tags.py            — Merge feature CSVs by Tag column
```

---

## Input

A CSV file with at minimum two columns:

| Tag     | sequence         |
|---------|------------------|
| variantA | MKTLLLT...      |
| variantB | MKTLLLT...      |

---

## Pipeline overview

### 1. Generate Boltz2 input YAMLs — `csv2yamls_w_molecules.py`

Creates one YAML per variant per ligand, with an `affinity` property block.

```bash
python3 csv2yamls_w_molecules.py data.csv \
    --output-dir output_yamls/ \
    --molecules '{"DOP": "C1=CC(=C(C=C1CCN)O)O", "5HT": "C1=CC2=C(C=C1O)C(=CN2)CCN"}'
```

`--molecules` accepts a JSON string or a path to a `.json` file.

### 2. Run Boltz2 predictions

The `pipeline.sh` calls `boltz predict` directly. The conda environment must have `boltz` installed. Adapt the `eval "$(conda shell.bash hook)"` / `conda activate` lines to your environment.

### 3. Extract features — `analysis.sh`

Runs all Python analysis scripts on a predictions directory and merges the results.

```bash
bash analysis.sh /path/to/predictions /path/to/output
```

Produces `volumes_variances_affinities_openess_clean.csv` with columns:

| Column | Source |
|--------|--------|
| `Tag` | variant name (suffix stripped) |
| `Ligand` | DOP or 5HT |
| `overlap_volume`, `overlap_w_pos_volume`, `overlap_w_neg_volume` | `batch_LigOverlapVol.py` |
| `ligand_pLDDT_avg/min/max` | `batch_LigOverlapVol.py` |
| `variance_avg`, `variance_pLDDT_w` | `batch_distanceMaps_variance.py` |
| `complex_PDE_avg/var/min/max`, `PAE_avg`, `PDE_avg` | `batch_distanceMaps_variance.py` |
| `affinity_pred_value`, `affinity_probability_binary` | `getAffinities.py` |
| `openess_avg/min/max/range` | `getOpenessDistances.py` |

### 4. Full pipeline — `pipeline.sh`

Orchestrates steps 1–3 for both DOP and 5HT separately via SLURM.
Edit `PROJECT_PATH` and `DATA_PATH` at the top of the file before submitting.

```bash
sbatch pipeline.sh
```

---

## Running individual scripts in Colab

Each Python script is self-contained and can be imported or called directly:

```python
# In Colab — clone the repo and call scripts from the repo directory
import subprocess, sys

subprocess.run([sys.executable, "csv2yamls_w_molecules.py",
    "data.csv", "--output-dir", "yamls/",
    "--molecules", '{"DOP": "C1=CC(=C(C=C1CCN)O)O"}'], check=True)
```

Or import functions directly:

```python
from batch_LigOverlapVol import process_pdb_files_in_subfolder
from getAffinities import extract_affinity_values
from merge_csv_tags import merge_csv_files
```

---

## Dependencies

```
biopython
numpy
scipy
matplotlib
pandas
pyyaml
```
