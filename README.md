# iDopa_featurization

[![Open run notebook in Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/adrishg/iDopa_featurization/blob/main/runiDopaBoltz2.ipynb)

Tools for iDopaSnFR Structure Ensemble featurization and experimental-feature prediction. The repository combines Boltz2 structure generation, ensemble-derived structural feature extraction, and trained prediction bundles for biosensor variant evaluation.

Starting from a CSV of variant sequences or the Colab notebook input form, the workflow produces Boltz2 structure ensembles, unified feature tables, and prediction-ready outputs for dopamine, serotonin, custom ligand, and APO runs.

---

## Repository structure

```
pipeline.sh                  — Full pipeline (YAML generation → Boltz2 → analysis), SLURM-ready
analysis.sh                  — Feature extraction + CSV merging for one prediction set
runiDopaBoltz2.ipynb              — Colab notebook for ligand/APO Boltz2 runs and repo-backed featurization
runiDopaBoltz2_SEprediction.ipynb — Colab notebook for Boltz2 + Structure Ensemble feature prediction

featurization/
  csv2yamls_w_molecules.py       — Generate Boltz2 YAML inputs from a CSV of sequences
  batch_LigOverlapVol.py         — Compute ligand overlap volume across ensemble models
  batch_distanceMaps_variance.py — Compute distance-map variance and confidence metrics
  getAffinities.py               — Extract affinity predictions from Boltz2 output files
  getOpenessDistances.py         — Measure PBP domain openness (Cα–Cα distance)
  getOpenessDistancesProp.py     — Openness with open/closed proportion statistics
  merge_csv_tags.py              — Merge feature CSVs by Tag column

SEmodels/
  prediction_bundle_*.zip         — Structure Ensemble prediction model bundles
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

### 0. Run from Colab

Open the Boltz2 run notebook directly:

[![Open runiDopaBoltz2 in Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/adrishg/iDopa_featurization/blob/main/runiDopaBoltz2.ipynb)

Open the Structure Ensemble prediction notebook directly:

[![Open runiDopaBoltz2_SEprediction in Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/adrishg/iDopa_featurization/blob/main/runiDopaBoltz2_SEprediction.ipynb)

By default, the run notebook starts DOP and 5HT ligand runs. It also supports:

- DOP and 5HT ligand runs with Boltz2 affinity requests.
- Custom ligand runs using a label and SMILES.
- Optional APO protein-only runs with no ligand chain and no affinity request.
- Optional repo-backed feature extraction from `featurization/`.

### 1. Generate Boltz2 input YAMLs — `featurization/csv2yamls_w_molecules.py`

Creates one YAML per variant per ligand, with an `affinity` property block.

```bash
python3 featurization/csv2yamls_w_molecules.py data.csv \
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

For ligand predictions, produces `volumes_variances_affinities_openess_clean.csv` with columns:

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

subprocess.run([sys.executable, "featurization/csv2yamls_w_molecules.py",
    "data.csv", "--output-dir", "yamls/",
    "--molecules", '{"DOP": "C1=CC(=C(C=C1CCN)O)O"}'], check=True)
```

Or import functions directly:

```python
from featurization.batch_LigOverlapVol import process_pdb_files_in_subfolder
from featurization.getAffinities import extract_affinity_values
from featurization.merge_csv_tags import merge_csv_files
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
