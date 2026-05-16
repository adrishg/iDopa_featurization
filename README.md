# iDopaSnFR Structural Ensemble Modeling & Prediction

[![Run Boltz2 Notebook](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/adrishg/iDopa_featurization/blob/main/runiDopaBoltz2.ipynb)
[![Run SE Prediction Notebook](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/adrishg/iDopa_featurization/blob/main/runiDopaBoltz2_SEprediction.ipynb)

Structure-ensemble-based prediction framework for engineering and evaluating **iDopaSnFR dopamine biosensor variants** using **Boltz-2 structural modeling**, ensemble-derived structural features, and machine learning prediction models.

This repository combines:

- Boltz-2 structure generation
- Structural ensemble (SE) feature extraction
- Automated visualization and analysis
- Experimental-property prediction models
- Future integration with protein language model (PLM) voting systems

The framework was designed for iterative biosensor engineering, where ensembles of predicted structures are used to estimate experimentally relevant biosensor properties before experimental screening.

> Small context: iDopaSnFR is a dopamine biosensor derived from a PBP/Venus-Flytrap-like scaffold engineered for neurotransmitter sensing.

---

# Overview

The repository currently contains **two main Colab notebooks**:

| Notebook | Purpose |
|---|---|
| `runiDopaBoltz2.ipynb` | Generate and visualize Boltz-2 structural models |
| `runiDopaBoltz2_SEprediction.ipynb` | Run structural ensemble feature extraction + experimental-property prediction |

The workflow supports:

- Dopamine (`DOP`) ligand modeling
- Serotonin (`5HT`) ligand modeling
- Optional APO runs
- Optional custom ligands via SMILES
- Structural ensemble featurization
- Prediction of experimental biosensor properties

---

# Repository Structure

```text
runiDopaBoltz2.ipynb
    Main Boltz-2 notebook for generating and visualizing models

runiDopaBoltz2_SEprediction.ipynb
    Structural Ensemble prediction workflow notebook

pipeline.sh
    Full SLURM-ready pipeline

analysis.sh
    Feature extraction and CSV merging

featurization/
    csv2yamls_w_molecules.py
    batch_LigOverlapVol.py
    batch_distanceMaps_variance.py
    getAffinities.py
    getOpenessDistances.py
    getOpenessDistancesProp.py
    merge_csv_tags.py

SEmodels/
    prediction_bundle_*.zip
    Best-performing trained prediction models

PLM/
    Placeholder for future Protein Language Model integration
```

---

# Structural Ensemble (SE) Workflow

The Structural Ensemble workflow uses multiple Boltz-2 models per sequence to estimate conformational variability and extract predictive structural descriptors.

Current standard workflow:

```text
Sequence
    ↓
Boltz-2 Ensemble Modeling
    ↓
Structural Ensemble Feature Extraction
    ↓
Prediction Models
    ↓
Experimental Property Predictions
```

Typical run:

- ~25 Boltz-2 models per ligand condition
- Dopamine ensemble
- Serotonin ensemble
- Optional APO ensemble

The extracted ensemble features are then used as input for trained prediction models.

---

# What the Structural Ensemble Features Capture

The SE featurization pipeline extracts metrics describing:

## Ligand confidence and consistency

- Ligand pLDDT
- Ligand positional variance
- Ensemble ligand pose agreement
- Confidence-weighted ligand overlap

## Predicted affinity metrics

- Boltz-2 affinity prediction values
- Binary affinity confidence estimates

## Conformational variability

- Distance-map variance
- PDE / PAE ensemble metrics
- Ensemble spread

## Biosensor opening/closing behavior

- Open-state distances
- Closed-state distances
- Openness range across ensemble

## Ensemble stability

- Cross-model structural consistency
- Confidence-weighted variability estimates

---

# Main Prediction Targets

The repository currently predicts several experimentally relevant biosensor properties independently.

Current prediction targets include:

| Property | Description |
|---|---|
| `DA_EC50` | Dopamine EC50 |
| `DA_F_F` | Dopamine ΔF/F |
| `5HT_EC50` | Serotonin EC50 |
| `5HT_F_F` | Serotonin ΔF/F |
| `F0` | Baseline fluorescence |
| `5HT_Selectivity` | Dopamine vs serotonin selectivity |

Predictions are currently returned as:

- Predicted class
- Human-readable label
- Confidence score

Example output:

```text
DA_EC50__pred_label           Low (<= 231.3)
DA_EC50__pred_conf            0.946

DA_F_F__pred_label            Mid (1.894–2.537)
DA_F_F__pred_conf             0.570

5HT_Selectivity__pred_label   Mid (1.705–9.937)
5HT_Selectivity__pred_conf    0.995
```

---

# Notebook 1 — Boltz-2 Modeling

## `runiDopaBoltz2.ipynb`

[![Run Notebook](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/adrishg/iDopa_featurization/blob/main/runiDopaBoltz2.ipynb)

This notebook is intended for:

- Running Boltz-2 predictions
- Visualizing generated structures
- Comparing ligand-bound states
- Rapid exploratory modeling

Supports:

- Dopamine (`DOP`)
- Serotonin (`5HT`)
- Optional APO runs
- Optional custom ligands via SMILES

Features:

- Interactive visualization
- Ensemble visualization
- Reference structure overlay
- Ligand comparison
- Quick structural inspection

---

# Notebook 2 — Structural Ensemble Prediction

## `runiDopaBoltz2_SEprediction.ipynb`

[![Run Notebook](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/adrishg/iDopa_featurization/blob/main/runiDopaBoltz2_SEprediction.ipynb)

This notebook performs the full SE workflow:

1. Generate Boltz-2 ensembles
2. Extract structural ensemble features
3. Visualize ensemble metrics
4. Compare against references
5. Run trained prediction models
6. Generate prediction tables

This notebook:

- Runs ~25 models per ligand condition
- Extracts all structural ensemble descriptors
- Automatically assembles feature tables
- Applies latest trained prediction bundles

---

# Example SE Feature Table

The generated feature tables include metrics such as:

| Feature Type | Examples |
|---|---|
| Ligand overlap | `overlap_volume` |
| Confidence-weighted overlap | `overlap_w_pos_volume` |
| Ensemble variance | `variance_avg` |
| PDE/PAE metrics | `complex_PDE_avg` |
| Openness metrics | `openess_avg` |
| Affinity metrics | `affinity_pred_value` |
| Ligand confidence | `ligand_pLDDT_avg` |

---

# Current Modeling Strategy

The current Structural Ensemble framework focuses on:

- Ensemble-based structural variability
- Confidence-aware structural descriptors
- Ligand consistency across models
- Multi-property prediction

Rather than relying on a single predicted structure, the workflow uses ensemble behavior across many Boltz-2 models as predictive signal.

---

# Future Directions

Planned additions include:

## Protein Language Models (PLM)

A parallel PLM-based prediction pipeline is currently under development.

Goals include:

- Combining SE and PLM predictions
- Building consensus/voting systems
- Prioritizing variants supported by multiple independent models

Future workflows may use:

- Structural Ensemble models
- Protein Language Models
- Consensus ranking systems

to prioritize variants for experimental testing.

---

# Notes

This repository is under active development and currently serves as both:

- A modeling framework
- An experimental biosensor engineering platform

Some modules and prediction systems are still evolving as additional experimental datasets become available.
