# iDopa_featurization
This repository provides a collection of scripts for the automated generation, modeling, and analysis of iDopaSnFR biosensor variants using Boltz2, an all-atom diffusion model for joint protein–ligand structure prediction.

The pipeline starts from a CSV file containing variant information and proceeds through YAML generation, Boltz2 modeling, and structural feature extraction to produce a unified dataset describing each biosensor variant and its interactions with small molecules.

Overview

The workflow is designed to:

1. Generate Boltz2 input YAMLs
From a CSV file listing the biosensor variants (columns Tag and sequence), the generate_yamls_from_csv.py script creates YAML configuration files pairing each protein variant with one or more ligands.

Ligands are defined via a molecules dictionary (SMILES strings).

Current ligands include dopamine (DOP) and serotonin (5HT).

Each YAML includes a properties block specifying the binder chain.

2. Run Boltz2 predictions
The SLURM script (runBzprediction.sh) submits Boltz2 jobs on the cluster.
Each job:

Runs 25 diffusion samples using 0 recycling steps by default (can be parameterized).

Produces structural models for both dopamine and serotonin bound states.

3. Extract structural and model-based features
The analysis phase gathers predicted and computed descriptors, including:

Model metrics: predicted ligand pLDDT, affinity, and confidence scores.

4. Structural features: openness of the PBP domain, overlap volume between predicted ligand positions, and per-variant variability across models and variance across the ensemble of 25 predictions.

5. Compile final tables
The output is a single CSV summarizing model features for all variants and ligands.
These tables serve as inputs for downstream analysis, ranking, or visualization of biosensor performance and conformational diversity.