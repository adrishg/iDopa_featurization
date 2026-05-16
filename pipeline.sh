#!/bin/bash
#SBATCH -p <partition>
#SBATCH --gres=gpu:1
#SBATCH -t 5-48:00:00
#SBATCH --job-name=boltz2-iDopa
#SBATCH --mem=125G
#SBATCH --mail-user=your-email@institution.edu
#SBATCH --mail-type=END

# ============================================================
# iDopa featurization pipeline
#
# Generates Boltz2 predictions for DOP and 5HT separately,
# then extracts structural features and merges into a final CSV.
#
# Fill in the CONFIGURE section before running.
# ============================================================

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
FEATURIZATION_DIR="$SCRIPT_DIR/featurization"

# ---- CONFIGURE ----
PROJECT_PATH=''    # root project directory
DATA_PATH=''       # input CSV (requires Tag and sequence columns)
NUM_MODELS=25
RECYCLES=0
# -------------------

# Environment — adapt to your cluster (conda env must have boltz installed)
eval "$(conda shell.bash hook)"
conda activate boltz2

# Derived paths
MODELS_DIR="$PROJECT_PATH/1-Models_Bz2"
ANALYZED_DIR="$PROJECT_PATH/3-AnalyzedData"

mkdir -p "$MODELS_DIR" "$ANALYZED_DIR"/{DOP,5HT}

# --- Step 1: Generate Boltz2 input YAMLs ---
python3 "$FEATURIZATION_DIR/csv2yamls_w_molecules.py" \
    "$DATA_PATH" \
    --output-dir "$PROJECT_PATH/iDopa_DOP" \
    --molecules '{"DOP": "C1=CC(=C(C=C1CCN)O)O"}'

python3 "$FEATURIZATION_DIR/csv2yamls_w_molecules.py" \
    "$DATA_PATH" \
    --output-dir "$PROJECT_PATH/iDopa_5HT" \
    --molecules '{"5HT": "C1=CC2=C(C=C1O)C(=CN2)CCN"}'

# --- Step 2: Run Boltz2 predictions ---
for LIGAND in DOP 5HT; do
    boltz predict "$PROJECT_PATH/iDopa_$LIGAND" \
        --output_format pdb \
        --use_msa_server \
        --out_dir "$MODELS_DIR/boltz_results_iDopa_$LIGAND" \
        --diffusion_samples "$NUM_MODELS" \
        --diffusion_samples_affinity "$NUM_MODELS" \
        --recycling_steps "$RECYCLES" \
        --use_potentials \
        --write_full_pae \
        --write_full_pde
done

# --- Step 3: Extract features ---
bash "$SCRIPT_DIR/analysis.sh" \
    "$MODELS_DIR/boltz_results_iDopa_DOP/predictions" \
    "$ANALYZED_DIR/DOP"

bash "$SCRIPT_DIR/analysis.sh" \
    "$MODELS_DIR/boltz_results_iDopa_5HT/predictions" \
    "$ANALYZED_DIR/5HT"
