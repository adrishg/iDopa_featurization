#!/bin/bash
#SBATCH -p gpu-vladimir        # Requesting to run on Anova-0
#SBATCH --gres=gpu:1           # Request 1 GPUs
#SBATCH -t 3-48:00:00
#SBATCH --job-name=boltz2-iDopa
#SBATCH --mem=125G
#SBATCH --mail-user=ahgonzalez@ucdavis.edu
#SBATCH --mail-type=END

# ============================================================
# Example usage:
# sbatch --wait runBzprediction.sh \
#     --input-dir /path/to/input \
#     --output-dir /path/to/output \
#     --num-models 25 \
#     --recycles 0
# ============================================================

INPUT_DIR=""
OUT_DIR=""
#By default 25 models 0 recycles
NUM_MODELS=25
RECYCLES=0

usage() {
  echo "Usage: sbatch [sbatch options] $0 --input-dir <dir> --output-dir <dir> [--num-models 25] [--recycles 0]"
  exit 1
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --input-dir)
      INPUT_DIR="${2:-}"; shift 2 ;;
    --output-dir)
      OUT_DIR="${2:-}"; shift 2 ;;
    --num-models)
      NUM_MODELS="${2:-}"; shift 2 ;;
    --recycles)
      RECYCLES="${2:-}"; shift 2 ;;
    -*|--*)
      echo "Unknown option: $1"; usage ;;
    *)
      echo "Unexpected positional argument: $1"; usage ;;
  esac
done

[[ -z "$INPUT_DIR" || -z "$OUT_DIR" ]] && usage

# -------------------------
# Environment setup
# -------------------------
source /share/yarovlab/ahgz/.bashrc
eval "$(conda shell.bash hook)"   # activates conda shell support
conda activate boltz2

# Just printing things in the slurm
echo "Using GPU (SLURM_JOB_GPUS): ${SLURM_JOB_GPUS:-<unset>}"
echo "CUDA_VISIBLE_DEVICES as seen by job: ${CUDA_VISIBLE_DEVICES:-<unset>}"
nvidia-smi 

mkdir -p "$OUT_DIR"

# -------------------------
# Run Boltz prediction
# -------------------------
echo "Running boltz2 prediction:"
echo "  Input dir   : $INPUT_DIR"
echo "  Output dir  : $OUT_DIR"
echo "  Num models  : $NUM_MODELS"
echo "  Recycles    : $RECYCLES"

CUDA_VISIBLE_DEVICES="${CUDA_VISIBLE_DEVICES:-0}" boltz predict "$INPUT_DIR" \
  --output_format pdb \
  --use_msa_server \
  --out_dir "$OUT_DIR" \
  --diffusion_samples "$NUM_MODELS" \
  --diffusion_samples_affinity "$NUM_MODELS" \
  --recycling_steps "$RECYCLES" \
  --use_potentials \
  --write_full_pae \
  --write_full_pde

echo "Boltz2 prediction complete."

