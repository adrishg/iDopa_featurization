#!/bin/bash
# ==========================================================
# Feature extraction pipeline for one set of Boltz2 predictions.
# Usage: bash analysis.sh <models_path> <output_path>
# ==========================================================

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
FEATURIZATION_DIR="$SCRIPT_DIR/featurization"

models_path="$1"
analyzed_path="$2"

if [ -z "$models_path" ] || [ -z "$analyzed_path" ]; then
    echo "Usage: bash analysis.sh <models_path> <output_path>"
    exit 1
fi

echo "=== Running analysis pipeline ==="
echo "Models:   $models_path"
echo "Output:   $analyzed_path"
echo

mkdir -p "$analyzed_path"

# === Step 1: Feature extraction ===
python3 "$FEATURIZATION_DIR/batch_LigOverlapVol.py" \
    --input_dir "$models_path" \
    --output_dir "$analyzed_path"

python3 "$FEATURIZATION_DIR/batch_distanceMaps_variance.py" \
    --input_dir "$models_path" \
    --output_dir "$analyzed_path"

python3 "$FEATURIZATION_DIR/getAffinities.py" \
    --input-dir "$models_path" \
    --output-csv "$analyzed_path/affinities.csv"

python3 "$FEATURIZATION_DIR/getOpenessDistances.py" \
    --parent-folder "$models_path" \
    --res1 40 --res2 389 --chain A \
    --output-csv "$analyzed_path/openess.csv"

# === Step 2: Merge feature CSVs ===
python3 "$FEATURIZATION_DIR/merge_csv_tags.py" \
    --primary_csv "$analyzed_path/overall_folder_summary.csv" \
    --secondary_csv "$analyzed_path/composite_variances.csv" \
    --output_csv "$analyzed_path/volumes_variances.csv" \
    --ref_column Tag \
    --columns_to_merge variance_avg variance_pLDDT_w complex_PDE_avg complex_PDE_var \
        complex_PDE_min complex_PDE_max PAE_avg PDE_avg

python3 "$FEATURIZATION_DIR/merge_csv_tags.py" \
    --primary_csv "$analyzed_path/volumes_variances.csv" \
    --secondary_csv "$analyzed_path/affinities.csv" \
    --output_csv "$analyzed_path/volumes_variances_affinities.csv" \
    --ref_column Tag \
    --columns_to_merge affinity_pred_value affinity_probability_binary

python3 "$FEATURIZATION_DIR/merge_csv_tags.py" \
    --primary_csv "$analyzed_path/volumes_variances_affinities.csv" \
    --secondary_csv "$analyzed_path/openess.csv" \
    --output_csv "$analyzed_path/volumes_variances_affinities_openess.csv" \
    --ref_column Tag \
    --columns_to_merge openess_avg openess_min openess_max openess_range

# === Step 3: Strip _DOP/_5HT suffix from Tag, add Ligand column ===
python3 - "$analyzed_path" <<'EOF'
import sys, re
import pandas as pd
from pathlib import Path

out_dir = Path(sys.argv[1])
csv_path = out_dir / "volumes_variances_affinities_openess.csv"

df = pd.read_csv(csv_path)
df["Ligand"] = df["Tag"].str.extract(r"_(DOP|5HT)$", expand=False)
df["Tag"] = df["Tag"].str.replace(r"_(DOP|5HT)$", "", regex=True)
out = csv_path.with_name(csv_path.stem + "_clean.csv")
df.to_csv(out, index=False)
print(f"Cleaned tag suffixes -> {out}")
EOF

echo "Results in: $analyzed_path"
