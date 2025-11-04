#!/bin/bash
# ==========================================================
# Simple pipeline orchestrator for biosensor model analysis
# Usage: bash analysis.sh <models_path> <output_path>
# ==========================================================

set -e  # stop on first error

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

# === Step 1: Run analysis scripts ===
python3 /share/yarovlab/ahgz/Biosensor/iDopaSnFr/2-AnalysisNew/batch_LigOverlapVol.py \
    --input_dir "$models_path" \
    --output_dir "$analyzed_path"

python3 /share/yarovlab/ahgz/Biosensor/iDopaSnFr/2-AnalysisNew/batch_distanceMaps_variance.py \
    --input_dir "$models_path" \
    --output_dir "$analyzed_path"

python3 /share/yarovlab/ahgz/Biosensor/iDopaSnFr/2-AnalysisNew/getAffinities.py \
    --input-dir "$models_path" \
    --output-csv "$analyzed_path/affinities.csv"

python3 /share/yarovlab/ahgz/Biosensor/iDopaSnFr/2-AnalysisNew/getOpenessDistances.py \
    --parent-folder "$models_path" \
    --res1 40 --res2 389 --chain A \
    --output-csv "$analyzed_path/openess.csv"

# === Step 2: Merge CSVs ===
python3 /share/yarovlab/ahgz/Biosensor/iDopaSnFr/2-AnalysisNew/merge_csv_tags.py \
    --primary_csv "$analyzed_path/overall_folder_summary.csv" \
    --secondary_csv "$analyzed_path/composite_variances.csv" \
    --output_csv "$analyzed_path/volumes_variances.csv" \
    --ref_column 'Tag' \
    --columns_to_merge variance_avg variance_pLDDT_w complex_PDE_avg complex_PDE_var complex_PDE_min complex_PDE_max PAE_avg PDE_avg

python3 /share/yarovlab/ahgz/Biosensor/iDopaSnFr/2-AnalysisNew/merge_csv_tags.py \
    --primary_csv "$analyzed_path/volumes_variances.csv" \
    --secondary_csv "$analyzed_path/affinities.csv" \
    --output_csv "$analyzed_path/volumes_variances_affinities.csv" \
    --ref_column 'Tag' \
    --columns_to_merge affinity_pred_value affinity_probability_binary

python3 /share/yarovlab/ahgz/Biosensor/iDopaSnFr/2-AnalysisNew/merge_csv_tags.py \
    --primary_csv "$analyzed_path/volumes_variances_affinities.csv" \
    --secondary_csv "$analyzed_path/openess.csv" \
    --output_csv "$analyzed_path/volumes_variances_affinities_openess.csv" \
    --ref_column 'Tag' \
    --columns_to_merge openess_avg openess_min openess_max openess_range

# === Step 3: Clean _DOP/_5HT suffix from Tag ===
python3 - <<'EOF'
import pandas as pd
from pathlib import Path

analyzed = Path("$analyzed_path")
csv_path = analyzed / "volumes_variances_affinities_openess.csv"

df = pd.read_csv(csv_path)
df["Ligand"] = df["Tag"].str.extract(r"_(DOP|5HT)$", expand=False)
df["Tag"] = df["Tag"].str.replace(r"_(DOP|5HT)$", "", regex=True)
out = csv_path.with_name(csv_path.stem + "_clean.csv")
df.to_csv(out, index=False)
print(f"leaned tag suffixes → {out}")

echo "Results in: $analyzed_path"
