#!/bin/bash
#SBATCH -p gpu-vladimir        # Requesting to run on Anova-0
#SBATCH -t 5-48:00:00
#SBATCH --job-name=boltz2-iDopa_missPred_p5HT
#SBATCH --mem=125G
#SBATCH --mail-user=ahgonzalez@ucdavis.edu
#SBATCH --mail-type=END

source /share/yarovlab/ahgz/.bashrc

#The latest gcc version in barbera:
module load gcc/13.2.0

repo_path='/share/yarovlab/ahgz/Biosensor/iDopaSnFr/scripts/'
project_path='/share/yarovlab/ahgz/Biosensor/iDopaSnFr/'

#Example csv data file, only needs Tag and sequence
data_path='/share/yarovlab/ahgz/Biosensor/iDopaSnFr/0-Data/25-05-29_DB-Biosensor.csv' 


mkdir -p "$project_path"/{1-Models_Bz2,2-Data,3-AnalyzedData}
models_path='$project_path/1-Models_Bz2/boltz_results_iDopa/predictions/'
analyzed_data_path='/share/yarovlab/ahgz/Biosensor/iDopaSnFr/3-AnalyzedData/'

#Generates the yamls from the input data csv (by default generates one for 5HT and one for DOP)
# to change edit dictionary of ligands in generate_yamls_from_csv.py 
python3 csv2yamls_w_molecules.py \
    "$data_path/" \
    --output-dir "$project_path/iDopa" \
    --molecules '{"DOP": "C1=CC(=C(C=C1CCN)O)O", "5HT": "C1=CC2=C(C=C1O)C(=CN2)CCN"}' #Also can take json file (:

# Run Boltz2x 
sbatch --wait runBzprediction.sh \
    --input-dir \
    --output-dir \
    --num-models 25 \
    --recycles 0

#Analyzed models:
python3 /share/yarovlab/ahgz/Biosensor/iDopaSnFr/2-AnalysisNew/batch_LigOverlapVol.py \
   --input_dir $models_path \
   --output_dir $analyzed_data_path

python3 /share/yarovlab/ahgz/Biosensor/iDopaSnFr/2-AnalysisNew/batch_distanceMaps_variance.py \
   --input_dir $models_path \
   --output_dir $analyzed_data_path

python3 /share/yarovlab/ahgz/Biosensor/iDopaSnFr/2-AnalysisNew/getAffinities.py \
   --input-dir $models_path \
   --output-csv "$analyzed_data_path/affinities.csv"

python3 /share/yarovlab/ahgz/Biosensor/iDopaSnFr/2-AnalysisNew/getOpenessDistances.py \
    --parent-folder $models_path \
    --res1 40 --res2 389 \
    --chain A \
    --output-csv "$analyzed_data_path/openess.csv"

#Merge all data in oe csv
python3 /share/yarovlab/ahgz/Biosensor/iDopaSnFr/2-AnalysisNew/merge_csv_tags.py \
    --primary_csv "$analyzed_data_path/overall_folder_summary.csv" \
    --secondary_csv "$analyzed_data_path/composite_variances.csv" \
    --output_csv "$analyzed_data_path/volumes_variances.csv" \
    --ref_column 'Tag' \
    --columns_to_merge 'variance_avg' 'variance_pLDDT_w' 'complex_PDE_avg' 'complex_PDE_var' 'complex_PDE_min' 'complex_PDE_max' 'PAE_avg' 'PDE_avg'

python3 /share/yarovlab/ahgz/Biosensor/iDopaSnFr/2-AnalysisNew/merge_csv_tags.py \
    --primary_csv "$analyzed_data_path/volumes_variances.csv" \
    --secondary_csv "$analyzed_data_path/affinities.csv" \
    --output_csv "$analyzed_data_path/volumes_variances_affinities.csv" \
    --ref_column 'Tag' \
    --columns_to_merge 'affinity_pred_value' 'affinity_probability_binary'

python3 /share/yarovlab/ahgz/Biosensor/iDopaSnFr/2-AnalysisNew/merge_csv_tags.py \
    --primary_csv "$analyzed_data_path/volumes_variances_affinities.csv" \
    --secondary_csv "$analyzed_data_path/openess.csv" \
    --output_csv "$analyzed_data_path/volumes_variances_affinities_openess.csv" \
    --ref_column 'Tag' \
    --columns_to_merge 'openess_avg' 'openess_min' 'openess_max' 'openess_range'
# Until this point Tag contains _5HT and _DOP need an step to get rif of it and considers that comes from _DOP and _5HT models


