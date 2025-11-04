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
data_path='/share/yarovlab/ahgz/Biosensor/iDopaSnFr/0-Data/25-11-04_Top-92-sequences.csv' 


mkdir -p "$project_path"/{1-Models_Bz2,2-Data,3-AnalyzedData}

analyzed_data_path='/share/yarovlab/ahgz/Biosensor/iDopaSnFr/3-AnalyzedData/'

#Generates the yamls from the input data csv (by default generates one for 5HT and one for DOP)
# to change edit dictionary of ligands in generate_yamls_from_csv.py 
python3 csv2yamls_w_molecules.py \
    "$data_path/" \
    --output-dir "$project_path/iDopa_DOP" \
    --molecules '{"DOP": "C1=CC(=C(C=C1CCN)O)O"}' #Also can take json file (:

python3 csv2yamls_w_molecules.py \
    "$data_path/" \
    --output-dir "$project_path/iDopa_5HT" \
    --molecules '{"5HT": "C1=CC2=C(C=C1O)C(=CN2)CCN"}' #Also can take json file (:


# Run Boltz2x 
sbatch --wait runBzprediction.sh \
    --input-dir "$project_path/iDopa_DOP" \
    --output-dir \
    --num-models 25 \
    --recycles 0

sbatch --wait runBzprediction.sh \
    --input-dir "$project_path/iDopa_5HT" \
    --output-dir \
    --num-models 25 \
    --recycles 0

models_path_DOP='$project_path/1-Models_Bz2/boltz_results_iDopa_DOP/predictions/'
models_path_5HT='$project_path/1-Models_Bz2/boltz_results_iDopa_5HT/predictions/'

analyzed_data_path_DOP='/share/yarovlab/ahgz/Biosensor/iDopaSnFr/3-AnalyzedData/DOP/'
analyzed_data_path_5HT='/share/yarovlab/ahgz/Biosensor/iDopaSnFr/3-AnalyzedData/5HT'

mkdir "$analyzed_data_path_DOP"
mkdir "$analyzed_data_path_5HT"

sbatch --wait "$repo_path/analysis.sh" $models_path_DOP $analyzed_data_path_DOP
sbatch --wait "$repo_path/analysis.sh" $models_path_5HT $analyzed_data_path_5HT
