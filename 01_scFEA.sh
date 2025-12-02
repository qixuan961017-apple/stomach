#!/bin/bash
#SBATCH --job-name=scFEA
#SBATCH -N 1
#SBATCH -n 8
#SBATCH --partition=Fnode2
#SBATCH --mem 800G
#SBATCH -o %A
#SBATCH -t 72:00:00

cd ~/miniconda3/bin/
source activate
conda activate scFEA
cd /public/home/b20213040320/scFEA/scFEA
python /public/home/b20213040320/scFEA/scFEA/src/scFEA.py \
    --data_dir /public/home/b20213040320/scFEA/scFEA/data \
    --input_dir /public/home/b20213040320/zhangyahui/01_scFEA \
    --moduleGene_file module_gene_m168.csv \
    --test_file count_matrix.csv \
    --cName_file cName_c70_m168.csv \
    --sc_imputation True \
    --stoichiometry_matrix cmMat_c70_m168.csv \
    --output_flux_file     /public/home/b20213040320/zhangyahui/01_scFEA/output/adj_flux.csv \
    --output_balance_file    /public/home/b20213040320/zhangyahui/01_scFEA/output/adj_balance.csv
