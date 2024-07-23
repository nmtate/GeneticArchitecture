#!/bin/bash -l
#SBATCH -t 24:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=6gb
#SBATCH --mail-type=ALL  
#SBATCH --mail-user=nmtate@umn.edu 
#SBATCH --job-name genearch
#SBATCH -o ./jobs/%j.genearch.out
#SBATCH -e ./jobs/%j.genearch.err
#SBATCH -p msismall,msilarge

# Check if the directory exists
if [ ! -d "/scratch.global/genearch_nmtate" ]; then
    # If it doesn't exist, create it
    mkdir -p /scratch.global/genearch_nmtate
    echo "Directory '/scratch.global/genearch_nmtate' created successfully."
else
    echo "Directory '/scratch.global/genearch_nmtate' already exists."
fi

conda activate heritability

snakemake -s genearch.smk \
    --profile slurm.genearch/ \
    --configfile config.yaml \
    --keep-going

