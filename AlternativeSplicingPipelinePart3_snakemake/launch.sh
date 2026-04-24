#!/bin/bash
#SBATCH --mail-type=ALL
#SBATCH --job-name=splicing_comparisons
#SBATCH --partition=cpu
#SBATCH --output=%x.log
#SBATCH --error=%x.err
#SBATCH --time=5-00:00:00
#SBATCH --cpus-per-task=5


conda activate ont-sc-splice
module load python 

snakemake --profile new_slurm_profile 
          --configfile config.yaml 
