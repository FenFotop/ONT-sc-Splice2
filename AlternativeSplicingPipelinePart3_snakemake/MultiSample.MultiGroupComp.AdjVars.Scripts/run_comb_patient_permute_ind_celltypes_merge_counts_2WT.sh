#!/bin/bash
#SBATCH --time=3-00:00:00
#SBATCH --mem=64g
#SBATCH --cpus-per-task=10
#SBATCH --output=logs/combined_permutations_%j.log 

#########################################
### This script submits the script that permutes the combined patient data by merging all counts before calculating log odds ratio 
### permutations are done within each patient for each cell type matrix that is submitted to give one output 

splitfiles=$1
cell_metadata=$2
comp_groups_column=$3 #e.g. genotype labels
cell_groups_column=$4 #e.g. column with cell type labels 
nperm=$5
samples=$6
outdir=$7
outfile=$8
runfiles=$9
groupa=${10}
groupb=${11}

echo $splitfiles
echo $cell_metadata


Rscript "$runfiles"/combined_patient_permute_merge_reads_ind_celltypes_2WT.R $splitfiles $cell_metadata $comp_groups_column $cell_groups_column $nperm $samples $outdir $outfile $groupa $groupb
