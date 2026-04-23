#!/bin/bash
#SBATCH --time=3-00:00:00
#SBATCH --mem=64g
#SBATCH --cpus-per-task=1
#SBATCH --output=logs/gtest_combined_%j.log

#########################################
### G-test version of the permutation script (splice junctions)
### Uses r2dtable null — single-threaded per SLURM job, B permutations only.
### Much faster than the logOR label-shuffle version: no per-permutation
### matrix resampling, just r2dtable calls on pre-aggregated marginals.

splitfiles=$1
cell_metadata=$2
comp_groups_column=$3
cell_groups_column=$4
nperm=$5           # interpreted as B (number of G-test permutations)
samples=$6
outdir=$7
outfile=$8
runfiles=$9
groupa=${10}
groupb=${11}

echo $splitfiles
echo $cell_metadata

Rscript "$runfiles"/combined_patient_gtest_ind_celltypes_2WT.R \
  $splitfiles $cell_metadata $comp_groups_column $cell_groups_column \
  $nperm $samples $outdir $outfile $groupa $groupb
