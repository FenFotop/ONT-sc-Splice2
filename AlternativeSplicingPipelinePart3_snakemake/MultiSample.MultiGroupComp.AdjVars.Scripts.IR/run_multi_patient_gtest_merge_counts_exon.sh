#!/bin/bash
#SBATCH --time=3-00:00:00
#SBATCH --mem=64g
#SBATCH --cpus-per-task=1
#SBATCH --output=logs/gtest_ir_%j.log

#########################################
### G-test version of the permutation script (intron retention / exon-based)
### Uses r2dtable null — single-threaded per SLURM job, B permutations only.

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

Rscript "$runfiles"/multi_patient_gtest_exon.R \
  $splitfiles $cell_metadata $comp_groups_column $cell_groups_column \
  $nperm $samples $outdir $outfile $groupa $groupb
