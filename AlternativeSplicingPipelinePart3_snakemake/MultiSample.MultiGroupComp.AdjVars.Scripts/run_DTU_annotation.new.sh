#!/bin/bash
#SBATCH --time=3-00:00:00
#SBATCH --mem=32g
#SBATCH --job-name=submit_annotation


workingdir=$1
run_files=$2

Rscript "$run_files"/DTU_junction_annotation.new.R $workingdir

