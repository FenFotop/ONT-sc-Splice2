#!/bin/bash

#SBATCH --mem=100g
#SBATCH --partition=bigmem,pe2
#SBATCH --cpus-per-task=10

module load R/3.6.0


### arguments input for this file

workdir=$1
run_files=$2
samples=$3
nperm=$4
job_name=$5
comp_groups_column=$6
cell_groups_column=$7
groupList=$8 #e.g. multiple cell types. Format: ""HSPC" "MEP""
outputdir=$9 #DTU.output.merged.sample.2UMI
junc_cov_thresh=${10:-5} ##If nothing is set here, set the default to 5 for the read coverage threshold 

#run_files="$run_files"/Splicing.Pipeline.All.scripts.AdjVars/MultiSample.MultiGroupComp.AdjVars.Scripts

counts="$workdir"/JuncCounts_matrix
cell_metadata="$workdir"/Cell_metadata
junc_metadata="$workdir"/Junction_metadata

## input folder should have three folders
# 1.JuncCounts_matrix - all counts matrix for all patients being combined
# 2.Cell_metadata - all cell metadata tables for all patients
# 3.Junction_metadata - strand adjusted junction metadata for all patients 

###############################################################
######## Step 1: Combine Metadata to make Common Cluster ID's ############
### Input: metadata output from leafcutter annotation for all patients
### Input: path to output files/ working directory 
### Output: Combined output

## format: Rscript strand_adjustment.R <metadata> <output directory> 

###############################################################

#echo $workdir 
cd $workdir 
mkdir $outputdir
cd ./$outputdir
mkdir combined_metadata

Rscript "$run_files"/combine_metadata_exon.R $junc_metadata "$workdir"/"$outputdir"/combined_metadata $samples


##############################################################
######### Step 2: Split clusters for differential transcript usage
### Input: full counts matrix from leafcutter junction calling
### Input: Genotype matrix
### Input: Strand adjusted metadata
### Input: Pattern (i.e. barcode pattern, "_1", "_2", or "_3")
### Input: Output directory 

## format: Rscript split_clusters_v2.R <counts> <genotype> <metadata> <pattern> <output>

###############################################################


for type in ${groupList[*]};
do
	mkdir "$type"
	cd ./"$type"

	mkdir split_cluster_files
	cd ./split_cluster_files

	for i in {1..100}
	do 
	
	mkdir split_"$i"
	mkdir split_"$i"/counts_files
	mkdir split_"$i"/data_tables
	
	done

	Rscript "$run_files"/split_clusters_exon.R $counts $cell_metadata $comp_groups_column $cell_groups_column $type $junc_metadata "$workdir"/"$outputdir"/"$type"/split_cluster_files $samples $junc_cov_thresh



###########################################################
######## Step 3: Batch submit each split cluster for differential analysis
### Input: path to split files 
### Input: path to genotype matrix 
### Input: Number of permutations
### Input: output directory 
### Input: output file name 
### output: differential transcript table for 3p and 5p in two separate folders

## format: sbatch run_split_perm_within_celltype_5p_3p.sh <path to split files> <genotype> <nperm> <output.dir> <output.file> 

###########################################################

	cd ..
	mkdir split_cluster_output
	mkdir split_cluster_output/alt_three_prime
	mkdir logs

	permute_jobids=()
	for i in {1..100}; do
	permute_jobids+=($(sbatch --job-name="$job_name" "$run_files"/run_multi_patient_permute_merge_counts_exon.sh "$workdir"/"$outputdir"/"$type"/split_cluster_files/split_"$i" $cell_metadata $comp_groups_column $cell_groups_column $nperm $samples "$workdir"/"$outputdir"/"$type"/split_cluster_output output_"$i" "$run_files"))
	done 


###########################################################
####### Step 4: Merge final output into one file and merge with all annotation information 
### Input: run_files path 
### Input: outputs directory where all files are stored
### Input: strand adjusted metadata 
### Input: Final outfile 

	mkdir merge_final_output

	merge=($(sbatch --dependency=singleton --job-name="$job_name" "$run_files"/run_merge_combine_output_exon.sh "$run_files" "$workdir"/"$outputdir"/"$type"/split_cluster_output "$workdir"/"$outputdir"/combined_metadata/combined_metadata.csv "$workdir"/"$outputdir"/"$type"/merge_final_output))

        #Don't have another layer of annotation here (already did IR annotation 
	#annotate=($(sbatch --dependency=singleton --job-name="$job_name" "$run_files"/run_DTU_annotation.new.sh "$workdir"/"$outputdir"/"$type"/merge_final_output))

cd ..
done

echo "Done!" 

