### merge outputs from parallelized permutation, combined patient
library(tidyverse)
library(tidyr)

args <- commandArgs(TRUE)
path.to.outputs <- args[1]
path.to.metadata <- args[2]
output <- args[3]


## test data
# path.to.outputs = "/gpfs/commons/groups/landau_lab/mariela/SRSF2_project/output/sicelore2/CF-1620.262447.364042.466631.merged.DTU/DTU.output.merged.sample.SRSF2.IMPvsEP/IMP_vs_EP/split_cluster_output"
# path.to.metadata = "/gpfs/commons/groups/landau_lab/mariela/SRSF2_project/output/sicelore2/CF-1620.262447.364042.466631.merged.DTU/DTU.output.merged.sample.SRSF2.IMPvsEP/combined_metadata/combined_metadata.csv"

path.to.three.outputs <- paste0(path.to.outputs, "/alt_three_prime/")
path.to.five.outputs <- paste0(path.to.outputs, "/alt_five_prime/")

## read in all outputs
setwd(path.to.three.outputs)
files <- list.files(path.to.three.outputs)
# files = list.files("/gpfs/commons/groups/landau_lab/mariela/SRSF2_project/output/sicelore2/CF-1620.262447.364042.466631.merged.DTU/DTU.output.merged.sample.SRSF2.IMPvsEP/IMP_vs_EP/split_cluster_output/alt_three_prime/")
three.output.list <- lapply(files, function(x) read.csv(file = x))
three.output <- do.call(rbind, three.output.list)

setwd(path.to.five.outputs)
files <- list.files(path.to.five.outputs)
five.output.list <- lapply(files, function(x) read.csv(file = x))
five.output <- do.call(rbind, five.output.list)
message("output files loaded")
#
# five.output = five.output %>% select(-five.group2.cluster.cov, -five.group1.cluster.cov)
# five.cluster.cov = five.output %>% group_by(three_prime_ID) %>% summarise(five.group1.cluster.cov = sum(obs.group1), five.group2.cluster.cov = sum(obs.group2))
# five.output = left_join(five.output, five.cluster.cov)

## load in metadata
metadata <- read.csv(path.to.metadata)
metadata$intron_junction <- paste(metadata$chr, metadata$start, metadata$end, sep = ":")

## remove the sample column and only focus on unique junctions
metadata <- metadata %>% select(-sample)
metadata <- distinct(metadata)
metadata$intron_junction <- paste(metadata$intron_junction, metadata$strand, sep = ":")
metadata <- metadata %>% select(-three_prime_ID, -five_prime_ID, -three_prime, -five_prime)

## merge output with metadata
three.merge.output <- inner_join(three.output, metadata)
five.merge.output <- inner_join(five.output, metadata)
message("Output merged with metadata")

## calculate PSI for group2 and group1

three.merge.output$group2.psi <- (three.merge.output$obs.group2 / three.merge.output$three.group2.cluster.cov) * 100
three.merge.output$group1.psi <- (three.merge.output$obs.group1 / three.merge.output$three.group1.cluster.cov) * 100
three.merge.output$dPSI <- three.merge.output$group1.psi - three.merge.output$group2.psi

five.merge.output$group2.psi <- (five.merge.output$obs.group2 / five.merge.output$five.group2.cluster.cov) * 100
five.merge.output$group1.psi <- (five.merge.output$obs.group1 / five.merge.output$five.group1.cluster.cov) * 100
five.merge.output$dPSI <- five.merge.output$group1.psi - five.merge.output$group2.psi

message("PSI calculated")

## write output
# setwd("/gpfs/commons/groups/landau_lab/SF3B1_splice_project/7.Permutation_junction_output/4.logOR.combined.analysis/1000_within_cell_type/MDS_P1_P2_P3")
setwd(output)
write.table(three.merge.output, "logOR_within_cell_type_ALT_3P_Junctions.txt")
write.table(five.merge.output, "logOR_within_cell_type_ALT_5P_Junctions.txt")
message("DONE!")
