### create combined metadata across all three samples 
library(tidyverse)

args = commandArgs(TRUE)
path.to.metadata = args[1]
output.dir = args[2]
sample.names = args[3]

sample.names = unlist(strsplit(sample.names, split = ","))
print(sample.names)

## test data
#path.to.metadata = "/gpfs/commons/groups/landau_lab/mariela/SRSF2_project/output/sicelore2/CF-1620.262447.364042.466631.merged.DTU/3.Junction_metadata"
#sample.names = c("CF-1620.262447","CF-1620.364042","CF-1620.466631")


## load in metadata
setwd(path.to.metadata)
files = grep(paste(sample.names, collapse = "|"), list.files(path.to.metadata), value = T)
metadata.list = lapply(files, function(x) read.csv(file = x, stringsAsFactors = F, sep=","))
positions <- sapply(sample.names, function(sample) grep(sample, files))
names(metadata.list)[positions] <- sample.names

## identify strand adjusted start and end site
for (sample in sample.names) {
  metadata.list[[sample]]$five_prime = "NA"
  metadata.list[[sample]]$three_prime = "NA"
  metadata.list[[sample]] = metadata.list[[sample]] %>% mutate(five_prime = ifelse(strand == "+", start, end),
                                      three_prime = ifelse(strand =="+", end,start))
}

## select only parts of metadata that we care about
for (sample in sample.names) {
  metadata.list[[sample]] = metadata.list[[sample]] %>% select(chr, start, end, strand, gene, five_prime, three_prime, verdict, startClass, endClass, fivep_class, threep_class, fivep_distance, threep_distance, relStartExon_skipping, relEndExon_skipping)
  metadata.list[[sample]]$sample <- sample
}

#combine into one list of all intron junctions and create new three prime and five prime cluster IDs 
metadata = do.call(rbind, metadata.list)
metadata$five_prime_ID = metadata %>% group_by(chr, five_prime, strand) %>% group_indices()
metadata$three_prime_ID = metadata %>% group_by(chr, three_prime, strand) %>% group_indices()

metadata$five_prime_ID = gsub("^", "clu_", metadata$five_prime_ID)
metadata$three_prime_ID = gsub("^", "clu_", metadata$three_prime_ID)
#metadata$sample = rownames(metadata)
#metadata$sample = gsub("\\..*", "", metadata$sample)

## write output
setwd(output.dir)
write.csv(metadata, "combined_metadata.csv", quote = FALSE, row.names = FALSE)

#write.csv(metadata, "/gpfs/commons/home/ahawkins/MDS_ONT_splice/Junction_matrix_processing/MDS_combined/6.combined_metadata/MDS_combined_metadata.txt", quote = FALSE)
#metadata = read.csv("/gpfs/commons/home/ahawkins/MDS_ONT_splice/Junction_matrix_processing/MDS_combined/6.combined_metadata/MDS_combined_metadata.txt", stringsAsFactors = FALSE, row.names = 1)

