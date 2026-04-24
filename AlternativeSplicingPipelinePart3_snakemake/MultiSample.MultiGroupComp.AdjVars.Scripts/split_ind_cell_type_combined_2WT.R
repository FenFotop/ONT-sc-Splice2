library(Matrix)
library(parallel)
library(tidyverse)
library(matrixStats)
library(data.table)

args <- commandArgs(TRUE)
path.to.count.matrix <- args[1]
path.to.cell.meta <- args[2]
comp.groups.column <- args[3] # e.g. genotype labels
cell.groups.column <- args[4] # e.g. column with cell type labels
cell.group <- args[5] ## specific cell type group e.g. one cell type
path.to.junc.meta <- args[6]
output.dir <- args[7]
sample.names <- args[8]
cov.thresh <- as.integer(args[9]) ## set to 5 as default
group1_name <- args[10]
group2_name <- args[11]
nchunk <- as.integer(args[12])

print(args)
# if(is.na(cov.thresh)){
#  cov.thresh = 5
# }

print("junc coverage threshold")
print(cov.thresh)

## ---- Test data ---- ##
# path.to.count.matrix = "/gpfs/commons/groups/landau_lab/mariela/SRSF2_project/output/sicelore2/CF-1620.262447.364042.466631.merged.DTU/1.JuncCounts_matrix/"
# path.to.cell.meta = "/gpfs/commons/groups/landau_lab/mariela/SRSF2_project/output/sicelore2/CF-1620.262447.364042.466631.merged.DTU/2.Cell_metadata"
# comp.groups.column = "HSCSvsEP"
# cell.groups.column = "major.cell.groups"
# cell.group = "HSCS_vs_EP"
# sample.names = c("CF-1620.262447","CF-1620.364042","CF-1620.466631")
# path.to.junc.meta = "/gpfs/commons/groups/landau_lab/mariela/SRSF2_project/output/sicelore2/CF-1620.262447.364042.466631.merged.DTU/3.Junction_metadata/"
# sample.names = c("CF-1620.262447")
# cov.thresh = 5

## -------- ##

chunk <- function(x, n) split(x, factor(sort(rank(x) %% n)))

sample.names <- unlist(strsplit(sample.names, split = ","))
print(sample.names)

## load in matrix
print(path.to.count.matrix)
setwd(path.to.count.matrix)
# files = list.files(path.to.count.matrix)
files <- grep(paste(sample.names, collapse = "|"), list.files(path.to.count.matrix), value = T)

mtx.list <- lapply(files, function(x) fread(file = x))
positions <- sapply(sample.names, function(sample) grep(sample, files))
names(mtx.list)[positions] <- sample.names
message("Matrix loaded")


print(path.to.junc.meta)
setwd(path.to.junc.meta)
# files = list.files(path.to.junc.meta)
files <- grep(paste(sample.names, collapse = "|"), list.files(path.to.junc.meta), value = T)

metadata.list <- lapply(files, function(x) read.csv(file = x, stringsAsFactors = F))
positions <- sapply(sample.names, function(sample) grep(sample, files))
names(metadata.list)[positions] <- sample.names

## Should we remove the junctions where the strand information doesn't match here???

## Strand - from the junction annotation process - i.e. the strand info for the gene based on what's in the reference
## Strand_1 - from the junction x cell count matrix rownames i.e. straight from the read information
# By using the opposite ones, it removes any junctions that don't match for their strand assingment!

for (sample in sample.names) {
  mtx.list[[sample]] <- setDF(mtx.list[[sample]])
  rownames(mtx.list[[sample]]) <- mtx.list[[sample]]$intron_junction
  # Remove extra junctions where the strand information doens't match up , are we looking counts from that??
  metadata.list[[sample]] <- metadata.list[[sample]][which(metadata.list[[sample]]$strand == metadata.list[[sample]]$strand_1), ]
  metadata.list[[sample]] <- metadata.list[[sample]] %>% mutate(intron_junction = paste(chr, start, end, strand, sep = ":"))
  shared <- intersect(rownames(mtx.list[[sample]]), metadata.list[[sample]]$intron_junction)
  mtx.list[[sample]] <- mtx.list[[sample]][shared, ]
  metadata.list[[sample]] <- metadata.list[[sample]][which(metadata.list[[sample]]$intron_junction %in% shared), ]
  print(dim(mtx.list[[sample]]))
  print(dim(metadata.list[[sample]]))
}


# load genotype information
# genotype.list = list()
cell.meta.list <- list()
setwd(path.to.cell.meta)

cell.meta.files <- grep(paste(sample.names, collapse = "|"), list.files(path.to.cell.meta), value = T)

print(cell.meta.files)

## Removing cell type and cell type pattern requirement - should I add in the permutation group requirement and only include cell


for (file in cell.meta.files) {
  cell.meta <- as.data.frame(read.table(file, stringsAsFactors = F, header = TRUE, fill = TRUE, sep = "\t"))
  id <- gsub("_metadata.txt", "", file)
  cell.meta <- cell.meta[, c(cell.groups.column, comp.groups.column)]
  print("Removing Pattern on cell barcode")
  rownames(cell.meta) <- unlist(lapply(strsplit(rownames(cell.meta), split = "_"), "[", 1))
  sum(rownames(cell.meta) %in% colnames(mtx.list[[id]]))
  cell.meta <- cell.meta[rownames(cell.meta) %in% colnames(mtx.list[[id]]), ]
  cell.meta.sub <- cell.meta[!is.na(cell.meta[, comp.groups.column]), ]

  cell.meta.list[[id]] <- cell.meta.sub
  print(cell.meta.list[[id]][1:5, ])
  print(dim(cell.meta.list[[id]]))
}

message("Genotype loaded")


## Create combined metadata with combined three prime ID and five prime ID
combined_data <- list()
for (sample in sample.names) {
  combined_data[[sample]] <- data.frame(intron_junction = rownames(mtx.list[[sample]]))
  combined_data[[sample]] <- combined_data[[sample]] %>%
    separate(intron_junction, into = c("chr", "start", "end", "strand_1"), sep = ":") %>%
    mutate(intron_junction = paste(chr, start, end, strand_1, sep = ":"))
  combined_data[[sample]]$sample <- sample
}


combined <- do.call(rbind, combined_data)

# Add sample same before hand
# combined$sample = rownames(combined)
# combined$sample = gsub("\\..*", "", combined$sample)

combined$five_prime <- "NA"
combined$three_prime <- "NA"
combined <- combined %>% mutate(
  five_prime = ifelse(strand_1 == "+", start, end),
  three_prime = ifelse(strand_1 == "+", end, start)
)
combined$five_prime_ID <- combined %>%
  group_by(chr, five_prime, strand_1) %>%
  group_indices()
combined$three_prime_ID <- combined %>%
  group_by(chr, three_prime, strand_1) %>%
  group_indices()
combined$five_prime_ID <- gsub("^", "clu_", combined$five_prime_ID)
combined$three_prime_ID <- gsub("^", "clu_", combined$three_prime_ID)

combined <- combined %>% mutate(unique_ID = paste(chr, start, end, strand_1, three_prime_ID, five_prime_ID, sep = ":"))

# filename = paste(output.dir, "/combined.data.txt", sep = "")
# write.table(combined, filename)
# dim(combined)

message("sample cluster ID's combined")


data.comb.list <- list()


for (sample in sample.names) {
  group1 <- rownames(cell.meta.list[[sample]])[cell.meta.list[[sample]][, comp.groups.column] == group1_name & cell.meta.list[[sample]][, cell.groups.column] == cell.group]
  group2 <- rownames(cell.meta.list[[sample]])[cell.meta.list[[sample]][, comp.groups.column] == group2_name & cell.meta.list[[sample]][, cell.groups.column] == cell.group]

  # mut = rownames(genotype.list[[sample]])[genotype.list[[sample]]$Final_Genotype == "MUT" & genotype.list[[sample]]$Cell.Assignment_2 == celltype]

  data <- combined[which(combined$sample == sample), ]
  data$obs.group1 <- rowSums(mtx.list[[sample]][, colnames(mtx.list[[sample]]) %in% group1, drop = FALSE])
  data$obs.group2 <- rowSums(mtx.list[[sample]][, colnames(mtx.list[[sample]]) %in% group2, drop = FALSE])
  data$total.reads.per.junction <- data$obs.group1 + data$obs.group2
  data.comb.list[[sample]] <- data
  data.comb.list[[sample]]$sample <- sample
}


data.comb <- bind_rows(data.comb.list)
## After joining the 2 count dataframes, we are now merging the counts into 1
data.comb.counts <- data.comb %>%
  group_by(intron_junction) %>%
  summarise(obs.group1 = sum(obs.group1), obs.group2 = sum(obs.group2), total.reads.per.junction = sum(total.reads.per.junction))
data.comb <- data.comb %>% select(-obs.group1, -obs.group2, -total.reads.per.junction, -sample)
## After merging the counts, we can now remove redundant rows (as junctions that appeared in multiple samples will be recorded miltiple times)
data.comb <- distinct(data.comb)
data.comb <- left_join(data.comb, data.comb.counts)

## total read per junction threhsold variable here
data.comb <- data.comb[which(data.comb$total.reads.per.junction >= cov.thresh), ]
dim(data.comb)

# Correct constitutive score:
metadata.list <- map(metadata.list, function(x) {
  if (class(x$constitutive.score) != "character") {
    x$constitutive.score <- ifelse(is.na(x$constitutive.score), ".", x$constitutive.score)
  }
  return(x)
})
comb.metadata <- distinct(Reduce(full_join, metadata.list) %>% select(intron_junction, startClass, endClass))
dim(comb.metadata)
# comb.metadata = distinct(Reduce(function(x, y) full_join(x, y, by = c("intron_junction", "startClass", "endClass")), metadata.list) %>%
#  select(intron_junction, startClass, endClass))

data.comb <- inner_join(data.comb, comb.metadata)
# idx = which(data.comb$startClass == "" & data.comb$endClass == "")
# data.comb[-idx,]
# dim(data.comb)

clusters <- data.comb %>%
  group_by(five_prime_ID) %>%
  tally()
filt.clusters <- clusters %>% filter(n != 1)
clust.three <- filt.clusters$five_prime_ID
# alt.clust.three = data.comb[which(data.comb$endClass == "not_main_3_prime" | data.comb$startClass == "not_main_3_prime"), "five_prime_ID"]
# alt.clust.three = intersect(alt.clust.three, clust.three)
alt.clust.three <- clust.three
print("Total 3' clusters after filter:")
length(alt.clust.three)

clusters <- data.comb %>%
  group_by(three_prime_ID) %>%
  tally()
filt.clusters <- clusters %>% filter(n != 1)
clust.five <- filt.clusters$three_prime_ID
# alt.clust.five = data.comb[which(data.comb$endClass == "not_main_5_prime" | data.comb$startClass == "not_main_5_prime"), "three_prime_ID"]
# alt.clust.five = intersect(alt.clust.five, clust.five)
alt.clust.five <- clust.five

print("Total 5' clusters after filter:")
length(alt.clust.five)


# Check the junction exist in all cells 
#shared_intron_junctions <- Reduce(intersect, sapply(mtx.list, rownames))
shared_intron_junctions <- if(length(mtx.list) == 1){
  rownames(mtx.list[[1]])
  }else{
    Reduce(intersect, sapply(mtx.list, rownames))
}
## filter for three prime junctions
data.filt.three <- data.comb[which(data.comb$five_prime_ID %in% alt.clust.three), ]

## filter for five prime junctions
data.filt.five <- data.comb[which(data.comb$three_prime_ID %in% alt.clust.five), ]



message("data filtered")

setwd(output.dir)

if (min(length(alt.clust.three),length(alt.clust.five)) == 0) {
  message("Not enough IR events to split.")
} else {
  # Determine how to split the clusters
  if (min(length(alt.clust.three) , length(alt.clust.five)) < nchunk) {
     split_clusters_three <- alt.clust.three
     split_clusters_five <- alt.clust.five

  } else {
    split_clusters_three <- chunk(alt.clust.three, nchunk)
    split_clusters_five <- chunk(alt.clust.five, nchunk)

  }
  
  for (i in 1:length(split_clusters_three)) {
    for (sample in sample.names) {
      data.split <- data.filt.three[which(data.filt.three$five_prime_ID %in% split_clusters_three[[i]]), ]
      split_junc <- unique(as.character(data.split$intron_junction))
      # prevent NAs 
      split_junc <- split_junc[split_junc %in% shared_intron_junctions] 
      data.split <- data.split %>% 
        filter(intron_junction %in% split_junc) %>% 
        distinct(intron_junction, .keep_all = TRUE)
      
      # Check if data.split is empty
      if (nrow(data.split) == 0) {
        next  # Skip to the next iteration of the outer loop
      }
      
      mtx.split <- mtx.list[[sample]][split_junc, ]
      workdir <- paste("./split_", i, "/three_prime/data_tables/", sep = "")
      if (!dir.exists(workdir)) dir.create(workdir, recursive = TRUE)
      
      filename <- paste(sample, "data.filt", i, "csv", sep = ".")
      filename <- paste(workdir, filename, sep = "")
      print(dim(data.split))
      write.csv(data.split, filename, quote = FALSE, row.names = FALSE)

      workdir <- paste("./split_", i, "/three_prime/counts_files/", sep = "")
      if (!dir.exists(workdir)) dir.create(workdir, recursive = TRUE)
      
      filename <- paste(sample, "mtx.filt", i, "txt", sep = ".")
      filename <- paste(workdir, filename, sep = "")
      write.table(mtx.split, filename, quote = FALSE, row.names = FALSE)
  }
  }
  for (i in 1:length(split_clusters_five)) {
    for (sample in sample.names) {
      data.split <- data.filt.five[which(data.filt.five$three_prime_ID %in% split_clusters_five[[i]]), ]
      split_junc <- unique(as.character(data.split$intron_junction))
      split_junc <- split_junc[split_junc %in% shared_intron_junctions] 
      data.split <- data.split %>% 
        filter(intron_junction %in% split_junc) %>% 
        distinct(intron_junction, .keep_all = TRUE)
      
      # Check if data.split is empty
      if (nrow(data.split) == 0) {
        next  # Skip to the next iteration of the outer loop
      }
      
      mtx.split <- mtx.list[[sample]][split_junc, ]
      workdir <- paste("./split_", i, "/five_prime/data_tables/", sep = "")
      if (!dir.exists(workdir)) dir.create(workdir, recursive = TRUE)
      
      filename <- paste(sample, "data.filt", i, "csv", sep = ".")
      filename <- paste(workdir, filename, sep = "")
      print(dim(data.split))
      
      write.csv(data.split, filename, quote = FALSE, row.names = FALSE)

      workdir <- paste("./split_", i, "/five_prime/counts_files/", sep = "")
      if (!dir.exists(workdir)) dir.create(workdir, recursive = TRUE)
      
      filename <- paste(sample, "mtx.filt", i, "txt", sep = ".")
      filename <- paste(workdir, filename, sep = "")
      write.table(mtx.split, filename, quote = FALSE, row.names = FALSE)
    }
  }
}

message("Done splitting clusters!")
