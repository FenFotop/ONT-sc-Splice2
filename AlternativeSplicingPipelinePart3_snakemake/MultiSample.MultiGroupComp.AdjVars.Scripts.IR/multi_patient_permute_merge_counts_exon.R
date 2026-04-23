## calculate log odds ratio for all IR clusters (exon-based)
## permute by shuffling genotype within each patient
## calculate one odds ratio by pseudobulking across all patients
## adds in regularization factor of 1e-5 to all junctions

library(Matrix)
library(tidyverse)
library(matrixStats)
library(data.table)

args <- commandArgs(TRUE)
options(dplyr.summarise.inform = FALSE)
path.to.split      <- args[1]
path.to.cell.meta  <- args[2]
comp.groups.column <- args[3]
cell.groups.column <- args[4]
nperm              <- as.numeric(args[5])
sample.names       <- args[6]
output.dir         <- args[7]
output.file        <- args[8]
group1_name        <- args[9]
group2_name        <- args[10]

print(args)

path.to.three.matrix <- paste0(path.to.split, "/counts_files")
path.to.three.data   <- paste0(path.to.split, "/data_tables")

sample.names <- unlist(strsplit(sample.names, split = ","))
print(sample.names)

# Read count matrices
setwd(path.to.three.matrix)
file.paths <- list.files(path.to.three.matrix)
suffix <- lapply(strsplit(file.paths[1], split = "mtx"), "[", 2)
files <- unlist(lapply(sample.names, function(s) paste0(s, ".mtx", suffix)))
three.mtx.list <- lapply(files, function(x) read.table(file = x, header = TRUE, stringsAsFactors = FALSE))
names(three.mtx.list) <- sample.names

message("Matrix loaded")

# Load cell metadata
cell.meta.list <- list()
cell.meta.files <- grep(paste0(sample.names, collapse = "|"),
                        list.files(path.to.cell.meta, full.names = TRUE), value = TRUE)
for (file in cell.meta.files) {
  id <- gsub("_metadata.txt", "", basename(file))
  cell.meta <- as.data.frame(read.table(file, stringsAsFactors = FALSE, sep = "\t"))
  cell.meta <- cell.meta[!is.na(cell.meta[, cell.groups.column]), ]
  cell.meta <- cell.meta[, c(cell.groups.column, comp.groups.column)]
  rownames(cell.meta) <- unlist(lapply(strsplit(rownames(cell.meta), split = "_"), "[", 1))
  cell.meta <- cell.meta[rownames(cell.meta) %in% colnames(three.mtx.list[[id]]), ]
  cell.meta.list[[id]] <- cell.meta[!is.na(cell.meta[, comp.groups.column]), ]
}
message("Cell meta data loaded")

setwd(path.to.three.data)
files <- grep(paste0(sample.names, collapse = "|"), list.files(path.to.three.data), value = TRUE)
three.data.comb <- read.csv(files[1])

message("All Data Loaded")

for (sample in sample.names) {
  three.mtx.list[[sample]]$sample <- sample
}

three.prime.clusters <- as.character(unique(three.data.comb$five_prime_ID))

# ---- Vectorized observed logOR ----
DT3 <- as.data.table(three.data.comb)
DT3[, g1all := sum(obs.group1), by = five_prime_ID]
DT3[, g2all := sum(obs.group2), by = five_prime_ID]
DT3[, logOR := log(((obs.group1 + 1e-5) / (g1all - obs.group1 + 1e-5)) /
                   ((obs.group2 + 1e-5) / (g2all - obs.group2 + 1e-5)))]
three.obs.ratio.num <- setNames(DT3$logOR, DT3$exon_coordinates)

three.data.comb$alt_three_prime_exon_coordinates <- paste(three.data.comb$exon_coordinates,
                                                           three.data.comb$five_prime_ID, sep = ":")
three.sample.output <- data.frame(three.obs.logOR.ratio = three.obs.ratio.num,
                                  exon_coordinates = names(three.obs.ratio.num))
three.sample.output <- left_join(three.sample.output, three.data.comb, by = "exon_coordinates")

message("Observed difference calculated")

# Junction-ID-only templates for permutation
three.shf.tmpl <- lapply(sample.names, function(s)
  data.frame(exon_coordinates = three.data.comb$exon_coordinates,
             five_prime_ID    = three.data.comb$five_prime_ID))
names(three.shf.tmpl) <- sample.names

library(dplyr)
library(parallel)
library(pbapply)

# SLURM-aware core selection
total_cores <- length(parallelly::availableWorkers())
slurm_cpus  <- suppressWarnings(as.integer(Sys.getenv("SLURM_CPUS_PER_TASK")))
use_cores   <- if (!is.na(slurm_cpus) && slurm_cpus > 0L) max(1L, slurm_cpus - 1L) else max(1L, total_cores - 1L)
message(sprintf("using %d cores", use_cores))

# ---- Permutation loop (vectorized) ----
results_foreach <- pblapply(1:nperm, function(x) {
  set.seed(x)

  shuffles <- lapply(sample.names, function(s) {
    cv <- as.character(cell.meta.list[[s]][, comp.groups.column])
    sh <- sample(cv, size = length(cv), replace = FALSE)
    names(sh) <- rownames(cell.meta.list[[s]])
    list(g1 = names(sh)[sh == group1_name],
         g2 = names(sh)[sh == group2_name])
  })
  names(shuffles) <- sample.names

  shf3 <- rbindlist(lapply(sample.names, function(s) {
    d <- three.shf.tmpl[[s]]
    m <- three.mtx.list[[s]]
    d$shf.group1 <- rowSums(m[, colnames(m) %in% shuffles[[s]]$g1, drop = FALSE])
    d$shf.group2 <- rowSums(m[, colnames(m) %in% shuffles[[s]]$g2, drop = FALSE])
    d
  }))
  shf3 <- shf3[, .(shf.group1 = sum(shf.group1), shf.group2 = sum(shf.group2)),
               by = .(exon_coordinates, five_prime_ID)]
  shf3[is.na(shf3)] <- 0
  shf3[, g1all := sum(shf.group1), by = five_prime_ID]
  shf3[, g2all := sum(shf.group2), by = five_prime_ID]
  shf3[, logOR := log(((shf.group1 + 1e-5) / (g1all - shf.group1 + 1e-5)) /
                      ((shf.group2 + 1e-5) / (g2all - shf.group2 + 1e-5)))]
  three.shf.ratio.num <- setNames(shf3$logOR, shf3$exon_coordinates)

  abs(three.obs.ratio.num) > abs(three.shf.ratio.num[names(three.obs.ratio.num)])
}, cl = use_cores)

pvals.three <- 1 - colSums(do.call(rbind, results_foreach)) / (nperm + 1)

message("Creating final data frames")
final.three <- data.frame(pvalue = pvals.three, three.obs.logOR.ratio = three.obs.ratio.num,
                          exon_coordinates = names(three.obs.ratio.num))
final.three <- left_join(final.three, three.sample.output)
three.cluster.cov <- final.three %>%
  group_by(five_prime_ID) %>%
  summarise(three.group1.cluster.cov = sum(obs.group1), three.group2.cluster.cov = sum(obs.group2))
final.three <- left_join(final.three, three.cluster.cov, by = "five_prime_ID")

message("Done with permutations!")
message("writing output")
setwd(output.dir)
write.csv(final.three, file = paste0("./alt_three_prime/", output.file, ".csv"), quote = FALSE, row.names = FALSE)
message("Done!!")
