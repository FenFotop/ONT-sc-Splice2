## calculate log odds ratio for all clusters with >=2 junctions per cluster and total reads in genotyped cells >= 5
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

path.to.three.matrix <- paste0(path.to.split, "/three_prime/counts_files")
path.to.five.matrix  <- paste0(path.to.split, "/five_prime/counts_files")
path.to.three.data   <- paste0(path.to.split, "/three_prime/data_tables")
path.to.five.data    <- paste0(path.to.split, "/five_prime/data_tables")

sample.names <- unlist(strsplit(sample.names, split = ","))
print(sample.names)

# Read count matrices
setwd(path.to.three.matrix)
file.paths <- list.files(path.to.three.matrix)
suffix <- lapply(strsplit(file.paths[1], split = "mtx"), "[", 2)
files <- unlist(lapply(sample.names, function(s) paste0(s, ".mtx", suffix)))
three.mtx.list <- lapply(files, function(x) read.table(file = x, header = TRUE, stringsAsFactors = FALSE, fill = TRUE))
names(three.mtx.list) <- sample.names

setwd(path.to.five.matrix)
files <- unlist(lapply(sample.names, function(s) paste0(s, ".mtx", suffix)))
five.mtx.list <- lapply(files, function(x) read.table(file = x, header = TRUE, stringsAsFactors = FALSE, fill = TRUE))
names(five.mtx.list) <- sample.names

message("Matrix loaded")

# Load cell metadata
cell.meta.list <- list()
setwd(path.to.cell.meta)
cell.meta.files <- grep(paste0(sample.names, collapse = "|"), list.files(path.to.cell.meta), value = TRUE)
for (file in cell.meta.files) {
  cell.meta <- as.data.frame(read.table(file, stringsAsFactors = FALSE, header = TRUE, fill = TRUE, sep = "\t"))
  id <- gsub("_metadata.txt", "", file)
  cell.meta <- cell.meta[, c(cell.groups.column, comp.groups.column)]
  rownames(cell.meta) <- unlist(lapply(strsplit(rownames(cell.meta), split = "_"), "[", 1))
  cell.meta <- cell.meta[rownames(cell.meta) %in% colnames(three.mtx.list[[id]]), ]
  cell.meta.list[[id]] <- cell.meta[!is.na(cell.meta[, comp.groups.column]), ]
}
message("Cell meta data loaded")

setwd(path.to.three.data)
files <- grep(paste0(sample.names, collapse = "|"), list.files(path.to.three.data), value = TRUE)
three.data.comb <- read.csv(files[1])

setwd(path.to.five.data)
files <- grep(paste0(sample.names, collapse = "|"), list.files(path.to.five.data), value = TRUE)
five.data.comb <- read.csv(files[1])

message("All Data Loaded")

for (sample in sample.names) {
  three.mtx.list[[sample]]$sample <- sample
  five.mtx.list[[sample]]$sample  <- sample
}

three.prime.clusters <- as.character(unique(three.data.comb$five_prime_ID))
five.prime.clusters  <- as.character(unique(five.data.comb$three_prime_ID))

# ---- Vectorized observed logOR (data.table, one pass per cluster direction) ----
DT3 <- as.data.table(three.data.comb)
DT3[, g1all := sum(obs.group1), by = five_prime_ID]
DT3[, g2all := sum(obs.group2), by = five_prime_ID]
DT3[, logOR := log(((obs.group1 + 1e-5) / (g1all - obs.group1 + 1e-5)) /
                   ((obs.group2 + 1e-5) / (g2all - obs.group2 + 1e-5)))]
three.obs.ratio.num <- setNames(DT3$logOR, DT3$intron_junction)

DT5 <- as.data.table(five.data.comb)
DT5[, g1all := sum(obs.group1), by = three_prime_ID]
DT5[, g2all := sum(obs.group2), by = three_prime_ID]
DT5[, logOR := log(((obs.group1 + 1e-5) / (g1all - obs.group1 + 1e-5)) /
                   ((obs.group2 + 1e-5) / (g2all - obs.group2 + 1e-5)))]
five.obs.ratio.num <- setNames(DT5$logOR, DT5$intron_junction)

three.data.comb$alt_three_prime_intron_junction <- paste(three.data.comb$intron_junction, three.data.comb$five_prime_ID, sep = ":")
five.data.comb$alt_five_prime_intron_junction   <- paste(five.data.comb$intron_junction,  five.data.comb$three_prime_ID, sep = ":")

three.sample.output <- data.frame(three.obs.logOR.ratio = three.obs.ratio.num, intron_junction = names(three.obs.ratio.num))
three.sample.output <- left_join(three.sample.output, three.data.comb, by = "intron_junction")

five.sample.output <- data.frame(five.obs.logOR.ratio = five.obs.ratio.num, intron_junction = names(five.obs.ratio.num))
five.sample.output <- left_join(five.sample.output, five.data.comb, by = "intron_junction")

message("Observed difference calculated")

# Junction-ID-only templates for building shuffled data inside the permutation loop
three.shf.tmpl <- lapply(sample.names, function(s)
  data.frame(intron_junction = three.data.comb$intron_junction,
             five_prime_ID   = three.data.comb$five_prime_ID))
names(three.shf.tmpl) <- sample.names

five.shf.tmpl <- lapply(sample.names, function(s)
  data.frame(intron_junction = five.data.comb$intron_junction,
             three_prime_ID  = five.data.comb$three_prime_ID))
names(five.shf.tmpl) <- sample.names

library(dplyr)
library(parallel)
library(pbapply)

# SLURM-aware core selection (same fix applied to Part 2)
total_cores <- length(parallelly::availableWorkers())
slurm_cpus  <- suppressWarnings(as.integer(Sys.getenv("SLURM_CPUS_PER_TASK")))
use_cores   <- if (!is.na(slurm_cpus) && slurm_cpus > 0L) max(1L, slurm_cpus - 1L) else max(1L, total_cores - 1L)
message(sprintf("using %d cores", use_cores))

# ---- Permutation loop (vectorized with data.table inside each iteration) ----
results_foreach <- pblapply(1:nperm, function(x) {
  set.seed(x)

  # Shuffle genotype labels once per sample; same shuffle used for 3' and 5'
  shuffles <- lapply(sample.names, function(s) {
    cv <- as.character(cell.meta.list[[s]][, comp.groups.column])
    sh <- sample(cv, size = length(cv), replace = FALSE)
    names(sh) <- rownames(cell.meta.list[[s]])
    list(g1 = names(sh)[sh == group1_name],
         g2 = names(sh)[sh == group2_name])
  })
  names(shuffles) <- sample.names

  # Three-prime: aggregate shuffled counts and compute logOR in one pass
  shf3 <- rbindlist(lapply(sample.names, function(s) {
    d <- three.shf.tmpl[[s]]
    m <- three.mtx.list[[s]]
    d$shf.group1 <- rowSums(m[, colnames(m) %in% shuffles[[s]]$g1, drop = FALSE])
    d$shf.group2 <- rowSums(m[, colnames(m) %in% shuffles[[s]]$g2, drop = FALSE])
    d
  }))
  shf3 <- shf3[, .(shf.group1 = sum(shf.group1), shf.group2 = sum(shf.group2)),
               by = .(intron_junction, five_prime_ID)]
  shf3[is.na(shf3)] <- 0
  shf3[, g1all := sum(shf.group1), by = five_prime_ID]
  shf3[, g2all := sum(shf.group2), by = five_prime_ID]
  shf3[, logOR := log(((shf.group1 + 1e-5) / (g1all - shf.group1 + 1e-5)) /
                      ((shf.group2 + 1e-5) / (g2all - shf.group2 + 1e-5)))]
  three.shf.ratio.num <- setNames(shf3$logOR, shf3$intron_junction)

  # Five-prime: same pattern
  shf5 <- rbindlist(lapply(sample.names, function(s) {
    d <- five.shf.tmpl[[s]]
    m <- five.mtx.list[[s]]
    d$shf.group1 <- rowSums(m[, colnames(m) %in% shuffles[[s]]$g1, drop = FALSE])
    d$shf.group2 <- rowSums(m[, colnames(m) %in% shuffles[[s]]$g2, drop = FALSE])
    d
  }))
  shf5 <- shf5[, .(shf.group1 = sum(shf.group1), shf.group2 = sum(shf.group2)),
               by = .(intron_junction, three_prime_ID)]
  shf5[is.na(shf5)] <- 0
  shf5[, g1all := sum(shf.group1), by = three_prime_ID]
  shf5[, g2all := sum(shf.group2), by = three_prime_ID]
  shf5[, logOR := log(((shf.group1 + 1e-5) / (g1all - shf.group1 + 1e-5)) /
                      ((shf.group2 + 1e-5) / (g2all - shf.group2 + 1e-5)))]
  five.shf.ratio.num <- setNames(shf5$logOR, shf5$intron_junction)

  three.shf.diff <- abs(three.obs.ratio.num) > abs(three.shf.ratio.num[names(three.obs.ratio.num)])
  five.shf.diff  <- abs(five.obs.ratio.num)  > abs(five.shf.ratio.num[names(five.obs.ratio.num)])
  c(three.shf.diff, five.shf.diff)
}, cl = use_cores)

pvals       <- 1 - colSums(do.call(rbind, results_foreach)) / (nperm + 1)
pvals.three <- pvals[1:length(three.obs.ratio.num)]
pvals.five  <- pvals[(length(three.obs.ratio.num) + 1):length(pvals)]

message("Creating final data frames")
final.three <- data.frame(pvalue = pvals.three, three.obs.logOR.ratio = three.obs.ratio.num,
                          intron_junction = names(three.obs.ratio.num))
final.three <- left_join(final.three, three.sample.output)
three.cluster.cov <- final.three %>%
  group_by(five_prime_ID) %>%
  summarise(three.group1.cluster.cov = sum(obs.group1), three.group2.cluster.cov = sum(obs.group2))
final.three <- left_join(final.three, three.cluster.cov, by = "five_prime_ID")

final.five <- data.frame(pvalue = pvals.five, five.obs.logOR.ratio = five.obs.ratio.num,
                         intron_junction = names(five.obs.ratio.num))
final.five <- left_join(final.five, five.sample.output)
five.cluster.cov <- final.five %>%
  group_by(three_prime_ID) %>%
  summarise(five.group1.cluster.cov = sum(obs.group1), five.group2.cluster.cov = sum(obs.group2))
final.five <- left_join(final.five, five.cluster.cov, by = "three_prime_ID")

message("Done with permutations!")
message("writing output")
setwd(output.dir)
write.csv(final.three, file = paste0("./alt_three_prime/", output.file, ".csv"), quote = FALSE, row.names = FALSE)
write.csv(final.five,  file = paste0("./alt_five_prime/",  output.file, ".csv"), quote = FALSE, row.names = FALSE)
message("Done!!")
