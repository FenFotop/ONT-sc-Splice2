## G-test version of the differential splicing permutation
##
## Statistical approach:
##   - Observed statistic: pooled G statistic across samples (per cluster)
##   - Null: r2dtable generates random contingency tables with the same row/column
##     marginals as observed (Patefield algorithm). Faster and more principled than
##     cell-label shuffling because it permutes at the pseudobulk level.
##   - B permutations suffice (default 3000 vs 100,000 for label-shuffling)
##   - p-value is at the cluster level; all junctions in a cluster share it.
##     Per-junction logOR is still computed for effect-size ranking.
## Output: same column layout as the logOR version — fully compatible with
##   merge_final_output_comb_patient_merge_counts.R downstream.

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
B                  <- as.integer(args[5])   # number of G-test permutations
sample.names       <- args[6]
output.dir         <- args[7]
output.file        <- args[8]
group1_name        <- args[9]
group2_name        <- args[10]

print(args)
message(sprintf("G-test mode: B = %d permutations", B))

path.to.three.matrix <- paste0(path.to.split, "/three_prime/counts_files")
path.to.five.matrix  <- paste0(path.to.split, "/five_prime/counts_files")
path.to.three.data   <- paste0(path.to.split, "/three_prime/data_tables")
path.to.five.data    <- paste0(path.to.split, "/five_prime/data_tables")

sample.names <- unlist(strsplit(sample.names, split = ","))

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

message("Matrices loaded")

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

# Pre-compute per-sample group aggregates ONCE (not repeated per permutation)
three_sample_mats <- lapply(sample.names, function(s) {
  m    <- three.mtx.list[[s]]
  meta <- cell.meta.list[[s]]
  g1   <- rownames(meta)[meta[, comp.groups.column] == group1_name]
  g2   <- rownames(meta)[meta[, comp.groups.column] == group2_name]
  cbind(group1 = rowSums(m[, colnames(m) %in% g1, drop = FALSE]),
        group2 = rowSums(m[, colnames(m) %in% g2, drop = FALSE]))
})
names(three_sample_mats) <- sample.names

five_sample_mats <- lapply(sample.names, function(s) {
  m    <- five.mtx.list[[s]]
  meta <- cell.meta.list[[s]]
  g1   <- rownames(meta)[meta[, comp.groups.column] == group1_name]
  g2   <- rownames(meta)[meta[, comp.groups.column] == group2_name]
  cbind(group1 = rowSums(m[, colnames(m) %in% g1, drop = FALSE]),
        group2 = rowSums(m[, colnames(m) %in% g2, drop = FALSE]))
})
names(five_sample_mats) <- sample.names

setwd(path.to.three.data)
files <- grep(paste0(sample.names, collapse = "|"), list.files(path.to.three.data), value = TRUE)
three.data.comb <- read.csv(files[1])

setwd(path.to.five.data)
files <- grep(paste0(sample.names, collapse = "|"), list.files(path.to.five.data), value = TRUE)
five.data.comb <- read.csv(files[1])

message("All Data Loaded")

three.prime.clusters <- as.character(unique(three.data.comb$five_prime_ID))
five.prime.clusters  <- as.character(unique(five.data.comb$three_prime_ID))

# ---- G-test helpers ----
G_from_M <- function(M) {
  N <- sum(M)
  if (N == 0) return(0)
  rs  <- rowSums(M)
  cs  <- colSums(M)
  E   <- outer(rs, cs) / N
  sel <- M > 0 & E > 0
  if (!any(sel)) return(0)
  2 * sum(M[sel] * log(M[sel] / E[sel]))
}

cluster_gtest <- function(cluster_junctions, sample_mats, B) {
  donor_mats <- lapply(sample_mats, function(mat) {
    sub <- mat[cluster_junctions, , drop = FALSE]
    sub <- sub[rowSums(sub) > 0, , drop = FALSE]
    if (nrow(sub) < 2 || any(colSums(sub) == 0)) return(NULL)
    sub
  })
  donor_mats <- Filter(Negate(is.null), donor_mats)
  if (length(donor_mats) == 0) return(list(stat = NA_real_, p = NA_real_))

  obs <- sum(vapply(donor_mats, G_from_M, numeric(1)))
  if (!is.finite(obs)) return(list(stat = NA_real_, p = NA_real_))

  null_stats <- numeric(B)
  for (b in seq_len(B)) {
    Gb <- 0
    for (M in donor_mats) {
      M0 <- tryCatch(stats::r2dtable(1L, rowSums(M), colSums(M))[[1L]],
                     error = function(e) NULL)
      if (!is.null(M0)) Gb <- Gb + G_from_M(M0)
    }
    null_stats[b] <- Gb
  }

  p <- (1 + sum(null_stats >= obs, na.rm = TRUE)) / (B + 1)
  list(stat = obs, p = p)
}

# ---- Per-junction logOR (effect size, same formula as logOR version) ----
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

message("Per-junction logOR computed (effect size)")

# ---- Run G-test for each cluster ----
message(sprintf("Running G-test for %d 3' clusters and %d 5' clusters",
                length(three.prime.clusters), length(five.prime.clusters)))

three_gtest_res <- rbindlist(lapply(three.prime.clusters, function(cl) {
  juncs <- as.character(three.data.comb[three.data.comb$five_prime_ID == cl, "intron_junction"])
  res   <- cluster_gtest(juncs, three_sample_mats, B)
  data.table(five_prime_ID = cl, gstat = res$stat, pvalue = res$p)
}))

five_gtest_res <- rbindlist(lapply(five.prime.clusters, function(cl) {
  juncs <- as.character(five.data.comb[five.data.comb$three_prime_ID == cl, "intron_junction"])
  res   <- cluster_gtest(juncs, five_sample_mats, B)
  data.table(three_prime_ID = cl, gstat = res$stat, pvalue = res$p)
}))

message("G-tests complete")

# ---- Assemble output (compatible with merge_final_output_comb_patient_merge_counts.R) ----
three.data.comb$alt_three_prime_intron_junction <- paste(three.data.comb$intron_junction, three.data.comb$five_prime_ID, sep = ":")
five.data.comb$alt_five_prime_intron_junction   <- paste(five.data.comb$intron_junction,  five.data.comb$three_prime_ID, sep = ":")

# Attach per-junction logOR and cluster-level p-value
final.three <- as.data.table(three.data.comb)
final.three[, three.obs.logOR.ratio := three.obs.ratio.num[intron_junction]]
final.three <- merge(final.three, three_gtest_res, by = "five_prime_ID")
three.cluster.cov <- final.three[, .(three.group1.cluster.cov = sum(obs.group1),
                                      three.group2.cluster.cov = sum(obs.group2)),
                                   by = five_prime_ID]
final.three <- merge(final.three, three.cluster.cov, by = "five_prime_ID")

final.five <- as.data.table(five.data.comb)
final.five[, five.obs.logOR.ratio := five.obs.ratio.num[intron_junction]]
final.five <- merge(final.five, five_gtest_res, by = "three_prime_ID")
five.cluster.cov <- final.five[, .(five.group1.cluster.cov = sum(obs.group1),
                                     five.group2.cluster.cov = sum(obs.group2)),
                                  by = three_prime_ID]
final.five <- merge(final.five, five.cluster.cov, by = "three_prime_ID")

message("writing output")
setwd(output.dir)
write.csv(as.data.frame(final.three), file = paste0("./alt_three_prime/", output.file, ".csv"), quote = FALSE, row.names = FALSE)
write.csv(as.data.frame(final.five),  file = paste0("./alt_five_prime/",  output.file, ".csv"), quote = FALSE, row.names = FALSE)
message("Done!!")
