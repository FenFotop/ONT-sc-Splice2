## G-test version of the IR differential splicing permutation (exon-based)
##
## Mirrors combined_patient_gtest_ind_celltypes_2WT.R but for intron retention:
##   - Only the 5'-cluster direction (five_prime_ID) is used (IR is donor-anchored)
##   - Uses exon_coordinates instead of intron_junction
##   - r2dtable null preserves per-sample marginals; B permutations suffice
## Output: compatible with merge_combine_output_exon.R downstream.

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
B                  <- as.integer(args[5])
sample.names       <- args[6]
output.dir         <- args[7]
output.file        <- args[8]
group1_name        <- args[9]
group2_name        <- args[10]

print(args)
message(sprintf("G-test IR mode: B = %d permutations", B))

path.to.three.matrix <- paste0(path.to.split, "/counts_files")
path.to.three.data   <- paste0(path.to.split, "/data_tables")

sample.names <- unlist(strsplit(sample.names, split = ","))

# Read count matrices
setwd(path.to.three.matrix)
file.paths <- list.files(path.to.three.matrix)
suffix <- lapply(strsplit(file.paths[1], split = "mtx"), "[", 2)
files <- unlist(lapply(sample.names, function(s) paste0(s, ".mtx", suffix)))
three.mtx.list <- lapply(files, function(x) read.table(file = x, header = TRUE, stringsAsFactors = FALSE))
names(three.mtx.list) <- sample.names

message("Matrices loaded")

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

# Pre-compute per-sample group aggregates ONCE
three_sample_mats <- lapply(sample.names, function(s) {
  m    <- three.mtx.list[[s]]
  meta <- cell.meta.list[[s]]
  g1   <- rownames(meta)[meta[, comp.groups.column] == group1_name]
  g2   <- rownames(meta)[meta[, comp.groups.column] == group2_name]
  cbind(group1 = rowSums(m[, colnames(m) %in% g1, drop = FALSE]),
        group2 = rowSums(m[, colnames(m) %in% g2, drop = FALSE]))
})
names(three_sample_mats) <- sample.names

setwd(path.to.three.data)
files <- grep(paste0(sample.names, collapse = "|"), list.files(path.to.three.data), value = TRUE)
three.data.comb <- read.csv(files[1])

message("All Data Loaded")

three.prime.clusters <- as.character(unique(three.data.comb$five_prime_ID))

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

cluster_gtest <- function(cluster_coords, sample_mats, B) {
  donor_mats <- lapply(sample_mats, function(mat) {
    sub <- mat[cluster_coords, , drop = FALSE]
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

# ---- Per-exon logOR (effect size) ----
DT3 <- as.data.table(three.data.comb)
DT3[, g1all := sum(obs.group1), by = five_prime_ID]
DT3[, g2all := sum(obs.group2), by = five_prime_ID]
DT3[, logOR := log(((obs.group1 + 1e-5) / (g1all - obs.group1 + 1e-5)) /
                   ((obs.group2 + 1e-5) / (g2all - obs.group2 + 1e-5)))]
three.obs.ratio.num <- setNames(DT3$logOR, DT3$exon_coordinates)

message("Per-exon logOR computed (effect size)")

# ---- Run G-test for each cluster ----
message(sprintf("Running G-test for %d IR clusters", length(three.prime.clusters)))

three_gtest_res <- rbindlist(lapply(three.prime.clusters, function(cl) {
  coords <- as.character(three.data.comb[three.data.comb$five_prime_ID == cl, "exon_coordinates"])
  res    <- cluster_gtest(coords, three_sample_mats, B)
  data.table(five_prime_ID = cl, gstat = res$stat, pvalue = res$p)
}))

message("G-tests complete")

# ---- Assemble output (compatible with merge_combine_output_exon.R) ----
three.data.comb$alt_three_prime_exon_coordinates <- paste(three.data.comb$exon_coordinates,
                                                           three.data.comb$five_prime_ID, sep = ":")

final.three <- as.data.table(three.data.comb)
final.three[, three.obs.logOR.ratio := three.obs.ratio.num[exon_coordinates]]
final.three <- merge(final.three, three_gtest_res, by = "five_prime_ID")
three.cluster.cov <- final.three[, .(three.group1.cluster.cov = sum(obs.group1),
                                      three.group2.cluster.cov = sum(obs.group2)),
                                   by = five_prime_ID]
final.three <- merge(final.three, three.cluster.cov, by = "five_prime_ID")

message("writing output")
setwd(output.dir)
write.csv(as.data.frame(final.three),
          file = paste0("./alt_three_prime/", output.file, ".csv"),
          quote = FALSE, row.names = FALSE)
message("Done!!")
