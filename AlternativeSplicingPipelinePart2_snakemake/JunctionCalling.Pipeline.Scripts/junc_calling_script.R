require(tidyverse)
require(magrittr)
require(irlba)
library(optparse)
library(dplyr)
library(data.table)
library(stringr)

Sys.setenv("VROOM_CONNECTION_SIZE" = 100000000)


option_parser=OptionParser(
  usage="%prog [options] <name>_counts_sc_txt.gz path/to/output_folder output_prefix path/to/annotation_code/prefix_"
)

parsed_args <- parse_args(option_parser, positional_arguments = 4)

input_matrix <- parsed_args$args[1]
output_folder <- parsed_args$args[2]
output_prefix <- parsed_args$args[3]
annotation_code <- parsed_args$args[4]

cat("Results to be saved in:",output_folder, "\n")

output <- paste0(output_folder,"/",output_prefix)

dat = bigreadr::big_fread2(input_matrix, data.table = F)
#print("Saving matrix as an R object")
#saveRDS(dat, version=2, paste0(output,"_counts_sc.rda"))
#print("Done")

#dat = readRDS("/gpfs/commons/home/pchamely/leafcutter_scripts/p1_output/p1_ont_v2_align_test/p1_align_smartseq_counts_sc.rda")

#junc_meta = dat %>% select("chrom", "strand", "start", "end")
junc_meta = dat[,1:4]
counts = dat[,5:ncol(dat)]

#Adding on the junction information to the counts matrix!
junc_meta$intron_junction <- paste(junc_meta$chrom, ":", junc_meta$start, ":", junc_meta$end,":",junc_meta$strand, sep = "")

#Save count matrix
perind_numbers.counts <- bind_cols("intron_junction" = junc_meta[,c("intron_junction")], counts)
print("Writing out full count matrix table")
fwrite(perind_numbers.counts, file= paste0(output ,"_perind_numbers.counts.txt") , row.names = F, col.names = T, quote = F, sep = " " )
print("Done")

#perind_numbers.counts[duplicated(perind_numers.counts$intron_junction),] %>% arrange(intron_junction)


## Junction Annotation:
print("Beginning Junction Annotation Process")

# annotation
#exon_file <- paste0(annotation_code, "_all_exons.txt.gz")
all_introns <- paste0(annotation_code,"_all_introns.bed.gz" )
threeprime_file <- paste0(annotation_code,"_threeprime.bed.gz")
fiveprime_file <- paste0(annotation_code,"_fiveprime.bed.gz")

cat("Using annotation at:", annotation_code,"\n")

##This is where it reads in the count matrix with all the junctions and cells
cat("Loading counts from", paste0(output ,"_perind_numers.counts.txt"), "\n")
#counts <- read.table(counts_file, check.names=FALSE)
counts <- perind_numbers.counts
counts <- as.data.frame(counts)
rownames(counts) <- counts$intron_junction
print("finished loading count matirx")

#For the test run
add_chr=function(chrs)
  if (!grepl("chr",chrs[1])) paste0("chr",chrs) else chrs

all.introns <- as.data.frame(str_split_fixed(rownames(counts), ":", 4), stringsAsFactors = FALSE )
names(all.introns) <- c("chr","start","end","strand_1")


##--------added in here------##
##Adjusting the meta-data taking strandedness into account
all.introns = all.introns %>% mutate(five_prime = ifelse(strand_1 =="+", start, end),
                                     three_prime = ifelse(strand_1 =="+", end, start))

#Five prime groups
#fp_groups = all.introns %>% group_indices(chr, five_prime)
#all.introns$clusterID_5p <- paste("clu_" ,fp_groups, sep="")

#Three prime groups
#tp_groups = all.introns %>% group_indices(chr, three_prime)
#all.introns$clusterID_3p <- paste("clu_" ,tp_groups, sep="")


setDT(all.introns)
all.introns[, clusterID_5p := paste0("clu_", .GRP), by = .(chr, five_prime)]
# Group by 'chr' and 'three_prime' using .GRP (equivalent to dplyr::group_indices)
all.introns[, clusterID_3p := paste0("clu_", .GRP), by = .(chr, three_prime)]
all.introns <- as.data.frame(all.introns)

##--------added in here------##


#all.introns$chr <- add_chr(all.introns$chr)
if(nrow(all.introns) == 0 ){
  stop("No intron regions found! Please check your input files.")
}

all.introns$start <- as.numeric(all.introns$start)
all.introns$end <- as.numeric(all.introns$end)
all.junctions <- dplyr::select(all.introns, chr, start, end, clusterID_5p)
colnames(all.junctions)[4] <- "clusterID"

##Making an intron database
#paste("zcat < ",all_introns)
intron_db <- fread(all_introns, data.table = FALSE)
colnames(intron_db)[1:4]=c("chr","start","end","gene")
intron_db$end <- (intron_db$end-1)
intron_db$chr <- add_chr(intron_db$chr)
print("finished making intron database")

##Making the 3' database
threeprime_db <- fread(threeprime_file, data.table = FALSE)
colnames(threeprime_db)[1:7]=c("chr","start","end","gene","gene_id","strand","transcript")
threeprime_db$end <- (threeprime_db$end-1)
threeprime_db$start <- (threeprime_db$start-1)
threeprime_db$chr <- add_chr(threeprime_db$chr)
print("finished making 3P database")

##Making the 5' database
#fiveprime_db <- fread(paste("zcat < ",fiveprime_file), data.table = FALSE)
fiveprime_db <- fread(fiveprime_file, data.table = FALSE)
colnames(fiveprime_db)[1:7]=c("chr","start","end","gene","gene_id","strand","transcript")
fiveprime_db$chr <- add_chr(fiveprime_db$chr)
print("finished making 5' database")

##Adjusting the all_introns matrix
all.introns$start_match <- NA
all.introns$end_match <- NA

print("calculate unique matches")
##try putting this before and if it works then see if you can create a new row with the adjusted start and end positions of the junctions!!!
# Convert to data.table
setDT(all.introns)
setDT(fiveprime_db)
setDT(threeprime_db)

# Initialize start_match and end_match as numeric columns
all.introns$start_match <- NA
all.introns$end_match <- NA

all.introns[, `:=` (start_match = as.integer(start_match), end_match = as.integer(end_match))]

# Function to get the closest match (vectorised binary search — O(m log m + n log n))
find_closest <- function(query, db) {
  if (length(db) == 0) return(rep(NA_integer_, length(query)))
  db_sorted <- sort(db)
  idx <- findInterval(query, db_sorted)
  lo  <- pmax(1L, idx)
  hi  <- pmin(length(db_sorted), idx + 1L)
  ifelse(abs(query - db_sorted[lo]) <= abs(query - db_sorted[hi]),
         db_sorted[lo], db_sorted[hi])
}

# Loop over chromosomes
for (chrom in unique(all.introns$chr)) {
  # Subset by chromosome
  introns_chr <- all.introns[chr == chrom]
  five_chr <- fiveprime_db[chr == chrom, .(start)]
  three_chr <- threeprime_db[chr == chrom, .(start)]
  
  # Find closest matches
  start_matches <- find_closest(introns_chr$start, five_chr$start)
  end_matches <- find_closest(introns_chr$end, three_chr$start)
  
  # Update the main data.table with matched values
  all.introns[chr == chrom, start_match := start_matches]
  all.introns[chr == chrom, end_match := end_matches]
}

print("successfull unique matches")

# Calculate differences
all.introns[, `:=` (start_diff = start - start_match, end_diff = end - end_match)]

# Strand adjustment
all.introns[, `:=` (fivep_diff = NA_real_, threep_diff = NA_real_)]
all.introns[strand_1 == "+", `:=` (fivep_diff = start_diff, threep_diff = end_diff)]
all.introns[strand_1 == "-", `:=` (fivep_diff = -end_diff, threep_diff = -start_diff)]

##Create all the intersection databases
all.introns_intersect = all.junctions %>%
  collapse::join(intron_db, how= "left", on=c("chr","start","end"))

#all.introns_intersect[which(all.introns_intersect$gene == "DYNLL1"),]

threeprime_intersect = all.junctions %>%
  select(chr, clusterID, start=end) %>%
  collapse::join(threeprime_db, how= "left", on=c("chr","start"))

fiveprime_intersect =  all.junctions %>%
  select(chr, clusterID, start) %>%
  collapse::join(fiveprime_db, how= "left", on=c("chr","start"))

print("Annotating junctions")

library(data.table)
library(foreach)
library(doParallel)


# Set up parallel backend — respect SLURM allocation to avoid oversubscription
total_cores <- length(parallelly::availableWorkers())
slurm_cpus  <- suppressWarnings(as.integer(Sys.getenv("SLURM_CPUS_PER_TASK")))
use_cores   <- if (!is.na(slurm_cpus) && slurm_cpus > 0L) max(1L, slurm_cpus - 1L) else max(1L, total_cores - 2L)
print(paste0("Using ", use_cores, " cores"))
cl <- makeCluster(use_cores)
registerDoParallel(cl)

# Convert data frames to data.tables
setDT(all.junctions)
setDT(fiveprime_intersect)
setDT(threeprime_intersect)
setDT(all.introns_intersect)

clusters <- unique( all.junctions$clusterID)

# Initialize lists to store results
verdict.list <- vector("list", length(clusters))
strand.list <- vector("list", length(clusters))
coord.list <- vector("list", length(clusters))
gene.list <- vector("list", length(clusters))
ensemblID.list <- vector("list", length(clusters))
transcripts.list <- vector("list", length(clusters))
constitutive.list <- vector("list", length(clusters))

# Parallelize the cluster loop
results <- foreach(clu = clusters, .packages = c("data.table", "dplyr")) %dopar% {
  # Subset cluster-specific data
  cluster <- all.junctions[clusterID == clu]
  fprimeClu <- fiveprime_intersect[clusterID == clu]
  tprimeClu <- threeprime_intersect[clusterID == clu]
  bothSSClu <- all.introns_intersect[clusterID == clu]
  
  # Join to find matching splice sites
  fprime <- cluster[fprimeClu, on = .(chr, start), allow.cartesian = TRUE]
  tprime <- cluster[tprimeClu, on = .(chr, end = start), allow.cartesian = TRUE]
  bothSS <- cluster[bothSSClu, on = .(chr, start, end), , allow.cartesian = TRUE]
  
  # Find the most represented gene and Ensembl ID
  cluster_gene <- names(sort(table(c(tprime$gene, fprime$gene)), decreasing = TRUE))[1]
  cluster_gene <- ifelse(is.null(cluster_gene), ".", cluster_gene)
  
  gene_strand <- if (cluster_gene != ".") {
    strands <- c(tprime$strand, fprime$strand)
    unique_strand <- unique(strands[!is.na(strands) & strands != "."])
    if (length(unique_strand) == 1) unique_strand else NA
  } else {
    NA
  }
  
  cluster_ensemblIDs <- names(sort(table(c(tprime$gene_id, fprime$gene_id)), decreasing = TRUE))
  cluster_ensemblID <- cluster_ensemblIDs[cluster_ensemblIDs != "."][1]
  cluster_ensemblID <- ifelse(length(cluster_ensemblID) == 0, ".", cluster_ensemblID)
  
  # Initialize vectors for the current cluster
  verdict <- rep("error", nrow(cluster))
  coord <- paste(cluster$chr, cluster$start, cluster$end)
  gene <- rep(cluster_gene, nrow(cluster))
  strand <- rep(gene_strand, nrow(cluster))
  ensemblID <- rep(cluster_ensemblID, nrow(cluster))
  transcripts <- vector("list", nrow(cluster))
  
  # Loop over introns to fill in details
  for (intron in seq_len(nrow(cluster))) {
    cur_start     <- cluster$start[intron]
    cur_end       <- cluster$end[intron]
    fprime_intron <- fprime[start == cur_start]
    tprime_intron <- tprime[end == cur_end]
    bothSS_intron <- bothSSClu[start == cur_start & end == cur_end]
    
    transcripts[[intron]] <- intersect(tprime_intron$transcript, fprime_intron$transcript)
    
    unknown_3p <- all(is.na(tprime_intron$gene))
    unknown_5p <- all(is.na(fprime_intron$gene))
    
    #if (is.na(gene_strand)) {
    #  verdict[intron] <- "unknown_strand"
    #} else {
    #  if (unknown_3p && unknown_5p) {
    #    verdict[intron] <- "cryptic_unanchored"
    #  } else if ((unknown_3p && all(!is.na(fprime_intron$gene)) && gene_strand == "+") ||
    #             (unknown_5p && all(!is.na(tprime_intron$gene)) && gene_strand == "-")) {
    #    verdict[intron] <- "cryptic_threeprime"
    #  } else if ((unknown_5p && all(!is.na(tprime_intron$gene)) && gene_strand == "+") ||
    #             (unknown_3p && all(!is.na(fprime_intron$gene)) && gene_strand == "-")) {
    #    verdict[intron] <- "cryptic_fiveprime"
    #  } else if (!is.na(gene_strand) && (any(!is.na(tprime_intron$gene)) || any(!is.na(fprime_intron$gene)))) {
    #   verdict[intron] <- "cryptic"
    #  } else if (all(!is.na(tprime_intron$gene)) && all(!is.na(fprime_intron$gene))) {
    #    verdict[intron] <- if (all(!is.na(bothSS_intron$gene))) "annotated" else "novel annotated pair"
    #  }
    
    
    if (is.na(gene_strand)) {
      verdict[intron] <- "unknown_strand"
    } else {
      if( all( is.na(tprime_intron$gene )) & all( is.na(fprime_intron$gene))){
        verdict[intron] <- "cryptic_unanchored"
      }
      if( (all( is.na(tprime_intron$gene )) & all( !is.na(fprime_intron$gene ) ) & all(gene_strand == "+") ) |
          ( all( is.na(fprime_intron$gene )) & all( !is.na(tprime_intron$gene ) ) & all(gene_strand == "-") )
      ){ verdict[intron] <- "cryptic_threeprime"
      }
      if(
        ( all( !is.na(tprime_intron$gene )) & all( is.na(fprime_intron$gene ) ) & all(gene_strand == "+") ) |
        ( all( !is.na(fprime_intron$gene )) & all( is.na(tprime_intron$gene ) ) & all(gene_strand == "-") )
      ){ verdict[intron] <- "cryptic_fiveprime"
      }
      if( is.na(gene_strand) & ( all( !is.na(tprime_intron$gene )) | all( !is.na(fprime_intron$gene ) ) ) ){
        verdict[intron] <- "cryptic"
      }
      if( # if both splice sites are annotated
        all( !is.na(tprime_intron$gene ) ) & all( !is.na(fprime_intron$gene ) )
      ){
        # test if the splice sites are paired in a known intron
        if( all( !is.na(bothSS_intron$gene )) ){
          verdict[intron] <- "annotated"
        }else{ # both are annotated but never in the same junction
          verdict[intron] <- "novel annotated pair"
        }
        
      }
    }
  }
  
  
  # Calculate constitutive score
  all_transcripts <- unique(unlist(transcripts))
  all_transcripts <- all_transcripts[all_transcripts != "."]
  constitutive <- sapply(transcripts, function(x) length(x[x != "."]) / length(all_transcripts))
  
  # Return results for the current cluster
  list(
    verdict = verdict,
    coord = coord,
    gene = gene,
    strand = strand,
    ensemblID = ensemblID,
    transcripts = sapply(transcripts, paste, collapse = "+"),
    constitutive = constitutive
  )
}

stopCluster(cl)

results <- rbindlist(results)

setDT(all.introns)
# Ensure 'coord' exists in both tables
all.introns[, coord := sprintf("%s %s %s", chr, start, end)]
# Set 'coord' as the key in both tables for efficient joining
setkey(all.introns, coord)
setkey(results, coord)

# Perform the merge (this is efficient due to indexing)
all.introns <- all.introns[results, nomatch = 0]  # nomatch=0 excludes non-matching rows
all.introns[, constitutive.score := signif(constitutive, digits = 2)]


# Final adjustments
all.introns[, gene := ifelse(is.na(gene), ".", gene)]
all.introns[, ensemblID := ifelse(is.na(ensemblID), ".", ensemblID)]
all.introns[, transcripts := ifelse(transcripts == "", ".", transcripts)]
all.introns[, constitutive.score := ifelse(is.na(constitutive.score), ".", constitutive.score)]

print("Summary Counts for each junc type:")
table(all.introns$verdict)
print("Total number of junctions pre-filtering:")
nrow(all.introns)

print("Saving Introns info")

all_introns_meta <- as.data.frame(all.introns)

save(all_introns_meta,version=2, file = paste0(output,"_all.introns.info.Rdata") )
write.table(all_introns_meta, file= paste0(output, "_all.introns.info.txt"), quote=F, sep="\t")

print("Done")
