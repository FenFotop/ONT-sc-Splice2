require(tidyverse)
require(magrittr)
require(irlba)
library(optparse)
library(dplyr)
library(data.table)
library(stringr)

Sys.setenv("VROOM_CONNECTION_SIZE" = 100000000)


option_parser=OptionParser(
  usage="%prog [options] <name>_exons.tsv.gz path/to/output_folder output_prefix path/to/annotation_code/prefix_"
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
saveRDS(dat, version=2, paste0(output,"_counts_sc.rda"))
#print("Done")

#dat = readRDS("/gpfs/commons/home/pchamely/leafcutter_scripts/p1_output/p1_ont_v2_align_test/p1_align_smartseq_counts_sc.rda")

#junc_meta = dat %>% select("chrom", "strand", "start", "end")
junc_meta = dat[,1:4]
counts = dat[,5:ncol(dat)]

#Adding on the junction information to the counts matrix!
junc_meta$exon_coordinates <- paste(junc_meta$chrom, ":", junc_meta$start, ":", junc_meta$end,":",junc_meta$strand, sep = "")
saveRDS(junc_meta, version=2, file=paste0(output,"_junc.meta.rda"))

#Save count matrix
perind_numbers.counts <- bind_cols("exon_coordinates" = junc_meta[,c("exon_coordinates")], counts)
print("Writing out full count matrix table")
fwrite(perind_numbers.counts, file= paste0(output ,"_per.exon_numbers.counts.txt") , row.names = F, col.names = T, quote = F, sep = " " )
print("Done")

#perind_numbers.counts[duplicated(perind_numers.counts$exon_junction),] %>% arrange(exon_junction)


## Junction Annotation:
print("Beginning Junction Annotation Process")

# annotation
#exon_file <- paste0(annotation_code, "_all_exons.txt.gz")
all_exons <- paste0(annotation_code,"_all_exons.txt.gz" )
threeprime_file <- paste0(annotation_code,"_threeprime.bed.gz")
fiveprime_file <- paste0(annotation_code,"_fiveprime.bed.gz")

cat("Using annotation at:", annotation_code,"\n")

##This is where it reads in the count matrix with all the exons and cells
cat("Loading counts from", paste0(output ,"_perind_numers.counts.txt"), "\n")
#counts <- read.table(counts_file, check.names=FALSE)
counts <- perind_numbers.counts
counts <- as.data.frame(counts)
rownames(counts) <- counts$exon_coordinates
print("finished loading count matirx")

#For the test run
add_chr=function(chrs)
  if (!grepl("chr",chrs[1])) paste0("chr",chrs) else chrs

all.exons <- as.data.frame(str_split_fixed(rownames(counts), ":", 4), stringsAsFactors = FALSE )
names(all.exons) <- c("chr","start","end","strand")


##--------added in here------##
##Adjusting the meta-data taking strandedness into account
all.exons = all.exons %>% mutate(five_prime = ifelse(strand =="+", start, end),
                                     three_prime = ifelse(strand =="+", end, start))

#Five prime groups
#fp_groups = all.exons %>% group_indices(chr, five_prime)
#all.exons$clusterID_5p <- paste("clu_" ,fp_groups, sep="")

#Three prime groups
#tp_groups = all.exons %>% group_indices(chr, three_prime)
#all.exons$clusterID_3p <- paste("clu_" ,tp_groups, sep="")


setDT(all.exons)
all.exons[, clusterID_5p := paste0("clu_", .GRP), by = .(chr, five_prime)]
# Group by 'chr' and 'three_prime' using .GRP (equivalent to dplyr::group_indices)
all.exons[, clusterID_3p := paste0("clu_", .GRP), by = .(chr, three_prime)]
all.exons <- as.data.frame(all.exons)

##--------added in here------##


#all.exons$chr <- add_chr(all.exons$chr)
if(nrow(all.exons) == 0 ){
  stop("No exon regions found! Please check your input files.")
}

all.exons$start <- as.numeric(all.exons$start)
all.exons$end <- as.numeric(all.exons$end)
all.junctions <- dplyr::select(all.exons, chr, start, end, clusterID_5p)
colnames(all.junctions)[4] <- "clusterID"

##Making exon database
exon_db <- fread(all_exons, data.table = FALSE)
colnames(exon_db)[1:5]=c("chr","start","end","strand","gene")
exon_db$start <- (exon_db$start-1)
exon_db$chr <- add_chr(exon_db$chr)
print("finished making exon database")
#print(exon_db[1:10,])


##Making the 3' database of intron ends!! - note the 3p end of an intron is the 5p (start) of an exon
threeprime_db <- fread(threeprime_file, data.table = FALSE)
colnames(threeprime_db)[1:7]=c("chr","start","end","gene","gene_id","strand","transcript")
threeprime_db$end <- (threeprime_db$end-1)
threeprime_db$start <- (threeprime_db$start-1)
threeprime_db$chr <- add_chr(threeprime_db$chr)
print("finished making 3P database")
#print(threeprime_db[1:10,])

##Making the 5' database of intron ends!! - note the 5p end of an intron is the 3p (end) of an exon
#fiveprime_db <- fread(paste("zcat < ",fiveprime_file), data.table = FALSE)
fiveprime_db <- fread(fiveprime_file, data.table = FALSE)
colnames(fiveprime_db)[1:7]=c("chr","start","end","gene","gene_id","strand","transcript")
fiveprime_db$chr <- add_chr(fiveprime_db$chr)
print("finished making 5p database")
#print(fiveprime_db[1:10,])

##Adjusting the all_exons matrix
all.exons$start_match <- NA
all.exons$end_match <- NA

print("calculate unique matches")
##try putting this before and if it works then see if you can create a new row with the adjusted start and end positions of the exons!!!
# Convert to data.table
setDT(all.exons)
setDT(fiveprime_db)
setDT(threeprime_db)

# Initialize start_match and end_match as numeric columns
##Adjusting the all_introns matrix
all.exons$start_match <- NA
all.exons$end_match <- NA

all.exons[, `:=` (start_match = as.integer(start_match), end_match = as.integer(end_match))]

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
for (chrom in unique(all.exons$chr)) {
  # Subset by chromosome
  exons_chr <- all.exons[chr == chrom]
  five_chr <- fiveprime_db[chr == chrom, .(start)]
  three_chr <- threeprime_db[chr == chrom, .(start)]
  
  # Find closest matches
  start_matches <- find_closest(exons_chr$start, five_chr$start)
  end_matches <- find_closest(exons_chr$end, three_chr$start)
  
  # Update the main data.table with matched values
  all.exons[chr == chrom, start_match := start_matches]
  all.exons[chr == chrom, end_match := end_matches]
}

print("successfull unique matches")

# Calculate differences
all.exons[, `:=` (start_diff = start - start_match, end_diff = end - end_match)]

# Strand adjustment
all.exons[, `:=` (fivep_diff = NA_real_, threep_diff = NA_real_)]
all.exons[strand == "+", `:=` (fivep_diff = start_diff, threep_diff = end_diff)]
all.exons[strand == "-", `:=` (fivep_diff = -end_diff, threep_diff = -start_diff)]

##Create all the intersection databases
all.exons_intersect = all.junctions %>%
  collapse::join(exon_db, how= "left", on=c("chr","start","end"))

#all.exons_intersect[which(all.exons_intersect$gene == "DYNLL1"),]

threeprime_intersect = all.junctions %>%
  select(chr, clusterID, start=end) %>%
  collapse::join(fiveprime_db, how= "left", on=c("chr","start"))

fiveprime_intersect =  all.junctions %>%
  select(chr, clusterID, start) %>%
  collapse::join(threeprime_db, how= "left", on=c("chr","start"))

save(all.exons, all.junctions, all.exons_intersect,threeprime_intersect, fiveprime_intersect, version=2, file=paste0(output,"_junc.annotation.intermediates.Rdata"))

print("Annotating exons")

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
setDT(all.exons_intersect)

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
  bothSSClu <- all.exons_intersect[clusterID == clu]
  
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
  
  # Loop over exons to fill in details
  for (exon in seq_len(nrow(cluster))) {
    cur_start    <- cluster$start[exon]
    cur_end      <- cluster$end[exon]
    fprime_exon  <- fprime[start == cur_start]
    tprime_exon  <- tprime[end == cur_end]
    bothSS_exon  <- bothSSClu[start == cur_start & end == cur_end]
    
    transcripts[[exon]] <- intersect(tprime_exon$transcript, fprime_exon$transcript)
    
    verdict[exon] <- "error"
    
    unknown_3p <- all(is.na(tprime_exon$gene))
    unknown_5p <- all(is.na(fprime_exon$gene))
 

    if (is.na(gene_strand)) {
      verdict[exon] <- "unknown_strand"
    } else {
      if( all( is.na(tprime_exon$gene )) & all( is.na(fprime_exon$gene))){
        verdict[exon] <- "cryptic_unanchored"
      }
      if( (all( is.na(tprime_exon$gene )) & all( !is.na(fprime_exon$gene ) ) & all(gene_strand == "+") ) |
          ( all( is.na(fprime_exon$gene )) & all( !is.na(tprime_exon$gene ) ) & all(gene_strand == "-") )
      ){ verdict[exon] <- "cryptic_threeprime"
      }
      if(
        ( all( !is.na(tprime_exon$gene )) & all( is.na(fprime_exon$gene ) ) & all(gene_strand == "+") ) |
        ( all( !is.na(fprime_exon$gene )) & all( is.na(tprime_exon$gene ) ) & all(gene_strand == "-") )
      ){ verdict[exon] <- "cryptic_fiveprime"
      }
      if( is.na(gene_strand) & ( all( !is.na(tprime_exon$gene )) | all( !is.na(fprime_exon$gene ) ) ) ){
        verdict[exon] <- "cryptic"
      }
      if( # if both splice sites are annotated
        all( !is.na(tprime_exon$gene ) ) & all( !is.na(fprime_exon$gene ) )
      ){
        # test if the splice sites are paired in a known exon
        if( all( !is.na(bothSS_exon$gene )) ){
          verdict[exon] <- "annotated"
        }else{ # both are annotated but never in the same junction
          verdict[exon] <- "novel annotated pair"
        }
        
      }
    }
  }
  
  if( exon == nrow(cluster)){ # only on final exon
      all_transcripts <- unique( unlist( transcripts ) )
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
  } else{
      list(
      verdict = verdict,
      coord = coord,
      gene = gene,
      strand = strand,
      ensemblID = ensemblID,
      transcripts = NA, 
      constitutive = NA
      )
  } 
 
}

stopCluster(cl)

results <- rbindlist(results, fill = T)

setDT(all.exons)
# Ensure 'coord' exists in both tables
all.exons[, coord := sprintf("%s %s %s", chr, start, end)]
# Set 'coord' as the key in both tables for efficient joining
setkey(all.exons, coord)
setkey(results, coord)

# Perform the merge (this is efficient due to indexing)
all.exons <- all.exons[results, nomatch = 0]  # nomatch=0 excludes non-matching rows
all.exons[, constitutive.score := signif(constitutive, digits = 2)]


# Final adjustments
all.exons[, gene := ifelse(is.na(gene), ".", gene)]
all.exons[, ensemblID := ifelse(is.na(ensemblID), ".", ensemblID)]
all.exons[, transcripts := ifelse(transcripts == "", ".", transcripts)]
all.exons[, constitutive.score := ifelse(is.na(constitutive.score), ".", constitutive.score)]

print("Summary Counts for each junc type:")
table(all.exons$verdict)
print("Total number of exons pre-filtering:")
nrow(all.exons)

print("Saving exons info")

all_exons_meta <- as.data.frame(all.exons)

save(all_exons_meta,version=2, file = paste0(output,"_all.exons.info.Rdata") )
write.table(all_exons_meta, file= paste0(output, "_all.exons.info.txt"), quote=F, sep="\t")

print("Done")
