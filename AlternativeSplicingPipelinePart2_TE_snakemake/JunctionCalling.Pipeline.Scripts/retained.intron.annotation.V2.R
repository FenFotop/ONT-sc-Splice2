library(data.table)
library(optparse)
library(tidyverse)
##Load in the all.exons file from previous annotation step

Sys.setenv("VROOM_CONNECTION_SIZE" = 100000000)
option_parser=OptionParser(
  usage="%prog [options] <name>_all.exons.info.txt path/to/output_folder output_prefix path/to/intron.bed [path/to/te.bed.gz]"
)

parsed_args <- parse_args(option_parser, positional_arguments = c(4, 5))

input.exon.meta <- parsed_args$args[1]
output_folder   <- parsed_args$args[2]
output_prefix   <- parsed_args$args[3]
intron_bed      <- parsed_args$args[4]
te_bed_file     <- if (length(parsed_args$args) >= 5) parsed_args$args[5] else "NULL"
#patient <- "SMMG001A"

#exons <- read.csv(paste0( "/gpfs/commons/groups/landau_lab/SF3B1_splice_project/21.exon.centric.calling/",patient,"/leafcutter_outputs/exon.meta/",patient,"_all.exons.info.txt"), sep="\t")
#exons <- read.csv(paste0( "/gpfs/commons/groups/landau_lab/rraviram/Suva_lab_GBM/Splicing_ONT/BAM_Paulina/",patient,"/leafcutter_outputs/exon.meta/",patient,"_all.exons.info.txt"), sep="\t")
exons <- unique(read.csv(input.exon.meta, sep="\t"))

##Making an intron database
#paste("zcat < ",all_introns)
intron_db <- fread(intron_bed, data.table = FALSE)
colnames(intron_db)[1:4]=c("chr","start","end","gene")
intron_db$end <- (intron_db$end-1)
#intron_db$chr <- add_chr(intron_db$chr)
print("finished making intron database")

# --- Load TE annotation (optional) ---
# Outputs added to exons:
#   IR_contains_TE : "yes"/"no" — retained intron region overlaps at least one TE
#   IR_TE_names    : "|"-separated TE names found within the retained intron (or ".")
has_te_bed <- !is.null(te_bed_file) && te_bed_file != "NULL" && file.exists(te_bed_file)
if (has_te_bed) {
  print(paste0("Loading TE annotation from: ", te_bed_file))
  te_db <- fread(te_bed_file, data.table = FALSE)
  colnames(te_db)[1:4] <- c("te_chr", "te_start", "te_end", "te_name")
} else {
  te_db <- NULL
  if (te_bed_file != "NULL") warning("TE BED file not found: ", te_bed_file)
}

exons$IR            <- "no"
exons$IR_contains_TE <- "no"
exons$IR_TE_names    <- "."

for (chrom in unique(exons$chr)){
  print(chrom)
  all.exons.sub  <- exons[which(exons$chr == chrom),]
  all.introns.sub <- intron_db[which(intron_db$chr == chrom),]
  te_db_sub <- if (has_te_bed) te_db[te_db$te_chr == chrom, ] else data.frame()

  for (row in 1:nrow(all.exons.sub)){
    e.start <- all.exons.sub[row,"start"]
    e.end   <- all.exons.sub[row,"end"]

    retained.intron <- all.introns.sub[which(all.introns.sub$start > e.start & all.introns.sub$end < e.end), ]

    #Here we re trying to figure out if there are 2 separate exons around this intron  - grab all the exons that have same start or same end
    exon.cluster.1 <- all.exons.sub[which(all.exons.sub$start == e.start & all.exons.sub$end %in% retained.intron$start),]
    exon.cluster.2 <- all.exons.sub[which(all.exons.sub$start %in% retained.intron$end & all.exons.sub$end == e.end), ]

    if (nrow(retained.intron)>0 & nrow(exon.cluster.1)>0 & nrow(exon.cluster.2)>0){
      all.exons.sub[row,"IR"] <- "yes"

      # Check if the retained intron region overlaps any TE
      if (has_te_bed && nrow(te_db_sub) > 0) {
        ir_start <- min(retained.intron$start)
        ir_end   <- max(retained.intron$end)
        te_hits  <- te_db_sub[te_db_sub$te_start < ir_end & te_db_sub$te_end > ir_start, ]
        if (nrow(te_hits) > 0) {
          all.exons.sub[row, "IR_contains_TE"] <- "yes"
          all.exons.sub[row, "IR_TE_names"]    <- paste(unique(te_hits$te_name), collapse = "|")
        }
      }
    }
  }

  exons[which(exons$chr == chrom), "IR"]             <- all.exons.sub$IR
  exons[which(exons$chr == chrom), "IR_contains_TE"] <- all.exons.sub$IR_contains_TE
  exons[which(exons$chr == chrom), "IR_TE_names"]    <- all.exons.sub$IR_TE_names
}

print("Taking strandedness into account")

#Taking strandedness into account 
exons$fivep_distance <- exons$fivep_diff
exons$threep_distance <- exons$threep_diff
#exons = exons %>% filter(strand != 'NA')
exons$exon_coordinates = paste(exons$chr, exons$start, exons$end, exons$strand, sep = ":")
exons$intron_junction <- exons$exon_coordinates

##Add in the exon counts:
#exon.counts <- fread(paste0("/gpfs/commons/groups/landau_lab/SF3B1_splice_project/21.exon.centric.calling/",patient,"/leafcutter_outputs/exon.meta/",patient,"_per.exon_numbers.counts.txt"), sep =" ")
exon.counts <- fread(paste0(output_folder,output_prefix,"_per.exon_numbers.counts.txt"), sep =" ")
exon.counts <- as.data.frame(exon.counts)
rownames(exon.counts) <- exon.counts$exon_coordinates 

#Remove duplicated records 
# Convert to data.table
setDT(exons)
# Sort by exon_coordinates and descending constitutive, then select the first row in each group
exons_result <- exons[order(exon_coordinates, -constitutive), .SD[1], by = exon_coordinates]
exons <- as.data.frame(exons_result)
# Optionally, exons_result is already 'ungrouped' since data.table doesn't explicitly group the data
#exon.counts[exons$exon_coordinates, -1]
print("Final step")
#exons$total.cov <- rowSums(exon.counts[exons$exon_coordinates, -1])
# Compute rowSums in chunks to avoid loading the full subset into memory

coords <- exons$exon_coordinates
n <- length(coords)
chunk_size <- 10000
total_cov <- numeric(n)

for (i in seq(1, n, by = chunk_size)) {
  idx <- i:min(i + chunk_size - 1, n)
  total_cov[idx] <- rowSums(exon.counts[coords[idx], -1, drop = FALSE])
}
exons$total.cov <- total_cov


#write.csv(exons, file=paste0("/path/SF3B1_splice_project/21.exon.centric.calling/",patient,"/leafcutter_outputs/exon.meta/",patient,"_all.exons.info.wRIannotation.csv"))
fwrite(exons, file=paste0(output_folder,output_prefix,"_all.exons.info.wRIannotation.csv"), sep=",")
