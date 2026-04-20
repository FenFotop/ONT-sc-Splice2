library(data.table)
library(optparse)
library(tidyverse)
##Load in the all.exons file from previous annotation step

Sys.setenv("VROOM_CONNECTION_SIZE" = 100000000)
option_parser=OptionParser(
  usage="%prog [options] <name>_all.exons.info.txt path/to/output_folder output_prefix"
)

parsed_args <- parse_args(option_parser, positional_arguments = 4)

input.exon.meta <- parsed_args$args[1]
output_folder <- parsed_args$args[2]
output_prefix <- parsed_args$args[3]
intron_bed <- parsed_args$args[4]
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

exons$IR <- "no"

for (chrom in unique(exons$chr)){
  print(chrom)
  #all.exons.sub <- exons[which(exons$chr == chrom & exons$gene == "RPL18A"),]
  #all.introns.sub <- intron_db[which(intron_db$chr == chrom & intron_db$gene == "RPL18A"),]
  all.exons.sub <- exons[which(exons$chr == chrom),]
  all.introns.sub <- intron_db[which(intron_db$chr == chrom),]
  #all.exons.sub
  #all.exons.sub %>% arrange(start)
  #all.introns.sub %>% arrange(start)
  
  for (row in 1:nrow(all.exons.sub)){
    e.start <- all.exons.sub[row,"start"]
    e.end <- all.exons.sub[row,"end"]
    
    retained.intron <- all.introns.sub[which(all.introns.sub$start > e.start & all.introns.sub$end < e.end), ]
    
    #exon.cluster <- all.exons[which(all.exons$start == e.start & all.exons$end < e.end ), ]
    
    #Here we re trying to figure out if there are 2 separate exons around this intron  - grab all the exons that have same start or same end
    exon.cluster.1 <- all.exons.sub[which(all.exons.sub$start == e.start & all.exons.sub$end %in% retained.intron$start),]
    exon.cluster.2 <- all.exons.sub[which(all.exons.sub$start %in% retained.intron$end & all.exons.sub$end == e.end), ]
    
    if (nrow(retained.intron)>0 & nrow(exon.cluster.1)>0 & nrow(exon.cluster.2)>0){
      all.exons.sub[row,"IR"] <- "yes"
    } 
    
  }
  
  exons[which(exons$chr == chrom),"IR"] <- all.exons.sub$IR
  
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
