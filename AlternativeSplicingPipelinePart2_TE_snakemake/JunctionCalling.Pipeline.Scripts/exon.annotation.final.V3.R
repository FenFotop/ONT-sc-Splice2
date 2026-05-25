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


#input_matrix <- paste0("/gpfs/commons/groups/landau_lab/SF3B1_splice_project/21.exon.centric.calling/",patient,"/leafcutter_outputs/",patient,"_exons.tsv")
#output_folder <- paste0("/gpfs/commons/groups/landau_lab/SF3B1_splice_project/21.exon.centric.calling/",patient,"/leafcutter_outputs/exon.meta")
#output_prefix <- "A10_S106"
#annotation_code <- "/gpfs/commons/groups/landau_lab/SF3B1_splice_project/ONT_Splice_Pipeline/bin/annotation_reference/leafviz"

cat("Results to be saved in:",output_folder, "\n")

output <- paste0(output_folder,"/",output_prefix)

dat = bigreadr::big_fread2(input_matrix)
print("Saving matrix as an R object")
#saveRDS(dat, version=2, paste0(output,"_counts_sc.rda"))
print("Done")

#dat = readRDS(paste0("/gpfs/commons/groups/landau_lab/SF3B1_splice_project/21.exon.centric.calling/",patient,"/output_files/leafcutter_outputs/exon.meta/",patient,"_counts_sc.rda"))

#junc_meta = dat %>% select("chrom", "strand", "start", "end")
junc_meta = dat[,1:4]
counts = dat[,5:ncol(dat)]

#Adding on the junction information to the counts matrix!
junc_meta$exon_coordinates <- paste(junc_meta$chrom, ":", junc_meta$start, ":", junc_meta$end,":",junc_meta$strand, sep = "")
saveRDS(junc_meta, version=2, file=paste0(output,"_junc.meta.rda"))

#Save count matrix
per.exon_numbers.counts <- bind_cols(junc_meta[,c("exon_coordinates")], counts)
print("Writing out full count matrix table")
write.table(per.exon_numbers.counts, file= paste0(output ,"_per.exon_numbers.counts.txt") , row.names = F, col.names = T, quote = F, sep = " " )
print("Done")

#perind_numbers.counts[duplicated(perind_numers.counts$exon_junction),] %>% arrange(exon_junction)
#print(per.exon_numbers.counts[1:10,1:5])

## Exon Annotation:
print("Beginning Exon Annotation Process")

# annotation
exon_file <- paste0(annotation_code, "_all_exons.txt.gz")
#all_introns <- paste0(annotation_code,"_all_introns.bed.gz" )
threeprime_file <- paste0(annotation_code,"_threeprime.bed.gz")
fiveprime_file <- paste0(annotation_code,"_fiveprime.bed.gz")

cat("Using annotation at:", annotation_code,"\n")

##This is where it reads in the count matrix with all the junctions and cells
cat("Loading counts from", paste0(output ,"_per.exon_numers.counts.txt"), "\n")
#counts <- read.table(counts_file, check.names=FALSE)
counts <- per.exon_numbers.counts
counts <- as.data.frame(counts)
rownames(counts) <- counts$exon_coordinates
print("finished loading count matrix")

#For the test run
add_chr=function(chrs)
  if (!grepl("chr",chrs[1])) paste0("chr",chrs) else chrs

all.exons <- as.data.frame(str_split_fixed(rownames(counts), ":", 4), stringsAsFactors = FALSE )
names(all.exons) <- c("chr","start","end","strand")

#print(all.exons[1:10,])
##--------added in here------##
##Adjusting the meta-data taking strandedness into account
all.exons = all.exons %>% mutate(five_prime = ifelse(strand =="+", start, end),
                                     three_prime = ifelse(strand =="+", end, start))

#Five prime groups
fp_groups = all.exons %>% group_indices(chr, five_prime)
all.exons$clusterID_5p <- paste("clu_" ,fp_groups, sep="")

#Three prime groups
tp_groups = all.exons %>% group_indices(chr, three_prime)
all.exons$clusterID_3p <- paste("clu_" ,tp_groups, sep="")

##--------added in here------##


#all.introns$chr <- add_chr(all.introns$chr)
if(nrow(all.exons) == 0 ){
  stop("No exon regions found! Please check your input files.")
}

all.exons$start <- as.numeric(all.exons$start)
all.exons$end <- as.numeric(all.exons$end)
all.junctions <- dplyr::select(all.exons, chr, start, end, clusterID_5p) ##this only looks at the clustering on one end
colnames(all.junctions)[4] <- "clusterID"

#print(all.exons[1:10,])


##Making an intron database
#paste("zcat < ",all_introns)
#intron_db <- fread(all_introns, data.table = FALSE)
#colnames(intron_db)[1:4]=c("chr","start","end","gene")
#intron_db$end <- (intron_db$end-1)
#intron_db$chr <- add_chr(intron_db$chr)
#print("finished making intron database")

##Making exon database
exon_db <- fread(exon_file, data.table = FALSE)
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

##Adjusting the all_exon matrix
all.exons$start_match <- NA
all.exons$end_match <- NA

##try putting this before and if it works then see if you can create a new row with the adjusted start and end positions of the junctions!!!
for (chrom in unique(all.exons$chr)){

  all.exons[which(all.exons$chr == chrom),]$end_match <- sapply(all.exons[which(all.exons$chr == chrom),]$end, function(x) fiveprime_db[which(fiveprime_db$chr == chrom),]$start[order(abs(x - fiveprime_db[which(fiveprime_db$chr == chrom),]$start))][1])
  all.exons[which(all.exons$chr == chrom),]$start_match <- sapply(all.exons[which(all.exons$chr == chrom),]$start, function(x) threeprime_db[which(threeprime_db$chr == chrom),]$start[order(abs(x - threeprime_db[which(threeprime_db$chr == chrom),]$start))][1])

}
print("finished prep step 1...")
#print(all.exons[1:10,])

all.exons$start_diff <- all.exons$start - all.exons$start_match
all.exons$end_diff <- all.exons$end - all.exons$end_match

##Strand adjust the 5p and 3p distances:
all.exons$fivep_diff <- NA
all.exons$threep_diff <- NA

all.exons[which(all.exons$strand == "+"), ]$fivep_diff <- all.exons[which(all.exons$strand == "+"), ]$start_diff
all.exons[which(all.exons$strand == "+"), ]$threep_diff <- all.exons[which(all.exons$strand == "+"), ]$end_diff

all.exons[which(all.exons$strand == "-"),]$fivep_diff <- (all.exons[which(all.exons$strand == "-"), ]$end_diff)*-1
all.exons[which(all.exons$strand == "-"),]$threep_diff <- (all.exons[which(all.exons$strand == "-"), ]$start_diff)*-1

print("finished prep step 2...")

##Create all the intersection databases
all.exons_intersect = all.junctions %>%
  left_join(exon_db, by=c("chr","start","end"))

#print(all.exons_intersect[which(all.exons_intersect$gene %in%  c("ERGIC3","MPO") ),])

threeprime_intersect = all.junctions %>%
  select(chr, clusterID, start=end) %>%
  left_join(fiveprime_db, by=c("chr","start"))

fiveprime_intersect =  all.junctions %>%
  select(chr, clusterID, start) %>%
  left_join(threeprime_db, by=c("chr","start"))

print("finished prep step 3...")
save(all.exons,all.exons_intersect,threeprime_intersect, fiveprime_intersect, version=2, file=paste0(output,"_junc.annotation.intermediates.Rdata"))


print("Annotating junctions")

verdict.list <- list()
strand.list <- list()
coord.list <- list()
gene.list <- list()
ensemblID.list <- list()
transcripts.list <- list()
constitutive.list <- list()

clusters <- unique(all.junctions$clusterID)

for( clu in clusters){
  # for each exon in the cluster, check for coverage of both
  # output a vector of string descriptions
  cluster <- all.junctions %>% filter( clusterID == clu )
  
  # first subset the intersected files to speed up later query - this uses the data.tables method
  fprimeClu <- fiveprime_intersect %>% filter( clusterID == clu )
  tprimeClu <- threeprime_intersect %>% filter( clusterID == clu )
  bothSSClu <- all.exons_intersect %>% filter( clusterID == clu )
  
  # for each exon in the cluster:
  # create vector of overlapping splice sites, indexed by the row of the intersect
  
  # five prime splice sites
  fprime=cluster %>% left_join(fprimeClu, by=c("chr","start"="start"))
  
  # three prime splice sites
  tprime=cluster %>% left_join(tprimeClu, by=c("chr"="chr","end"="start"))
  
  # both splice sites
  bothSS=cluster %>% left_join(bothSSClu, by=c("chr","start","end"))
  
  # find gene and ensemblID by the most represented gene among all the splice sites - lazy
  cluster_gene <- names(sort(table(c(tprime$gene,fprime$gene)), decreasing = TRUE ))[1]
  
  # if no cluster gene found then leave as "."
  if( is.null(cluster_gene) ){
    cluster_gene <- "."
  }
  
  gene_strand <- NA
  if( cluster_gene != "." ){
    # get strand the same way - would prefer to use the strand of the junction
    strands <- c(tprime$strand, fprime$strand)
    # hope that all junctions align to the same gene on the same strand
    gene_strand <- unique( strands[ strands != "." & !is.na(strands) ])
    if( all(is.na(gene_strand)) | length(gene_strand) != 1 ){
      gene_strand <- NA
    }
  }
  # do the same for EnsemblID
  cluster_ensemblIDs <- names(sort(table( c(tprime$gene_id,fprime$gene_id)), decreasing = TRUE ))
  cluster_ensemblID <- cluster_ensemblIDs[ cluster_ensemblIDs != "." ][1]
  
  if( length( cluster_ensemblID ) == 0 ){
    cluster_ensemblID == "."
  }
  
  verdict <- c()
  coord <- c()
  gene <- c()
  strand <- c()
  ensemblID <- c()
  transcripts <- list()
  
  for( exon in 1:nrow(cluster) ){
    
    coord[exon] <- paste(cluster[exon,]$chr,cluster[exon,]$start, cluster[exon,]$end )
    strand[exon] <- gene_strand
    gene[exon] <- cluster_gene
    ensemblID[exon] <- cluster_ensemblID
    
    fprime_exon=cluster[exon,] %>% left_join(fprime, by=c("chr","start"))
    tprime_exon=cluster[exon,] %>% left_join(tprime, by=c("chr","end"))
    bothSS_exon=cluster[exon,] %>% left_join(bothSSClu, by=c("chr","start","end"))
    
    # for each exon create vector of all transcripts that contain both splice sites
    transcripts[[exon]] <- intersect( tprime_exon$transcript,fprime_exon$transcript )
    
    verdict[exon] <- "error"
    
    unknown_3p=all( is.na(tprime_exon$gene) )
    unknown_5p=all( is.na(fprime_exon$gene) )
    
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
    
    verdict.list[[clu]] <- verdict
    coord.list[[clu]] <- coord
    gene.list[[clu]] <- gene
    strand.list[[clu]] <- strand
    ensemblID.list[[clu]] <- ensemblID
    
    # once all the transcripts for all the exons are found, go back and work out how many constitutive each junction is. Does the junction appear in every transcript?
    
    if( exon == nrow(cluster)){ # only on final exon
      all_transcripts <- unique( unlist( transcripts ) )
      # remove "." - non-existent transcripts
      all_transcripts <- all_transcripts[ all_transcripts != "." ]
      
      constitutive <- lapply( transcripts, FUN = function(x) {
        # for each exon how many transcripts is it seen in?
        x <- x[ x != "." ]
        length(x) / length( all_transcripts)
        
      })
      
      constitutive.list[[clu]] <- constitutive
      
      # collapse all.exons transcripts for each exon into a single string
      transcripts.list[[clu]] <- lapply(transcripts, FUN = function(x) paste( x, collapse = "+" ) )
      
    }
    
  }
  
}

print("Preparing results")

# match the lists together
all.exons$strand.2 <- unlist(strand.list)[ match( paste(all.exons$chr, all.exons$start, all.exons$end ), unlist(coord.list)) ]

all.exons$verdict <- unlist(verdict.list)[ match( paste(all.exons$chr, all.exons$start, all.exons$end ), unlist(coord.list)) ]

all.exons$gene <- unlist(gene.list)[ match( paste( all.exons$chr, all.exons$start, all.exons$end ), unlist(coord.list)) ]

all.exons$ensemblID <- unlist(ensemblID.list)[ match( paste( all.exons$chr, all.exons$start, all.exons$end ), unlist(coord.list)) ]

all.exons$transcripts <- unlist( transcripts.list )[ match( paste( all.exons$chr, all.exons$start, all.exons$end ), unlist(coord.list)) ]

all.exons$constitutive.score <-  unlist( constitutive.list )[ match( paste( all.exons$chr, all.exons$start, all.exons$end ), unlist(coord.list)) ]

#all.exons$prediction <-  unlist( classification.list )[ match( paste( all.exons$chr, all.exons$start, all.exons$end ), unlist(coord.list)) ]

# replace NA values/missing transcripts with "."      ??? SHould i add the start_match and end_match to this list - should not have any missing values

all.exons %<>% mutate( gene=ifelse(is.na(gene), ".", gene),
                         ensemblID=ifelse(is.na(ensemblID), ".", ensemblID),
                         transcripts=ifelse(transcripts == "", ".", transcripts),
                         constitutive.score=signif(constitutive.score, digits = 2))


print("Summary Counts for each junc type:")
table(all.exons$verdict)
print("Total number of junctions pre-filtering:")
nrow(all.exons)

print("Saving exons info")

all_exons_meta <- all.exons

save(all_exons_meta,version=2, file = paste0(output,"_all.exons.info.Rdata") )
write.table(all_exons_meta, file= paste0(output, "_all.exons.info.txt"), quote=F, sep="\t")

print("Done")
