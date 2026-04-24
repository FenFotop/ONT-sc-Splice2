## add in annotations for categories for each junction
## columns being added in are as follows:
## Final_Verdict - Canonical, Cryptic Threeprime, Cryptic Fiveprime, or Cryptic Unanchored
## dPSI_threshold_0 - Cryptic threeprime/ fiveprime + dPSI > 0
## dPSI_threshold_5 - Cryptic threeprime/ fiveprime + dPSI >= 5
## pvalue_threshold - Cryptic threeprime/ fiveprime + pvalue < 0.05
## dPSI_0_pvalue_threshold - Cryptic threeprime/ fiveprime + dPSI > 0, pvalue < 0.05
## dPSI_5_pvalue_threshold - Cryptic threeprime/ fiveprime + dPSI >= 5, pvalue < 0.05

library(tidyverse)

args <- commandArgs(TRUE)
workdir <- args[1]
intronbedfile <- args[2]

path.to.three.data <- paste0(workdir, "/logOR_within_cell_type_ALT_3P_Junctions.txt")
# path.to.five.data = paste0(workdir,"/logOR_within_cell_type_ALT_5P_Junctions.txt")
output <- workdir

## test data
# path.to.three.data = "/gpfs/commons/groups/landau_lab/mariela/U2AF1_ST_project/output/splicing_results/output_dir/genotyped/merge_final_output_IR/logOR_within_cell_type_ALT_3P_Junctions.txt"
# path.to.five.data = "/gpfs/commons/groups/landau_lab/SF3B1_splice_project/7.Permutation_junction_output/7.Comb_patients_ind_celltype_merge_counts_2WT/MDS_P5_P6/NP/logOR_within_cell_type_ALT_5P_Junctions.txt"
# three.path = "/gpfs/commons/groups/landau_lab/mariela/U2AF1_ST_project/output/splicing_results/output_dir/genotyped/merge_final_output_IR/"
three.data <- read.table(path.to.three.data)


## add Final_Verdict column
three.data$Final_Verdict.w.skipping <- ifelse(three.data$IR == "yes", "Intron.Retention", "NA")
table(three.data$Final_Verdict)



three.data <- three.data %>% mutate(
  dPSI_threshold_0 = ifelse(abs(dPSI) > 0, "yes", "no"),
  dPSI_threshold_5 = ifelse(abs(dPSI) >= 5, "yes", "no"),
  pvalue_threshold = ifelse(pvalue < 0.05, "yes", "no")
)


## add in dPSI and pvalue double threshold columns
three.data <- three.data %>% mutate(
  dPSI_0_pvalue_threshold = ifelse(pvalue < 0.05 & abs(dPSI) > 0, "yes", "no"),
  dPSI_5_pvalue_threshold = ifelse(pvalue < 0.05 & abs(dPSI) >= 5, "yes", "no")
)

three.data <- three.data %>%
  filter(Final_Verdict.w.skipping == "Intron.Retention")


introns_bed <- data.table::fread(intronbedfile)
library(GenomicRanges)
colnames(introns_bed) <- c("chr", "start", "end", "gene", "ensembl_gene_id", "strand", "ensembl_transcript_id", "rank", "biotype", "annotation")


introns <- unique(GRanges(introns_bed))

tmp_exon <- gsub("(chr.+):([0-9]+):([0-9]+):(.)", "\\1:\\2-\\3:\\4", three.data$exon_coordinates)


ovps <- findOverlaps(introns, GRanges(tmp_exon), type = "within")

ovps <- as.data.frame(ovps)
ovps$width <- width(introns[ovps$queryHits])

ovps <- ovps %>%
  arrange(subjectHits, -width) %>%
  mutate(subjectHits = as.character(tmp_exon[subjectHits])) %>%
  mutate(queryHits = as.character(introns[queryHits]))

ovps <- ovps[!duplicated(ovps$subjectHits), ]

ovps$queryHits <- gsub("(chr.+):([0-9]+)\\-([0-9]+):(.)", "\\1:\\2:\\3:\\4", ovps$queryHits)
ovps$subjectHits <- gsub("(chr.+):([0-9]+)\\-([0-9]+):(.)", "\\1:\\2:\\3:\\4", ovps$subjectHits)
colnames(ovps) <- c("intron_junction", "exon_coordinates", "width")
three.data <- left_join(three.data, ovps, by = "exon_coordinates")

write.table(three.data, paste0(output, "/logOR_within_cell_type_IR_Junctions_with_threshold_info.txt"))
