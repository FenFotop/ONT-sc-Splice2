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


path.to.three.data <- paste0(workdir, "/logOR_within_cell_type_ALT_3P_Junctions.txt")
path.to.five.data <- paste0(workdir, "/logOR_within_cell_type_ALT_5P_Junctions.txt")
output <- workdir

## test data
# path.to.three.data = "/gpfs/commons/groups/landau_lab/SF3B1_splice_project/7.Permutation_junction_output/7.Comb_patients_ind_celltype_merge_counts_2WT/MDS_P5_P6/NP/logOR_within_cell_type_ALT_3P_Junctions.txt"
# path.to.five.data = "/gpfs/commons/groups/landau_lab/SF3B1_splice_project/7.Permutation_junction_output/7.Comb_patients_ind_celltype_merge_counts_2WT/MDS_P5_P6/NP/logOR_within_cell_type_ALT_5P_Junctions.txt"

three.data <- read.table(path.to.three.data)%>% 
  filter(!is.na(dPSI)) # remove NAs
five.data <- read.table(path.to.five.data) %>% 
  filter(!is.na(dPSI)) # remove NAs

## add Final_Verdict column
three.data$Final_Verdict <- NA
three.data[which(three.data$startClass == "main" & three.data$endClass == "main"), "Final_Verdict"] <- "Canonical"
three.data[which(((three.data$startClass == "main" & three.data$endClass == "not_main_3_prime") | (three.data$startClass == "not_main_3_prime" & three.data$endClass == "main")) &
  (three.data$threep_distance > (-100) & three.data$threep_distance < 0)), "Final_Verdict"] <- "Cryptic_threeprime"
three.data[which(((three.data$startClass == "main" & three.data$endClass == "not_main_5_prime") | (three.data$startClass == "not_main_5_prime" & three.data$endClass == "main")) &
  (three.data$fivep_distance > (-100) & three.data$fivep_distance < 0)), "Final_Verdict"] <- "Cryptic_fiveprime"
three.data[which(three.data$startClass == "not_main_3_prime" & three.data$endClass == "not_main_5_prime"), "Final_Verdict"] <- "cryptic_unanchored"
three.data[which(three.data$startClass == "not_main_5_prime" & three.data$endClass == "not_main_3_prime"), "Final_Verdict"] <- "cryptic_unanchored"
three.data[which(((three.data$startClass == "not_main_3_prime" & three.data$endClass == "main") | (three.data$startClass == "main" & three.data$endClass == "not_main_3_prime")) &
  (three.data$threep_distance < (-100) | three.data$threep_distance > 0)), "Final_Verdict"] <- "Alternative_threeprime"
three.data[which(((three.data$startClass == "not_main_5_prime" & three.data$endClass == "main") | (three.data$startClass == "main" & three.data$endClass == "not_main_5_prime")) &
  (three.data$fivep_distance < (-100) | three.data$fivep_distance > 0)), "Final_Verdict"] <- "Alternative_fiveprime"

table(three.data$Final_Verdict)

three.data$num.skipped.exons <- three.data$relStartExon_skipping - three.data$relEndExon_skipping
three.data$exon.skip <- "no"
three.data[which(abs(three.data$num.skipped.exons) > 0), ]$exon.skip <- "yes"
three.data$Final_Verdict.w.skipping <- paste(three.data$Final_Verdict)
three.data[which(three.data$exon.skip == "yes"), ]$Final_Verdict.w.skipping <- "exon.skip"


five.data$Final_Verdict <- NA
five.data[which(five.data$startClass == "main" & five.data$endClass == "main"), "Final_Verdict"] <- "Canonical"
five.data[which(((five.data$startClass == "main" & five.data$endClass == "not_main_3_prime") | (five.data$startClass == "not_main_3_prime" & five.data$endClass == "main")) &
  (five.data$threep_distance > (-100) & five.data$threep_distance < 0)), "Final_Verdict"] <- "Cryptic_threeprime"
five.data[which(((five.data$startClass == "main" & five.data$endClass == "not_main_5_prime") | (five.data$startClass == "not_main_5_prime" & five.data$endClass == "main")) &
  (five.data$fivep_distance > (-100) & five.data$fivep_distance < 0)), "Final_Verdict"] <- "Cryptic_fiveprime"
five.data[which(five.data$startClass == "not_main_3_prime" & five.data$endClass == "not_main_5_prime"), "Final_Verdict"] <- "cryptic_unanchored"
five.data[which(five.data$startClass == "not_main_5_prime" & five.data$endClass == "not_main_3_prime"), "Final_Verdict"] <- "cryptic_unanchored"
five.data[which(((five.data$startClass == "not_main_3_prime" & five.data$endClass == "main") | (five.data$startClass == "main" & five.data$endClass == "not_main_3_prime")) &
  (five.data$threep_distance < (-100) | five.data$threep_distance > 0)), "Final_Verdict"] <- "Alternative_threeprime"
five.data[which(((five.data$startClass == "not_main_5_prime" & five.data$endClass == "main") | (five.data$startClass == "main" & five.data$endClass == "not_main_5_prime")) &
  (five.data$fivep_distance < (-100) | five.data$fivep_distance > 0)), "Final_Verdict"] <- "Alternative_fiveprime"
table(five.data$Final_Verdict)

five.data$num.skipped.exons <- five.data$relStartExon_skipping - five.data$relEndExon_skipping
five.data$exon.skip <- "no"
five.data[which(abs(five.data$num.skipped.exons) > 0), ]$exon.skip <- "yes"
five.data$Final_Verdict.w.skipping <- paste(five.data$Final_Verdict)
five.data[which(five.data$exon.skip == "yes"), ]$Final_Verdict.w.skipping <- "exon.skip"

## add in dPSI and pvalue single threshold columns
three.data$dPSI_threshold_0 <- NA
five.data$dPSI_threshold_0 <- NA
three.data$dPSI_threshold_5 <- NA
five.data$dPSI_threshold_5 <- NA
three.data$pvalue_threshold <- NA
five.data$pvalue_threshold <- NA
three.data <- three.data %>% mutate(
  dPSI_threshold_0 = ifelse(abs(dPSI) > 0, "yes", "no"),
  dPSI_threshold_5 = ifelse(abs(dPSI) >= 5, "yes", "no"),
  pvalue_threshold = ifelse(pvalue < 0.05, "yes", "no")
)
five.data <- five.data %>% mutate(
  dPSI_threshold_0 = ifelse(abs(dPSI) > 0, "yes", "no"),
  dPSI_threshold_5 = ifelse(abs(dPSI) >= 5, "yes", "no"),
  pvalue_threshold = ifelse(pvalue < 0.05, "yes", "no")
)


## add in dPSI and pvalue double threshold columns
three.data$dPSI_0_pvalue_threshold <- "no"
five.data$dPSI_0_pvalue_threshold <- "no"
three.data$dPSI_5_pvalue_threshold <- "no"
five.data$dPSI_5_pvalue_threshold <- "no"
if (nrow(three.data[three.data$pvalue < 0.05 & abs(three.data$dPSI) > 0, ]) > 0) {
  three.data[three.data$pvalue < 0.05 & abs(three.data$dPSI) > 0, "dPSI_0_pvalue_threshold"] <- "yes"
}

if (nrow(three.data[three.data$pvalue < 0.05 & abs(three.data$dPSI) >= 5, ]) > 0) {
  three.data[three.data$pvalue < 0.05 & abs(three.data$dPSI) >= 5, "dPSI_5_pvalue_threshold"] <- "yes"
}

if (nrow(five.data[five.data$pvalue < 0.05 & abs(five.data$dPSI) > 0, ]) > 0) {
  five.data[five.data$pvalue < 0.05 & abs(five.data$dPSI) > 0, "dPSI_0_pvalue_threshold"] <- "yes"
}

if (nrow(five.data[five.data$pvalue < 0.05 & abs(five.data$dPSI) >= 5, ]) > 0) {
  five.data[five.data$pvalue < 0.05 & abs(five.data$dPSI) >= 5, "dPSI_5_pvalue_threshold"] <- "yes"
}


## filter for only cryptic sites
three.data.cryptic <- three.data %>% filter(Final_Verdict == "Cryptic_threeprime")
five.data.cryptic <- five.data %>% filter(Final_Verdict == "Cryptic_fiveprime")

## make a plot to look at cryptic sites falling into different thresholds
dpsi_0_threshold <- as.data.frame(table(three.data.cryptic$dPSI_threshold_0))
dpsi_5_threshold <- as.data.frame(table(three.data.cryptic$dPSI_threshold_5))
pvalue_threshold <- as.data.frame(table(three.data.cryptic$pvalue_threshold))
dpsi_0_pvalue_threshold <- as.data.frame(table(three.data.cryptic$dPSI_0_pvalue_threshold))
dpsi_5_pvalue_threshold <- as.data.frame(table(three.data.cryptic$dPSI_5_pvalue_threshold))

three.df.list <- list(dpsi_0_threshold, dpsi_5_threshold, pvalue_threshold, dpsi_0_pvalue_threshold, dpsi_5_pvalue_threshold)
names(three.df.list) <- c("dPSI_0", "dPSI_5", "pvalue_0.05", "dPSI_0_pvalue_0.05", "dPSI_5_pvalue_0.05")

three.df <- bind_rows(three.df.list, .id = "Category")


dpsi_0_threshold <- as.data.frame(table(five.data.cryptic$dPSI_threshold_0))
dpsi_5_threshold <- as.data.frame(table(five.data.cryptic$dPSI_threshold_5))
pvalue_threshold <- as.data.frame(table(five.data.cryptic$pvalue_threshold))
dpsi_0_pvalue_threshold <- as.data.frame(table(five.data.cryptic$dPSI_0_pvalue_threshold))
dpsi_5_pvalue_threshold <- as.data.frame(table(five.data.cryptic$dPSI_5_pvalue_threshold))

five.df.list <- list(dpsi_0_threshold, dpsi_5_threshold, pvalue_threshold, dpsi_0_pvalue_threshold, dpsi_5_pvalue_threshold)
names(five.df.list) <- c("dPSI_0", "dPSI_5", "pvalue_0.05", "dPSI_0_pvalue_0.05", "dPSI_5_pvalue_0.05")

five.df <- bind_rows(five.df.list, .id = "Category")




## save plots and files
library(forcats)
library(RColorBrewer)
setwd(output)

write.table(three.data, "logOR_within_cell_type_ALT_3P_Junctions_with_threshold_info.txt")
write.table(five.data, "logOR_within_cell_type_ALT_5P_Junctions_with_threshold_info.txt")
write.table(three.data.cryptic, "logOR_within_cell_type_ONLY_CRYPTIC_3P_Junctions_with_threshold_info.txt")
write.table(five.data.cryptic, "logOR_within_cell_type_ONLY_CRYPTIC_5P_Junctions_with_threshold_info.txt")

threep_plot <- three.df %>%
  mutate(Var1 = fct_reorder(Var1, Category)) %>%
  ggplot(aes(x = Category, y = Freq, fill = Var1)) +
  geom_bar(stat = "identity", position = position_dodge2(width = 0.9, preserve = "single")) +
  theme_classic() +
  scale_fill_brewer(palette = "Set1") +
  labs(fill = "Pass the threshold?") +
  ggtitle("3' Cryptic Junctions") +
  theme(plot.title = element_text(hjust = 0.5), axis.text.x = element_text(angle = 45, hjust = 1))

fivep_plot <- five.df %>%
  mutate(Var1 = fct_reorder(Var1, Category)) %>%
  ggplot(aes(x = Category, y = Freq, fill = Var1)) +
  geom_bar(stat = "identity", position = position_dodge2(width = 0.9, preserve = "single")) +
  theme_classic() +
  scale_fill_brewer(palette = "Set1") +
  labs(fill = "Pass the threshold?") +
  ggtitle("5' Cryptic Junctions") +
  theme(plot.title = element_text(hjust = 0.5), axis.text.x = element_text(angle = 45, hjust = 1))

library(patchwork)
pdf("splicing_threshold_distribution.pdf", height = 5, width = 10)
threep_plot + fivep_plot
dev.off()
