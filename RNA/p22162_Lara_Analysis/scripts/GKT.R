library(tidyverse)
library(readr)
library(openxlsx)
library(ggforce)
library(ggpubr)


setwd("/yerkes-cifs/runs/Analysis/2022_Analyses/p22162_Lara/p22162_Lara_Analysis/")

res_p22162_Lara <- read.xlsx("/yerkes-cifs/runs/Analysis/2022_Analyses/p22162_Lara/p22162_Lara_Analysis/p22162_Lara_Analysis_DESeq2_results.xlsx")
wb_p22162_Lara <- loadWorkbook("/yerkes-cifs/runs/Analysis/2022_Analyses/p22162_Lara/p22162_Lara_Analysis/p22162_Lara_Analysis_DESeq2_results.xlsx")
worksheetOrder(wb_p22162_Lara)
getSheetNames("/yerkes-cifs/runs/Analysis/2022_Analyses/p22162_Lara/p22162_Lara_Analysis/p22162_Lara_Analysis_DESeq2_results.xlsx")
sample_table_p22162_Lara <- read_tsv("/yerkes-cifs/runs/Analysis/2022_Analyses/p22162_Lara/p22162_Lara_Analysis/p22162_Lara_exp_design.txt")
rlog_forGSEA_p22172_Lara <- read_tsv("/yerkes-cifs/runs/Analysis/2022_Analyses/p22162_Lara/p22162_Lara_Analysis/rlog_forGSEA_p22162_Lara_Analysis.txt")
rlog_selected <- filter(rlog_forGSEA_p22172_Lara, Name %in% c("NECTIN4","CD274","TACSTD2","EFNB2")) %>%
  mutate(gene = Name,.before = Description,.keep = "unused") %>%
  select(!Description) %>%
  rename_with(~ str_replace(.x, "Emory.","Emory-"))
rlog_selected <- rlog_selected[-2,] # remove duplicate NECTIN4 row with pactically no expression
plot_table <- pivot_longer(rlog_selected,
                           cols = starts_with("Emory-"),
                           names_to = "SampleID",
                           values_to = "log2_Exp") %>%
  left_join(select(sample_table_p22162_Lara,SampleID,Group)) %>%
  filter(!str_detect(Group,"Mixed"))

facet(ggboxplot(plot_table,x = "Group",y="log2_Exp",add = c("jitter"), outlier.shape = NA), facet.by = "gene") +
  stat_compare_means(comparisons = combn(unique(plot_table$Group),2, simplify = FALSE)[1:4]) +
  theme_bw()


