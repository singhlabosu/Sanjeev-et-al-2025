library(tidyverse)
library(ggsignif)
###Making Supplementary figure
setwd("~/OneDrive - The Ohio State University/Singhlab/Bioinfo/PYMpaperFiguresUpdate/Figure5/")

UTR3quartiles<-read_tsv("ncEJCquartiles/HEKgene_groups.txt")
table(UTR3quartiles$annotation_KD_method)

SMGnorm<-UTR3quartiles%>%filter(annotation_KD_method=="50_all_SMG_normcounts")
colors <- c("#CCCCCC", "#B07C7C", "#A05252", "#801818")

##Getting last exon/UTR lengths of quartiles
#MANE isoforms
MANEisoforms<-read_tsv("../MANEexons.tsv")
MANEisoforms<-MANEisoforms%>%filter(ensembl_gene_id %in% SMGnorm$gene)

#Getting transcript info 
library(biomaRt)
listMarts(host='https://apr2020.archive.ensembl.org')
ensembl100=useMart(host='https://apr2020.archive.ensembl.org', biomart='ENSEMBL_MART_ENSEMBL', dataset='hsapiens_gene_ensembl')
attributes <- listAttributes(ensembl100)
filters <- listFilters(ensembl100)


UTRlengths<-getBM(attributes = c("ensembl_gene_id",
                                 "external_gene_name",
                                 "ensembl_transcript_id",
                                 "ensembl_exon_id",
                                 "chromosome_name",
                                 "exon_chrom_start",
                                 "exon_chrom_end",
                                 "rank",
                                 "transcript_length",
                                 "cds_length",
                                 "cds_start",
                                 "cds_end",
                                 "3_utr_start",
                                 "3_utr_end",
                                 "strand"),
                  filters = "ensembl_exon_id",
                  values = MANEisoforms$ensembl_exon_id,
                  mart = ensembl100,
                  useCache = FALSE)

UTRlengths<-UTRlengths%>%filter(ensembl_transcript_id %in% MANEisoforms$ensembl_transcript_id)%>%
  mutate(lastexonlength=abs(exon_chrom_start - exon_chrom_end),
         stopwithin = (cds_length - cds_end) == 0)
table(UTRlengths$stopwithin)

#Few entries are last exons without any CDS
list<-which(rowSums(is.na(UTRlengths)) > 0)
UTRlengths[list,]
UTRlengths[-list,]

##Calculating UTRlength
UTRlengths<-UTRlengths%>%mutate(UTR3lengthfromcoord = `3_utr_end` - `3_utr_start`,
                                UTR3lengthfromcds = lastexonlength - cds_end + cds_start,
                                UTR3match = (UTR3lengthfromcoord - UTR3lengthfromcds + 1) == 0)
table(UTRlengths$UTR3match)
#303/5872 genes have 3UTR lengths that dont match 

#plotting UTR lengths of quartiles
names(UTRlengths)[1]<-"gene"
SMGnorm<-inner_join(SMGnorm, UTRlengths[,c(1:3,19)])

#as boxplot
medians <- SMGnorm %>%
  group_by(quartile) %>%
  dplyr::summarize(median_value = median(UTR3lengthfromcds, na.rm = TRUE),
                   count = n())



s1<-SMGnorm %>%
  ggplot(aes(x = quartile,  
             y = UTR3lengthfromcds)) +
  geom_boxplot(aes(fill = quartile), varwidth = TRUE, outlier.shape = 1,outlier.size = 0.5) +
  theme_linedraw() +
  geom_signif(comparisons = list(c("Q1", "Q2")), map_signif_level = FALSE,
              y_position = 19000, tip_length = 0.005, size = 0.3,textsize=1.5) +
  geom_signif(comparisons = list(c("Q1", "Q3")), map_signif_level = FALSE,
              y_position = 21000, tip_length = 0.005, size = 0.3,textsize=1.5) +
  geom_signif(comparisons = list(c("Q1", "Q4")), map_signif_level = FALSE,
              y_position = 23000, tip_length = 0.005, size = 0.3,textsize=1.5) +
  geom_signif(comparisons = list(c("Q2", "Q3")), map_signif_level = FALSE,
              y_position = 25000, tip_length = 0.005, size = 0.3,textsize=1.5) +
  geom_signif(comparisons = list(c("Q2", "Q4")), map_signif_level = FALSE,
              y_position = 27000, tip_length = 0.005, size = 0.3,textsize=1.5) +
  geom_signif(comparisons = list(c("Q3", "Q4")), map_signif_level = FALSE,
              y_position = 29000, tip_length = 0.005, size = 0.3,textsize=1.5) +
  #coord_cartesian(ylim = c(-5, 6.5)) +
  scale_fill_manual(values = colors) +
  labs( x = "Quartiles", y= "3'UTR length") +
  theme(panel.grid.major = element_line(linetype = "blank"), 
        panel.grid.minor = element_line(linetype = "blank"),
        legend.position = "none",aspect.ratio = 2) +
  # Add median values below each boxplot
  geom_text(data = medians, aes(x = quartile, y = -400, 
                                label = sprintf("%.2f", median_value), 
                                color = quartile), 
            inherit.aes = FALSE, size = 1.5, vjust = 1) +
  scale_color_manual(values = colors)

s1
library(scales)
colors <- c("#CCCCCC", "#B07C7C", "#A05252", "#801818")
library(rstatix)
pair_df <- SMGnorm %>%
  wilcox_test(UTR3lengthfromcds ~ quartile, p.adjust.method = "bonferroni") %>%
  mutate(
    y_position = seq(4.5, 4.5 + 0.15*5, by = 0.15),
    xmin       = group1,
    xmax       = group2,
    annotation = sprintf("%.2e", p.adj) )         # format as 1.23e‑04
  

s2<-SMGnorm %>%
  ggplot(aes(x = quartile,  
             y = UTR3lengthfromcds)) +
  geom_boxplot(aes(fill = quartile), varwidth = TRUE, outlier.shape = 1,outlier.size = 0.5) +
  theme_linedraw() +
  geom_signif(
    data = pair_df,
    mapping = aes(xmin = xmin, xmax = xmax, y_position = y_position, annotations = annotation),
    manual = TRUE,
    inherit.aes = FALSE,
    tip_length = 0.005,
    size = 0.3,
    textsize = 1.5) +
  scale_y_log10(labels = label_number())+
  scale_fill_manual(values = colors) +
  labs( x = "Quartiles", y= "3'UTR length") +
  theme(panel.grid.major = element_line(linetype = "blank"), 
        panel.grid.minor = element_line(linetype = "blank"),
        legend.position = "none",aspect.ratio = 2) +
  # Add median values below each boxplot
  geom_text(data = medians, aes(x = quartile, y = 40, 
                                label = sprintf("%.2f", median_value), 
                                color = quartile), 
            inherit.aes = FALSE, size = 1.5, vjust = 1) +
  geom_text(data = medians, aes(x = quartile, y = 30, 
                                label = sprintf("%s", count)), 
            inherit.aes = FALSE, size = 1.5, vjust = 1, color="black") +
  scale_color_manual(values = colors)
s2

#Comparing with expression levels in HEKcells
SMG67kdvsControlDESeq<-read_tsv("~/OneDrive - The Ohio State University/Singhlab/Bioinfo/CASC3_RNPS1_SMG67_comparison/SMG67kdvsControlDESeq.tsv")
names(SMG67kdvsControlDESeq)[1]<-"ensembl_transcript_id"
SMGnorm<-inner_join(SMGnorm,SMG67kdvsControlDESeq)


medians <- SMGnorm %>%
  group_by(quartile) %>%
  dplyr::summarize(median_value = median(baseMean, na.rm = TRUE),
                   count = n())
pair_df <- SMGnorm %>%
  wilcox_test(baseMean ~ quartile, p.adjust.method = "bonferroni") %>%
  mutate(
    y_position = seq(5.2, 5.2 + 0.25*5, by = 0.25),
    xmin       = group1,
    xmax       = group2,
    annotation = sprintf("%.2e", p.adj) )  


  
s3<-SMGnorm %>%
  ggplot(aes(x = quartile,  
             y = baseMean)) +
  geom_boxplot(aes(fill = quartile), varwidth = TRUE, outlier.shape = 1,outlier.size = 0.5) +
  theme_linedraw()+
  geom_signif(
    data = pair_df,
    mapping = aes(xmin = xmin, xmax = xmax, y_position = y_position, annotations = annotation),
    manual = TRUE,
    inherit.aes = FALSE,
    tip_length = 0.005,
    size = 0.3,
    textsize = 1.5) +
  scale_y_log10(labels = label_number())+
  scale_fill_manual(values = colors) +
  labs( x = "Quartiles", y= "HEK293  baseMean") +
  theme(panel.grid.major = element_line(linetype = "blank"), 
        panel.grid.minor = element_line(linetype = "blank"),
        legend.position = "none",aspect.ratio = 2) +
  # Add median values below each boxplot
  geom_text(data = medians, aes(x = quartile, y = 4, 
                                label = sprintf("%.2f", median_value), 
                                color = quartile), 
            inherit.aes = FALSE, size = 1.5, vjust = 1) +
  geom_text(data = medians, aes(x = quartile, y = 3, 
                                label = sprintf("%s", count)), 
            inherit.aes = FALSE, size = 1.5, vjust = 1, color="black")+
  scale_color_manual(values = colors)
s3
s2

##plotting for HCTcells

UTR3quartiles2<-read_tsv("ncEJCquartiles/HCTgene_groups.txt")
table(UTR3quartiles2$annotation_KD_method)

UPF3norm<-UTR3quartiles2%>%filter(annotation_KD_method=="50_all_UPF3_normcounts")
colors <- colorRampPalette(c("#CCCCCC", "darkblue"))(4)
colors

# Define quartile counts
quartile_counts <- UPF3norm %>%
  group_by(quartile) %>%
  dplyr::summarise(n = n()) %>%
  mutate(label = paste0(quartile, " (n=", n, ")"))

# Map counts to color labels
quartile_labels <- setNames(quartile_counts$label, quartile_counts$quartile)

# Wilcoxon test p-values
p_values <- pairwise.wilcox.test(
  UPF3norm$log2foldchange, 
  UPF3norm$quartile, 
  p.adjust.method = "bonferroni"
)$p.value

p_Q1_Q2 <- p_values[[1]]
p_Q2_Q3 <-p_values[2,2]
p_Q3_Q4 <-p_values[3,3]
  
library(ggtext)  # make sure it's installed

# Rich text labels with color-coded Q1–Q4
labels <- data.frame(
  x = -3.2,
  y = c(1,0.95, 0.90),
  label = c(
    paste0("<span style='color:", colors[1], "'>Q1</span> vs <span style='color:", colors[2], "'>Q2</span>: p = ", signif(p_Q1_Q2, 3)),
    paste0("<span style='color:", colors[2], "'>Q2</span> vs <span style='color:", colors[3], "'>Q3</span>: p = ", signif(p_Q2_Q3, 3)),
    paste0("<span style='color:", colors[3], "'>Q3</span> vs <span style='color:", colors[4], "'>Q4</span>: p = ", signif(p_Q3_Q4, 3))
  )
)

# Plot
s4<-UPF3norm %>%
  ggplot(aes(log2foldchange, color = quartile)) +
  geom_vline(aes(xintercept = 0), linewidth = 0.02) +
  geom_hline(aes(yintercept = 0.5), linewidth = 0.02) +
  stat_ecdf(geom = "line", linewidth = 0.8) +
  coord_cartesian(xlim = c(-3, 3)) +
  scale_color_manual(values = colors, labels = quartile_labels, name = "Quartile") +
  geom_richtext(
    data = labels,
    aes(x = x, y = y, label = label),
    inherit.aes = FALSE,
    fill = NA, label.color = NA,  # Remove box and border
    size = 1.5, hjust = 0
  ) +
  theme_linedraw() +
  theme(
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    aspect.ratio = 1,
    legend.position = c(0.95, 0.05),
    legend.justification = c(1, 0),
    legend.title = element_text(size = 6), 
    legend.text = element_text(size = 4),
    legend.key = element_blank(),
    legend.box.background = element_blank(),
    legend.background = element_blank(),
    legend.key.size = unit(2, 'mm')
  ) +
  labs(
    y = "Cumulative frequency", 
    x = "log2FC SMG7KO+SMG6kd vs WT"
  )
s4


library(ggsignif)

# s4<-UPF3norm %>%
#   ggplot(aes(x = quartile,  
#              y = log2foldchange)) +
#   geom_boxplot(aes(fill = quartile), varwidth = TRUE, outlier.shape = 1,outlier.size = 0.5) +
#   theme_linedraw() +
#   geom_signif(comparisons = list(c("Q1", "Q2")), map_signif_level = FALSE,
#               y_position = 3.0, tip_length = 0.005, size = 0.3,textsize=1.5) +
#   geom_signif(comparisons = list(c("Q1", "Q3")), map_signif_level = FALSE,
#               y_position = 3.5, tip_length = 0.005, size = 0.3,textsize=1.5) +
#   geom_signif(comparisons = list(c("Q1", "Q4")), map_signif_level = FALSE,
#               y_position = 4.0, tip_length = 0.005, size = 0.3,textsize=1.5) +
#   geom_signif(comparisons = list(c("Q2", "Q3")), map_signif_level = FALSE,
#               y_position = 4.5, tip_length = 0.005, size = 0.3,textsize=1.5) +
#   geom_signif(comparisons = list(c("Q2", "Q4")), map_signif_level = FALSE,
#               y_position = 5.0, tip_length = 0.005, size = 0.3,textsize=1.5) +
#   geom_signif(comparisons = list(c("Q3", "Q4")), map_signif_level = FALSE,
#               y_position = 5.5, tip_length = 0.005, size = 0.3,textsize=1.5) +
#   coord_cartesian(ylim = c(-5, 6.5)) +
#   scale_fill_manual(values = colors) +
#   labs( x = "Quartiles", y= "log2FC UPF3dKO vs WT") +
#   theme(panel.grid.major = element_line(linetype = "blank"), 
#         panel.grid.minor = element_line(linetype = "blank"),
#         legend.position = "none",aspect.ratio = 2) +
#   # Add median values below each boxplot
#   geom_text(data = medians, aes(x = quartile, y = -3.87, 
#                                 label = sprintf("%.2f", median_value), 
#                                 color = quartile), 
#             inherit.aes = FALSE, size = 1.5, vjust = 1) +
#   scale_color_manual(values = colors)
# 
# s4

#plotting UTR lengths of quartiles
UPF3norm<-inner_join(UPF3norm, UTRlengths[,c(1,2,19)])

#as boxplot
medians <- UPF3norm %>%
  group_by(quartile) %>%
  dplyr::summarize(median_value = median(UTR3lengthfromcds, na.rm = TRUE),
                   count = n())



s5<-UPF3norm %>%
  ggplot(aes(x = quartile,  
             y = UTR3lengthfromcds)) +
  geom_boxplot(aes(fill = quartile), varwidth = TRUE, outlier.shape = 1,outlier.size = 0.5) +
  theme_linedraw() +
  geom_signif(comparisons = list(c("Q1", "Q2")), map_signif_level = FALSE,
              y_position = 19000, tip_length = 0.005, size = 0.3,textsize=1.5) +
  geom_signif(comparisons = list(c("Q1", "Q3")), map_signif_level = FALSE,
              y_position = 21000, tip_length = 0.005, size = 0.3,textsize=1.5) +
  geom_signif(comparisons = list(c("Q1", "Q4")), map_signif_level = FALSE,
              y_position = 23000, tip_length = 0.005, size = 0.3,textsize=1.5) +
  geom_signif(comparisons = list(c("Q2", "Q3")), map_signif_level = FALSE,
              y_position = 25000, tip_length = 0.005, size = 0.3,textsize=1.5) +
  geom_signif(comparisons = list(c("Q2", "Q4")), map_signif_level = FALSE,
              y_position = 27000, tip_length = 0.005, size = 0.3,textsize=1.5) +
  geom_signif(comparisons = list(c("Q3", "Q4")), map_signif_level = FALSE,
              y_position = 29000, tip_length = 0.005, size = 0.3,textsize=1.5) +
  #coord_cartesian(ylim = c(-5, 6.5)) +
  scale_fill_manual(values = colors) +
  labs( x = "Quartiles", y= "3'UTR length") +
  theme(panel.grid.major = element_line(linetype = "blank"), 
        panel.grid.minor = element_line(linetype = "blank"),
        legend.position = "none",aspect.ratio = 2) +
  # Add median values below each boxplot
  geom_text(data = medians, aes(x = quartile, y = -400, 
                                label = sprintf("%.2f", median_value), 
                                color = quartile), 
            inherit.aes = FALSE, size = 1.5, vjust = 1) +
  geom_text(data = medians, aes(x = quartile, y = -990, 
                                label = sprintf("%s", count)), 
            inherit.aes = FALSE, size = 1.5, vjust = 1, color="black")+
  scale_color_manual(values = colors)
s5

pair_df <- UPF3norm %>%
  wilcox_test(UTR3lengthfromcds ~ quartile, p.adjust.method = "bonferroni") %>%
  mutate(
    y_position = seq(4.5, 4.5 + 0.15*5, by = 0.15),
    xmin       = group1,
    xmax       = group2,
    annotation = sprintf("%.2e", p.adj))          # format as 1.23e‑04


s6<-UPF3norm %>%
  ggplot(aes(x = quartile,  
             y = UTR3lengthfromcds)) +
  geom_boxplot(aes(fill = quartile), varwidth = TRUE, outlier.shape = 1,outlier.size = 0.5) +
  theme_linedraw() +
  geom_signif(
    data = pair_df,
    mapping = aes(xmin = xmin, xmax = xmax, y_position = y_position, annotations = annotation),
    manual = TRUE,
    inherit.aes = FALSE,
    tip_length = 0.005,
    size = 0.3,
    textsize = 1.5)+
  scale_y_log10(labels = label_number())+
  scale_fill_manual(values = colors) +
  labs( x = "Quartiles", y= "3'UTR length") +
  theme(panel.grid.major = element_line(linetype = "blank"), 
        panel.grid.minor = element_line(linetype = "blank"),
        legend.position = "none",aspect.ratio = 2) +
  # Add median values below each boxplot
  geom_text(data = medians, aes(x = quartile, y = 85, 
                                label = sprintf("%.2f", median_value), 
                                color = quartile), 
            inherit.aes = FALSE, size = 1.5, vjust = 1) +
  geom_text(data = medians, aes(x = quartile, y = 70, 
                                label = sprintf("%s", count)), 
            inherit.aes = FALSE, size = 1.5, vjust = 1, color="black")+
  scale_color_manual(values = colors)
s6

#Comparing with expression levels in HCTcells
HctUPF1DESEq<-read_csv("~/Downloads/GSE179843_RF_RUVDESeq_WTsiUPF1vWTsiNC_kallisto_TPM_v100.csv")
names(HctUPF1DESEq)[2]<-"gene"
names(HctUPF1DESEq)[3]<-"baseMeanHCT"
HctUPF1DESEq<-HctUPF1DESEq%>%mutate(gene= str_split(gene, "\\.", simplify = T)[,1])
HctUPF1DESEq<-HctUPF1DESEq%>%mutate(tx_id= str_split(tx_id, "\\.", simplify = T)[,1])

HctUPF1DESEq<-HctUPF1DESEq%>%filter(tx_id %in% MANEisoforms$ensembl_transcript_id)


UPF3norm<-left_join(UPF3norm, HctUPF1DESEq[,c(2,3)])
medians <- UPF3norm %>%
  group_by(quartile) %>%
  dplyr::summarize(median_value = median(baseMeanHCT, na.rm = TRUE),
                   count = n())

pair_df <- UPF3norm %>%
  wilcox_test(baseMeanHCT ~ quartile, p.adjust.method = "bonferroni") %>%
  mutate(
    y_position = seq(5.2, 5.2 + 0.25*5, by = 0.25),
    xmin       = group1,
    xmax       = group2,
    annotation = sprintf("%.2e", p.adj))          # format as 1.23e‑04


s7<-UPF3norm %>%
  ggplot(aes(x = quartile,  
             y = baseMeanHCT)) +
  geom_boxplot(aes(fill = quartile), varwidth = TRUE, outlier.shape = 1,outlier.size = 0.5) +
  theme_linedraw()+
  geom_signif(
    data = pair_df,
    mapping = aes(xmin = xmin, xmax = xmax, y_position = y_position, annotations = annotation),
    manual = TRUE,
    inherit.aes = FALSE,
    tip_length = 0.005,
    size = 0.3,
    textsize = 1.5)+
  scale_y_log10(labels = label_number())+
  scale_fill_manual(values = colors) +
  labs( x = "Quartiles", y= "HCT116  baseMean") +
  theme(panel.grid.major = element_line(linetype = "blank"), 
        panel.grid.minor = element_line(linetype = "blank"),
        legend.position = "none",aspect.ratio = 2) +
  # Add median values below each boxplot
  geom_text(data = medians, aes(x = quartile, y = 4, 
                                label = sprintf("%.2f", median_value), 
                                color = quartile), 
            inherit.aes = FALSE, size = 1.5, vjust = 1) +
  geom_text(data = medians, aes(x = quartile, y = 3, 
                                label = sprintf("%s", count)), 
            inherit.aes = FALSE, size = 1.5, vjust = 1, color="black")+
  
  scale_color_manual(values = colors)

s7

colors <- c("#CCCCCC", "#B07C7C", "#A05252", "#801818")
E117RvsWT<-read_csv("~/OneDrive - The Ohio State University/Singhlab/Bioinfo/PYMpaperFiguresUpdate/DESeq2 for PYM samples/all_NT_res.csv")

E117RvsWT<- E117RvsWT%>%dplyr::rename(gene = ...1)%>%na.omit()
E117RvsWT<-inner_join(E117RvsWT, SMGnorm[,2:4])



#as boxplot
medians <- E117RvsWT %>%
  group_by(quartile) %>%
  dplyr::summarize(median_value = median(log2FoldChange, na.rm = TRUE),
                   count = n())

pair_df <- E117RvsWT %>%
  wilcox_test(log2FoldChange ~ quartile, p.adjust.method = "bonferroni") %>%
  mutate(
    y_position = seq(3.1, 3.1 + 0.4*5, by = 0.4),
    xmin       = group1,
    xmax       = group2,
    annotation = sprintf("%.2e", p.adj))          # format as 1.23e‑04

s8<-E117RvsWT %>%
  ggplot(aes(x = quartile,  
             y = log2FoldChange)) +
  geom_boxplot(aes(fill = quartile), varwidth = TRUE, outlier.shape = 1,outlier.size = 0.5) +
  theme_linedraw() +
  geom_signif(
    data = pair_df,
    mapping = aes(xmin = xmin, xmax = xmax, y_position = y_position, annotations = annotation),
    manual = TRUE,
    inherit.aes = FALSE,
    tip_length = 0.005,
    size = 0.3,
    textsize = 1.5) +
  coord_cartesian(ylim = c(-4, 5.3)) +
  scale_fill_manual(values = colors) +
  labs( x = "Quartiles", y= "log2FC E117R vs WT") +
  theme(panel.grid.major = element_line(linetype = "blank"), 
        panel.grid.minor = element_line(linetype = "blank"),
        legend.position = "none",aspect.ratio = 2) +
  # Add median values below each boxplot
  geom_text(data = medians, aes(x = quartile, y = -3.7, 
                                label = sprintf("%.2f", median_value), 
                                color = quartile), 
            inherit.aes = FALSE, size = 2, vjust = 1) +
  geom_text(data = medians, aes(x = quartile, y = -4, 
                                label = sprintf("%s", count)), 
            inherit.aes = FALSE, size = 1.5, vjust = 1, color="black")+
  scale_color_manual(values = colors)

s8

# 
# 
# ##plotting for HeLacells
# 
# UTR3quartiles3<-read_tsv("ncEJCquartiles/HeLagene_groups.txt")
# table(UTR3quartiles3$annotation_KD_method)
# 
# UPF1norm<-UTR3quartiles3%>%filter(annotation_KD_method=="50_all_UPF1_normcounts")
# colors <- colorRampPalette(c("#CCCCCC", "darkgreen"))(4)
# colors
# 
# # Define quartile counts
# quartile_counts <- UPF1norm %>%
#   group_by(quartile) %>%
#   dplyr::summarise(n = n()) %>%
#   mutate(label = paste0(quartile, " (n=", n, ")"))
# 
# # Map counts to color labels
# quartile_labels <- setNames(quartile_counts$label, quartile_counts$quartile)
# 
# # Wilcoxon test p-values
# p_Q1_Q2 <- wilcox.test(log2foldchange ~ quartile, 
#                        data = UPF1norm %>% filter(quartile %in% c("Q1", "Q2")))$p.value
# p_Q2_Q3 <- wilcox.test(log2foldchange ~ quartile, 
#                        data = UPF1norm %>% filter(quartile %in% c("Q2", "Q3")))$p.value
# p_Q3_Q4 <- wilcox.test(log2foldchange ~ quartile, 
#                        data = UPF1norm %>% filter(quartile %in% c("Q3", "Q4")))$p.value
# 
# library(ggtext)  # make sure it's installed
# 
# # Rich text labels with color-coded Q1–Q4
# labels <- data.frame(
#   x = -3.2,
#   y = c(1,0.95, 0.90),
#   label = c(
#     paste0("<span style='color:", colors[1], "'>Q1</span> vs <span style='color:", colors[2], "'>Q2</span>: p = ", signif(p_Q1_Q2, 3)),
#     paste0("<span style='color:", colors[2], "'>Q2</span> vs <span style='color:", colors[3], "'>Q3</span>: p = ", signif(p_Q2_Q3, 3)),
#     paste0("<span style='color:", colors[3], "'>Q3</span> vs <span style='color:", colors[4], "'>Q4</span>: p = ", signif(p_Q3_Q4, 3))
#   )
# )
# 
# # Plot
# s9<-UPF1norm %>%
#   ggplot(aes(log2foldchange, color = quartile)) +
#   geom_vline(aes(xintercept = 0), linewidth = 0.02) +
#   geom_hline(aes(yintercept = 0.5), linewidth = 0.02) +
#   stat_ecdf(geom = "line", linewidth = 0.8) +
#   coord_cartesian(xlim = c(-1, 1)) +
#   scale_color_manual(values = colors, labels = quartile_labels, name = "Quartile") +
#   geom_richtext(
#     data = labels,
#     aes(x = x, y = y, label = label),
#     inherit.aes = FALSE,
#     fill = NA, label.color = NA,  # Remove box and border
#     size = 2, hjust = 0
#   ) +
#   theme_linedraw() +
#   theme(
#     panel.grid.major = element_blank(), 
#     panel.grid.minor = element_blank(),
#     aspect.ratio = 1,
#     legend.position = c(0.95, 0.05),
#     legend.justification = c(1, 0),
#     legend.title = element_text(size = 6), 
#     legend.text = element_text(size = 4),
#     legend.key = element_blank(),
#     legend.box.background = element_blank(),
#     legend.background = element_blank(),
#     legend.key.size = unit(2, 'mm')
#   ) +
#   labs(
#     y = "Cumulative frequency", 
#     x = "log2FC UPF1kd vs WT"
#   )
# s9
# 
# medians <- UPF1norm %>%
#   group_by(quartile) %>%
#   dplyr::summarize(median_value = median(log2foldchange, na.rm = TRUE))
# 
# s10<-UPF1norm %>%
#   ggplot(aes(x = quartile,
#              y = log2foldchange)) +
#   geom_boxplot(aes(fill = quartile), varwidth = TRUE, outlier.shape = 1,outlier.size = 0.5) +
#   theme_linedraw() +
#   geom_signif(comparisons = list(c("Q1", "Q2")), map_signif_level = FALSE,
#               y_position = 3.0, tip_length = 0.005, size = 0.3,textsize=1.5) +
#   geom_signif(comparisons = list(c("Q1", "Q3")), map_signif_level = FALSE,
#               y_position = 3.5, tip_length = 0.005, size = 0.3,textsize=1.5) +
#   geom_signif(comparisons = list(c("Q1", "Q4")), map_signif_level = FALSE,
#               y_position = 4.0, tip_length = 0.005, size = 0.3,textsize=1.5) +
#   geom_signif(comparisons = list(c("Q2", "Q3")), map_signif_level = FALSE,
#               y_position = 4.5, tip_length = 0.005, size = 0.3,textsize=1.5) +
#   geom_signif(comparisons = list(c("Q2", "Q4")), map_signif_level = FALSE,
#               y_position = 5.0, tip_length = 0.005, size = 0.3,textsize=1.5) +
#   geom_signif(comparisons = list(c("Q3", "Q4")), map_signif_level = FALSE,
#               y_position = 5.5, tip_length = 0.005, size = 0.3,textsize=1.5) +
#   coord_cartesian(ylim = c(-5, 7)) +
#   scale_fill_manual(values = colors) +
#   labs( x = "Quartiles", y= "log2FC UPF1kd vs WT") +
#   theme(panel.grid.major = element_line(linetype = "blank"),
#         panel.grid.minor = element_line(linetype = "blank"),
#         legend.position = "none",aspect.ratio = 2) +
#   # Add median values below each boxplot
#   geom_text(data = medians, aes(x = quartile, y = -3.87,
#                                 label = sprintf("%.2f", median_value),
#                                 color = quartile),
#             inherit.aes = FALSE, size = 1.5, vjust = 1) +
#   scale_color_manual(values = colors)
# 
# s10

library(patchwork)
layout <- c(
  area(1,1,2,2),area(1, 3, 2,4),area(1,5,2,7),
  area(3,1,4,2),area(3,3,4,4),
  area(3,5,4,6),area(5,4,8,6))
plot(layout)

ggsave("S5patch1.pdf",
       s2+s3+s4+s6+s7+s8+
         plot_layout(design = layout)+
         plot_annotation(title = 'FIGURE S5', tag_levels = "a")  &
         theme(plot.tag = element_text(face = 'bold')),
       width = 8.5, 
       height = 11, 
       units = "in")


