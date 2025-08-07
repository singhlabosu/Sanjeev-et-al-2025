###Figure5 plots

library(tidyverse)
library(ggrepel)
library(patchwork)
library(tximport)
library(DESeq2)
library(ggsignif)



setwd("~/OneDrive - The Ohio State University/Singhlab/Bioinfo/PYMpaperFiguresUpdate/Figure5/")

#(NMD figure)
##Plotting PYMkd-NMD and other plots#####
#making a dataframe containing some metadata
# sampleinfo<-as.data.frame(list.files(path = "Kallisto/"))
# colnames(sampleinfo)<-"Files"
# sampleinfo$Sample<-rep(c("Control", "PYM KD","PYM OE","PYMΔΝ"), each=3)
# 
# 
# #importing kallisto output data into R
# files1 <- file.path("~/OneDrive - The Ohio State University/Singhlab/Bioinfo/PYMpaperFigures/Figure4/Kallisto", sampleinfo$Files, "abundance.h5")
# names(files1) <- sampleinfo$Files
# files1
# 
# txi.kallisto.tsv <- tximport(files1, type = "kallisto", txOut=TRUE, countsFromAbundance = "lengthScaledTPM", ignoreAfterBar = TRUE)
# #head(txi.kallisto.tsv$counts)
# counts<-as.data.frame(txi.kallisto.tsv$counts)
# 
# filter <- apply(counts, 1, function(x) length(x[x>0.1])>=12)
# counts<-counts[filter,]
# 
# #Differential expression analysis####
# #sampleinfo as dataframe
# sampleinfo$Sample <- as.factor(sampleinfo$Sample)
# 
# 
# ddscounts.table <- DESeqDataSetFromTximport(txi.kallisto.tsv, sampleinfo, formula(~ Sample))
# ddscounts.table = ddscounts.table[rownames(ddscounts.table) %in% rownames(counts)]
# ddscounts <- DESeq(ddscounts.table)

# #PCAplots
# vsd <- vst(ddscounts, blind = FALSE)
# plotPCA(vsd, intgroup = "Sample", ntop=500)
# 
# #PYMkd DEseq#####
# #Results
# #####1st pairwise comparison
# res_counts1 <- results(ddscounts, contrast = c("Sample","PYM KD","Control"))
# PYMkdvsControlDESeq<-as.data.frame(res_counts1)
# PYMkdvsControlDESeq<-na.omit(PYMkdvsControlDESeq)
# plotMA(res_counts1, ylim=c(-7,7))
# PYMkdvsControlDESeq$ensembl_transcript_id<-rownames(PYMkdvsControlDESeq)
# write_tsv(PYMkdvsControlDESeq, "DESeqResults/PYMkdvsControlDESeq_tpm0.1.tsv")
PYMkdvsControlDESeq<-read_tsv("DESeqResults/PYMkdvsControlDESeq_tpm0.1.tsv")


#PYMkd_ all PTC cdf plot####
#analysis with robert's PTC+ and PTC- list
PTClist <- read.delim("~/OneDrive - The Ohio State University/Singhlab/Bioinfo/CASC3/ENST_PTC-EPI-TFG.txt")

names(PYMkdvsControlDESeq)[7] <- "ENST.ID"
PYMkdvsControlDESeq$ENST.ID <- sub("\\..*", "", PYMkdvsControlDESeq$ENST.ID)
length(intersect(PYMkdvsControlDESeq$ENST.ID,PTClist$ENST.ID))

PYMkdvsControlDESeq<-inner_join(PYMkdvsControlDESeq, PTClist, by = "ENST.ID")

library(biomaRt)
listMarts(host='https://apr2020.archive.ensembl.org')
ensembl100=useMart(host='https://apr2020.archive.ensembl.org', biomart='ENSEMBL_MART_ENSEMBL', dataset='hsapiens_gene_ensembl')
attributes <- listAttributes(ensembl100)
filters <- listFilters(ensembl100)

PTCgenes <- getBM(attributes = c("ensembl_gene_id","ensembl_transcript_id",
                                 "transcript_biotype"),
                  filters = "ensembl_transcript_id",
                  values = PYMkdvsControlDESeq$ENST.ID,
                  mart = ensembl100,
                  useCache = FALSE)
table(PTCgenes$transcript_biotype)
names(PTCgenes)[2]<-"ENST.ID"
PYMkdvsControlDESeq<-inner_join(PYMkdvsControlDESeq, PTCgenes, by = "ENST.ID")

PYMkdvsControlDESeq <- PYMkdvsControlDESeq %>%
  dplyr::group_by(ensembl_gene_id) %>%                      # Group by gene ID
  filter(any(PTC.Status == TRUE) & any(PTC.Status == FALSE)) %>%  # Keep only genes with both TRUE & FALSE PTC statuses
  ungroup()



res <- wilcox.test(log2FoldChange ~ PTC.Status, data = PYMkdvsControlDESeq,
                   exact = FALSE, alternative = "less")
res$p.value
table(PYMkdvsControlDESeq$PTC.Status)

#Plot:PYMkdvsControl_allPTC####
# calculate pairwise p-values (for example, Wilcoxon test)
p_values <- pairwise.wilcox.test(
  PYMkdvsControlDESeq$log2FoldChange, 
  PYMkdvsControlDESeq$PTC.Status, 
  p.adjust.method = "bonferroni"
)$p.value

# Convert p-values to a readable format for display
# Extract and format p-values correctly based on row and column names
pval_text <- paste("p-values(Wilcoxon test):",
                   paste(format(p_values, scientific = TRUE, digits = 2)),  # "c vs l"
                   sep = "\n")

# Calculate the number of genes in each group
group_counts <- PYMkdvsControlDESeq %>%
  group_by(PTC.Status) %>%
  dplyr::summarize(count = n(), .groups = 'drop')

# Create formatted labels with counts
labels_with_counts <- paste(group_counts$PTC.Status, 
                            " (", group_counts$count, ")", sep = "")
labels_with_counts <-str_replace(labels_with_counts,"TRUE","PTC+")
labels_with_counts <-str_replace(labels_with_counts,"FALSE","PTC-")

# plotting as before
PYMkdvsControlDESeq %>%
  ggplot(aes(log2FoldChange, color = PTC.Status)) +
  geom_vline(aes(xintercept = 0), linewidth = 0.02)+
  geom_hline(aes(yintercept = 0.5), linewidth = 0.02)+
  stat_ecdf(geom = "line", linewidth = 0.8) +
  scale_color_manual(values = c("FALSE" = "grey45", "TRUE" = "red3"), 
                     labels = labels_with_counts) +  # Change legend text
  coord_cartesian(xlim = c(-2, 2)) +
  theme_linedraw() +
  theme(
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    aspect.ratio = 1,
    legend.position = c(0.95, 0.05),         # Move legend to bottom right
    legend.justification = c(1, 0)     # Ensure correct alignment at bottom right
  ) +
  labs(title = "PYM1 kd vs Control",
       y = "Cumulative frequency", x = "log2FC PYM1kd/Control") +
  annotate("text", x = -1.9, y = 0.92, label = pval_text, hjust = 0, size = 3, color = "black")  # Add p-values on top left


#Comparing with SMG6/7kd#####

SMG67kdvsControlDESeq<-read_tsv("~/OneDrive - The Ohio State University/Singhlab/Bioinfo/CASC3_RNPS1_SMG67_comparison/SMG67kdvsControlDESeq.tsv")

SMG67kdvsControlDESeq<-inner_join(SMG67kdvsControlDESeq,PTClist, by = "ENST.ID")
SMG67kdvsControlDESeq<-na.omit(SMG67kdvsControlDESeq)
table(SMG67kdvsControlDESeq$PTC.Status)

#Filtering for upregulated PTC+ and unchanged PTC-
PTCgenes_filtered <- SMG67kdvsControlDESeq %>%
  filter((PTC.Status == "FALSE") & log2FoldChange >= -0.58 & log2FoldChange <= 0.58 |
           (PTC.Status == "TRUE" & log2FoldChange > 0.58))
table(PTCgenes_filtered$PTC.Status)

##Making sure PTC+ genes that passed filtering has isoforms in PTC-

#Getting transcript info
# library(biomaRt)
# listMarts(host='https://apr2020.archive.ensembl.org')
# ensembl100=useMart(host='https://apr2020.archive.ensembl.org', biomart='ENSEMBL_MART_ENSEMBL', dataset='hsapiens_gene_ensembl')
# attributes <- listAttributes(ensembl100)
# filters <- listFilters(ensembl100)

PTCgenes <- getBM(attributes = c("ensembl_gene_id","ensembl_transcript_id",
                                 "transcript_biotype"),
                  filters = "ensembl_transcript_id",
                  values = PTCgenes_filtered$ENST.ID,
                  mart = ensembl100,
                  useCache = FALSE)

names(PTCgenes)[2]<-"ENST.ID"
table(PTCgenes$transcript_biotype)
PTCgenes_filtered<-inner_join(PTCgenes_filtered,PTCgenes,by = "ENST.ID")

# table(PTCgenes_filtered$PTC.Status)
# table(PTCgenes_filtered$transcript_biotype)
# table(PTCgenes_filtered$ensembl_gene_id)
# 
# # Create a cross-tabulated table of counts
# summary_counts <- as.data.frame(table(PTCgenes_filtered$transcript_biotype, PTCgenes_filtered$PTC.Status))
# summary_counts
PTCgenes_filtered <- PTCgenes_filtered %>%
  dplyr::group_by(ensembl_gene_id) %>%                      # Group by gene ID
  filter(any(PTC.Status == TRUE) & any(PTC.Status == FALSE)) %>%  # Keep only genes with both TRUE & FALSE PTC statuses
  ungroup()

write_tsv(PTCgenes_filtered,"DESeqResults/PTCgenes_SMGfiltered.tsv")


#PYM kd Filtered PTC+/- CDF#####

#Removing transcriptids that dont pass SMG6/7 filtering
PTCgenes_filtered<-read_tsv("DESeqResults/PTCgenes_SMGfiltered.tsv")
glimpse(PTCgenes_filtered)
glimpse(PYMkdvsControlDESeq)
FilteredPYMkdvsControlDESeq <- PYMkdvsControlDESeq%>%filter(ENST.ID %in% PTCgenes_filtered$ENST.ID)
gene_counts <- as.data.frame(table(FilteredPYMkdvsControlDESeq$ensembl_gene_id))
FilteredPYMkdvsControlDESeq <- FilteredPYMkdvsControlDESeq %>%
  dplyr::group_by(ensembl_gene_id) %>%                      # Group by gene ID
  filter(any(PTC.Status == TRUE) & any(PTC.Status == FALSE)) %>%  # Keep only genes with both TRUE & FALSE PTC statuses
  ungroup()
write_tsv(FilteredPYMkdvsControlDESeq,"DESeqResults/FilteredPYMkdvsControlDESeq.tsv")
FilteredPYMkdvsControlDESeq<-read_tsv("DESeqResults/FilteredPYMkdvsControlDESeq.tsv")

# calculate pairwise p-values (for example, Wilcoxon test)
p_values <- pairwise.wilcox.test(
  FilteredPYMkdvsControlDESeq$log2FoldChange, 
  FilteredPYMkdvsControlDESeq$PTC.Status, 
  p.adjust.method = "bonferroni"
)$p.value

# Convert p-values to a readable format for display
# Extract and format p-values correctly based on row and column names
pval_text <- paste("p-values(Wilcoxon test):",
                   paste(format(p_values, digits = 2)),  # "c vs l"
                   sep = "\n")

# Calculate the number of genes in each group
group_counts <- FilteredPYMkdvsControlDESeq %>%
  group_by(PTC.Status) %>%
  dplyr::summarize(count = n(), .groups = 'drop')

# Create formatted labels with counts
labels_with_counts <- paste(group_counts$PTC.Status, 
                            " (", group_counts$count, ")", sep = "")
labels_with_counts <-str_replace(labels_with_counts,"TRUE","PTC+")
labels_with_counts <-str_replace(labels_with_counts,"FALSE","PTC-")

# plotting as before
a<-FilteredPYMkdvsControlDESeq %>%
  ggplot(aes(log2FoldChange, color = PTC.Status)) +
  geom_vline(aes(xintercept = 0), linewidth = 0.02)+
  geom_hline(aes(yintercept = 0.5), linewidth = 0.02)+
  stat_ecdf(geom = "line", linewidth = 0.8) +
  scale_color_manual( "PTC Status",values = c("FALSE" = "grey45", "TRUE" = "red3"), 
                     labels = labels_with_counts) +  # Change legend text
  coord_cartesian(xlim = c(-2, 2)) +
  theme_linedraw() +
  theme(
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    aspect.ratio = 1,
    legend.position = c(0.98, 0.02),         # Move legend to bottom right
    legend.justification = c(1, 0),    # Ensure correct alignment at bottom right
    legend.title=element_text(size=7), 
    legend.text=element_text(size=5),
    legend.key.size =unit(2, 'mm')     
  ) +
  labs(y = "Cumulative frequency", x = "log2FC PYM1kd/Control") +
  annotate("text", x = -1.9, y = 0.92, label = pval_text, hjust = 0, size = 2, color = "black")  # Add p-values on top left
a+labs(title = "PYM1kd vs Control",  caption = "SMG6/7 validated")


res <- wilcox.test(log2FoldChange ~ PTC.Status, data = FilteredPYMkdvsControlDESeq,
                   exact = FALSE, alternative = "less")

res$p.value

table(FilteredPYMkdvsControlDESeq$PTC.Status)

#Second comparison: OE vs control#####


# res_counts2 <- results(ddscounts, contrast = c("Sample","PYM OE","Control"))
# 
# PYMOEvsControlDESeq<-as.data.frame(res_counts2)
# PYMOEvsControlDESeq<-na.omit(PYMOEvsControlDESeq)
# plotMA(res_counts2, ylim=c(-7,7))
# 
# glimpse(PYMOEvsControlDESeq)
# PYMOEvsControlDESeq$ENST.ID<-rownames(PYMOEvsControlDESeq)
# write_tsv(PYMOEvsControlDESeq, "DESeqResults/PYMOEvsControlDESeq_tpm0.1.tsv")
#PYM OE all PTC cdf plot####
PYMOEvsControlDESeq<-read_tsv("DESeqResults/PYMOEvsControlDESeq_tpm0.1.tsv")
PYMOEvsControlDESeq$ENST.ID <- sub("\\..*", "", PYMOEvsControlDESeq$ENST.ID)
PYMOEvsControlDESeq<-inner_join(PYMOEvsControlDESeq, PTClist, by = "ENST.ID")

PTCgenes <- getBM(attributes = c("ensembl_gene_id","ensembl_transcript_id",
                                 "transcript_biotype"),
                  filters = "ensembl_transcript_id",
                  values = PYMOEvsControlDESeq$ENST.ID,
                  mart = ensembl100,
                  useCache = FALSE)
table(PTCgenes$transcript_biotype)
names(PTCgenes)[2]<-"ENST.ID"
PYMOEvsControlDESeq<-inner_join(PYMOEvsControlDESeq, PTCgenes, by = "ENST.ID")


PYMOEvsControlDESeq <- PYMOEvsControlDESeq %>%
  dplyr::group_by(ensembl_gene_id) %>%                      # Group by gene ID
  filter(any(PTC.Status == TRUE) & any(PTC.Status == FALSE)) %>%  # Keep only genes with both TRUE & FALSE PTC statuses
  ungroup()

#Plot:PYMOEallPTC####
# calculate pairwise p-values (for example, Wilcoxon test)
p_values <- pairwise.wilcox.test(
  PYMOEvsControlDESeq$log2FoldChange, 
  PYMOEvsControlDESeq$PTC.Status, 
  p.adjust.method = "bonferroni"
)$p.value

# Convert p-values to a readable format for display
# Extract and format p-values correctly based on row and column names
pval_text <- paste("p-values(Wilcoxon test):",
                   paste(format(p_values, digits = 2)),  # "c vs l"
                   sep = "\n")

# Calculate the number of genes in each group
group_counts <- PYMOEvsControlDESeq %>%
  group_by(PTC.Status) %>%
  dplyr::summarize(count = n(), .groups = 'drop')

# Create formatted labels with counts
labels_with_counts <- paste(group_counts$PTC.Status, 
                            " (", group_counts$count, ")", sep = "")
labels_with_counts <-str_replace(labels_with_counts,"TRUE","PTC+")
labels_with_counts <-str_replace(labels_with_counts,"FALSE","PTC-")

# plotting as before
PYMOEvsControlDESeq %>%
  ggplot(aes(log2FoldChange, color = PTC.Status)) +
  geom_vline(aes(xintercept = 0), linewidth = 0.02)+
  geom_hline(aes(yintercept = 0.5), linewidth = 0.02)+
  stat_ecdf(geom = "line", linewidth = 0.8) +
  scale_color_manual(values = c("FALSE" = "grey45", "TRUE" = "red3"), 
                     labels = labels_with_counts) +  # Change legend text
  coord_cartesian(xlim = c(-2, 2)) +
  theme_linedraw() +
  theme(
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    aspect.ratio = 1,
    legend.position = c(0.95, 0.05),         # Move legend to bottom right
    legend.justification = c(1, 0)     # Ensure correct alignment at bottom right
  ) +
  labs(title = "PYM1 OE vs Control",
       y = "Cumulative frequency", x = "log2FC PYM1OE/Control") +
  annotate("text", x = -1.9, y = 0.92, label = pval_text, hjust = 0, size = 4, color = "black")  # Add p-values on top left
res <- wilcox.test(log2FoldChange ~ PTC.Status, data = PYMOEvsControlDESeq,
                   exact = FALSE, alternative = "less")
res$p.value
table(PYMOEvsControlDESeq$PTC.Status)

#PYMOE filtered cdf####
glimpse(PTCgenes_filtered)
glimpse(PYMOEvsControlDESeq)
FilteredPYMOEvsControlDESeq <- PYMOEvsControlDESeq%>%filter(ENST.ID %in% PTCgenes_filtered$ENST.ID)

FilteredPYMOEvsControlDESeq <- FilteredPYMOEvsControlDESeq %>%
  dplyr::group_by(ensembl_gene_id) %>%                      # Group by gene ID
  filter(any(PTC.Status == TRUE) & any(PTC.Status == FALSE)) %>%  # Keep only genes with both TRUE & FALSE PTC statuses
  ungroup()


# calculate pairwise p-values (for example, Wilcoxon test)
p_values <- pairwise.wilcox.test(
  FilteredPYMOEvsControlDESeq$log2FoldChange, 
  FilteredPYMOEvsControlDESeq$PTC.Status, 
  p.adjust.method = "bonferroni"
)$p.value

# Convert p-values to a readable format for display
# Extract and format p-values correctly based on row and column names
pval_text <- paste("p-values(Wilcoxon test):",
                   paste(format(p_values, digits = 2)),  # "c vs l"
                   sep = "\n")

# Calculate the number of genes in each group
group_counts <- FilteredPYMOEvsControlDESeq %>%
  group_by(PTC.Status) %>%
  dplyr::summarize(count = n(), .groups = 'drop')

# Create formatted labels with counts
labels_with_counts <- paste(group_counts$PTC.Status, 
                            " (", group_counts$count, ")", sep = "")
labels_with_counts <-str_replace(labels_with_counts,"TRUE","PTC+")
labels_with_counts <-str_replace(labels_with_counts,"FALSE","PTC-")


# plotting as before
c<-FilteredPYMOEvsControlDESeq %>%
  ggplot(aes(log2FoldChange, color = PTC.Status)) +
  geom_vline(aes(xintercept = 0), linewidth = 0.02)+
  geom_hline(aes(yintercept = 0.5), linewidth = 0.02)+
  stat_ecdf(geom = "line", linewidth = 0.8) +
  scale_color_manual(values = c("FALSE" = "grey45", "TRUE" = "red3"), 
                     labels = labels_with_counts) +  # Change legend text
  coord_cartesian(xlim = c(-2, 2)) +
  theme_linedraw() +
  theme(
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    aspect.ratio = 1,
    legend.position = c(0.98, 0.02),         # Move legend to bottom right
    legend.justification = c(1, 0),    # Ensure correct alignment at bottom right
    legend.title=element_text(size=7), 
    legend.text=element_text(size=5),
    legend.key.size =unit(2, 'mm')
  ) +
  labs(y = "Cumulative frequency", x = "log2FC PYM1 OE/Control") +
  annotate("text", x = -1.9, y = 0.92, label = pval_text, hjust = 0, size = 2, color = "black")  # Add p-values on top left
c+labs(title = "PYM1 OE vs Control",
      caption = "SMG6/7 validated")


res <- wilcox.test(log2FoldChange ~ PTC.Status, data = FilteredPYMOEvsControlDESeq,
                   exact = FALSE, alternative = "less")

res$p.value
table(FilteredPYMOEvsControlDESeq$PTC.Status)
#Third comparison: DelN vs control#####

#res_counts1 <- results(ddscounts, contrast = c("Sample","PYM KD","Control"))
res_counts3 <- results(ddscounts, contrast = c("Sample","PYMΔΝ","Control"))

PYMdelNOEvsControlDESeq<-as.data.frame(res_counts3)
PYMdelNOEvsControlDESeq<-na.omit(PYMdelNOEvsControlDESeq)
plotMA(res_counts3, ylim=c(-7,7))


glimpse(PYMdelNOEvsControlDESeq)
PYMdelNOEvsControlDESeq$ENST.ID <- rownames(PYMdelNOEvsControlDESeq)

write_tsv(PYMdelNOEvsControlDESeq, "DESeqResults/PYMdelNOEvsControlDESeq_tpm0.1.tsv")
#PYM del N all PTC cdf#####

PYMdelNOEvsControlDESeq<-read_tsv("DESeqResults/PYMdelNOEvsControlDESeq_tpm0.1.tsv")
PYMdelNOEvsControlDESeq$ENST.ID <- sub("\\..*", "", PYMdelNOEvsControlDESeq$ENST.ID)
PYMdelNOEvsControlDESeq<-inner_join(PYMdelNOEvsControlDESeq, PTClist, by = "ENST.ID")

PTCgenes <- getBM(attributes = c("ensembl_gene_id","ensembl_transcript_id",
                                 "transcript_biotype"),
                  filters = "ensembl_transcript_id",
                  values = PYMdelNOEvsControlDESeq$ENST.ID,
                  mart = ensembl100,
                  useCache = FALSE)
names(PTCgenes)[2]<-"ENST.ID"
PYMdelNOEvsControlDESeq<-inner_join(PYMdelNOEvsControlDESeq, PTCgenes, by = "ENST.ID")

PYMdelNOEvsControlDESeq <- PYMdelNOEvsControlDESeq %>%
  dplyr::group_by(ensembl_gene_id) %>%                      # Group by gene ID
  filter(any(PTC.Status == TRUE) & any(PTC.Status == FALSE)) %>%  # Keep only genes with both TRUE & FALSE PTC statuses
  ungroup()


# calculate pairwise p-values (for example, Wilcoxon test)
p_values <- pairwise.wilcox.test(
  PYMdelNOEvsControlDESeq$log2FoldChange, 
  PYMdelNOEvsControlDESeq$PTC.Status, 
  p.adjust.method = "bonferroni"
)$p.value

# Convert p-values to a readable format for display
# Extract and format p-values correctly based on row and column names
pval_text <- paste("p-values(Wilcoxon test):",
                   paste(format(p_values, digits = 2)),  # "c vs l"
                   sep = "\n")

# Calculate the number of genes in each group
group_counts <- PYMdelNOEvsControlDESeq %>%
  group_by(PTC.Status) %>%
  dplyr::summarize(count = n(), .groups = 'drop')

# Create formatted labels with counts
labels_with_counts <- paste(group_counts$PTC.Status, 
                            " (", group_counts$count, ")", sep = "")
labels_with_counts <-str_replace(labels_with_counts,"TRUE","PTC+")
labels_with_counts <-str_replace(labels_with_counts,"FALSE","PTC-")


# plotting as before
PYMdelNOEvsControlDESeq %>%
  ggplot(aes(log2FoldChange, color = PTC.Status)) +
  geom_vline(aes(xintercept = 0), linewidth = 0.02)+
  geom_hline(aes(yintercept = 0.5), linewidth = 0.02)+
  stat_ecdf(geom = "line", linewidth = 0.8) +
  scale_color_manual(values = c("FALSE" = "grey45", "TRUE" = "red3"), 
                     labels = labels_with_counts) +  # Change legend text
  coord_cartesian(xlim = c(-2, 2)) +
  theme_linedraw() +
  theme(
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    aspect.ratio = 1,
    legend.position = c(0.95, 0.05),         # Move legend to bottom right
    legend.justification = c(1, 0)     # Ensure correct alignment at bottom right
  ) +
  labs(title = "PYM1 DelN Overexpression vs Control",
       y = "Cumulative frequency", x = "log2FC PYM1 DelN OE/Control") +
  annotate("text", x = -1.9, y = 0.92, label = pval_text, hjust = 0, size = 4, color = "black")  # Add p-values on top left
res <- wilcox.test(log2FoldChange ~ PTC.Status, data = PYMdelNOEvsControlDESeq,
                   exact = FALSE, alternative = "less")
res$p.value
table(PYMdelNOEvsControlDESeq$PTC.Status)

##Filtered DelN CDF####
FilteredPYMdelNOEvsControlDESeq <- PYMdelNOEvsControlDESeq%>%filter(ENST.ID %in% PTCgenes_filtered$ENST.ID)
FilteredPYMdelNOEvsControlDESeq <- FilteredPYMdelNOEvsControlDESeq %>%
  dplyr::group_by(ensembl_gene_id) %>%                      # Group by gene ID
  filter(any(PTC.Status == TRUE) & any(PTC.Status == FALSE)) %>%  # Keep only genes with both TRUE & FALSE PTC statuses
  ungroup()


# calculate pairwise p-values (for example, Wilcoxon test)
p_values <- pairwise.wilcox.test(
  FilteredPYMdelNOEvsControlDESeq$log2FoldChange, 
  FilteredPYMdelNOEvsControlDESeq$PTC.Status, 
  p.adjust.method = "bonferroni"
)$p.value

# Convert p-values to a readable format for display
# Extract and format p-values correctly based on row and column names
pval_text <- paste("p-values(Wilcoxon test):",
                   paste(format(p_values, digits = 2)),  # "c vs l"
                   sep = "\n")

# Calculate the number of genes in each group
group_counts <- FilteredPYMdelNOEvsControlDESeq %>%
  group_by(PTC.Status) %>%
  dplyr::summarize(count = n(), .groups = 'drop')

# Create formatted labels with counts
labels_with_counts <- paste(group_counts$PTC.Status, 
                            " (", group_counts$count, ")", sep = "")
labels_with_counts <-str_replace(labels_with_counts,"TRUE","PTC+")
labels_with_counts <-str_replace(labels_with_counts,"FALSE","PTC-")


# plotting as before
d<-FilteredPYMdelNOEvsControlDESeq %>%
  ggplot(aes(log2FoldChange, color = PTC.Status)) +
  geom_vline(aes(xintercept = 0), linewidth = 0.02)+
  geom_hline(aes(yintercept = 0.5), linewidth = 0.02)+
  stat_ecdf(geom = "line", linewidth = 0.8) +
  scale_color_manual(values = c("FALSE" = "grey45", "TRUE" = "red3"), 
                     labels = labels_with_counts) +  # Change legend text
  coord_cartesian(xlim = c(-2, 2)) +
  theme_linedraw() +
  theme(
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    aspect.ratio = 1,
    legend.position = c(0.98, 0.02),         # Move legend to bottom right
    legend.justification = c(1, 0),    # Ensure correct alignment at bottom right
    legend.title=element_text(size=7), 
    legend.text=element_text(size=5),
    legend.key.size =unit(2, 'mm')
  ) +
  labs(y = "Cumulative frequency", x = "log2FC PYM1 DelN OE/Control") +
  annotate("text", x = -1.9, y = 0.92, label = pval_text, hjust = 0, size = 2, color = "black")  # Add p-values on top left

d+labs(title = "PYM1 DelN Overexpression vs Control", caption = "SMG6/7 validated")
res <- wilcox.test(log2FoldChange ~ PTC.Status, data = FilteredPYMdelNOEvsControlDESeq,
                   exact = FALSE, alternative = "less")

res$p.value

table(FilteredPYMdelNOEvsControlDESeq$PTC.Status)
#####




#QPCR plot#####


qPCRdata<-read_csv("qPCRdata.csv")

#analyzing
table(qPCRdata$Target)

##Replacing PTC- with NOR
qPCRdata$Target<-str_replace_all(qPCRdata$Target,"PTC-", "NOR")


qPCRdata<-qPCRdata%>%filter(!Target %in% c("B-actin","TBP"))
qPCRdata%>%filter(Target == "GAPDH")

qPCRdata$Gene<-str_split(qPCRdata$Target,"\\_",simplify = T)[,1]
table(qPCRdata$Gene)

qPCRdata%>%filter(Gene=="RPS9")
print(n = 81,
      qPCRdata%>%filter(Gene=="RPS9")%>%arrange(File,Sample,Target))

qPCRdata<-qPCRdata%>%mutate(Replicate = str_split(File, "Q", simplify = T)[,1])
print(n = 181,
      qPCRdata%>%group_by(File,Sample,Target)%>%dplyr::count())

print(n = 181,
      qPCRdata%>%group_by(Replicate,Sample,Target)%>%dplyr::count())


###Data wrangling
qPCRdata<-qPCRdata%>% mutate(NewColumn = File)



qPCRdata<-qPCRdata%>% mutate(NewColumn = case_when(Replicate == "R2" ~ "R2",
                                                   Replicate == "R3" ~ "R3",
                                                   TRUE ~ NewColumn))
table(qPCRdata$NewColumn)
table(qPCRdata$Replicate)
qPCRdata<-qPCRdata%>% mutate(RepSample = paste0(NewColumn, "_",Sample))

qPCRdata<-qPCRdata%>% mutate(FileSample = paste0(File, "_",Sample))

qPCRdata%>%filter(File=="R1Q1" & Target =="GAPDH")

library(Rmisc)
CqSummary <- summarySE(qPCRdata, measurevar="Cq", groupvars=c("FileSample","Target"))
#CqSummary <- summarySE(qPCRdata, measurevar="Cq", groupvars=c("RepSample","Target"))

##Plotting GAPDH Cq values
CqSummary%>%filter(Target=="RSRC2_PTC+")%>%ggplot(aes(FileSample, Cq))+geom_col()+
  geom_errorbar(aes(ymin=Cq-sd,ymax=Cq+sd),width=.2)+
  coord_cartesian(ylim = c(18,27)) + 
  theme(axis.text = element_text(hjust = 1,angle = 30))

CqSummary%>%filter(Target=="RSRC2_NOR")%>%ggplot(aes(FileSample, Cq))+geom_col()+
  geom_errorbar(aes(ymin=Cq-sd,ymax=Cq+sd),width=.2)+
  coord_cartesian(ylim = c(16,27)) + 
  theme(axis.text = element_text(hjust = 1,angle = 30))


delcqdf<-tibble()
for(i in unique(CqSummary$FileSample)){
  print(i)
  filteredCqsummary<-CqSummary%>%filter(FileSample==i)%>%mutate(
    Delcq = Cq - filter(CqSummary,FileSample==i&Target=="GAPDH")[,4] )
  delcqdf<-rbind(delcqdf,filteredCqsummary)
}
delcqdf<-mutate(delcqdf, Rel.levels = 2^-Delcq)


#Normalizing to respective Nluc
#create a column that says Nluc instead of Nluc_2 or Nluc_3
delcqdf$Sampletype<-str_split(delcqdf$FileSample,"_",simplify = T)[,2]
#Create a column that says replicate number instead of Nluc_2 or Nluc_3
delcqdf$RepNo<-str_split(delcqdf$FileSample,"_",simplify = T)[,1]

finalNormdf<-tibble()
for(j in unique(delcqdf$RepNo)){
  print(j)
  for(i in unique(delcqdf$Target)){
    print(i)
    finaldf<-delcqdf%>%filter(RepNo==j&Target==i)
    finaldf<-finaldf%>%mutate(
      Rel.Norm.Levels = Rel.levels/filter(finaldf,Sampletype=="siNC"&Target==i)[,9] )
    finalNormdf<-rbind(finalNormdf,finaldf)
  }}

names(finaldf)

###modified plot
library(wesanderson)
library(RColorBrewer)

table(finalNormdf$RepNo)
table(finalNormdf$Sampletype)
table(finalNormdf$Target)
str_replace(finalNormdf$Target, "NOR", "PTC-")
finalNormdf<-finalNormdf%>%mutate(Target = str_replace(Target, "NOR", "PTC-"))
finalNormdf%>%
  ggplot(aes(Target,Rel.Norm.Levels,fill=factor(Sampletype, levels = c('siNC', 'si PYM1', 'si eIF4A3', 'si UPF1'))))+
  stat_summary(fun = mean, geom = "bar", position = "dodge", width= 0.6)+
  stat_summary(fun.data = mean_cl_boot, geom = "errorbar",width = 0.3, position = position_dodge(width = 0.6))+
  scale_fill_manual(values = wes_palette("GrandBudapest2"))+
  #geom_point(position = position_dodge(width = 0.6), size = 0.5)+
  theme_classic()+
  labs(title="Normalized RNA levels - NMD factor depletion",
       subtitle = "3 Bio replicates, normalized to GAPDH and siNC")+
  guides(fill=guide_legend(title="Sample"))+
  theme(axis.text = element_text(hjust = 1,angle = 30),
        legend.position = c(0.01, 1),
        legend.justification = c(0, 1))

finalNormdf<-finalNormdf%>%mutate(Targets=str_split(Target,"\\_", simplify = T)[,1])
finalNormdf<-finalNormdf%>%mutate(PTC_status=str_split(Target,"\\_", simplify = T)[,2])

b<-finalNormdf%>%filter(Sampletype!="si UPF1")%>%filter(!Target %in%c("RPS9_PTC+","RPS9_PTC-"))%>%
  ggplot(aes(factor(Target, levels = c("GAPDH","GAS5","ZFAS1","RPS9_PTC-","RPS9_PTC+","RSRC2_PTC-","RSRC2_PTC+",
                                       "SRSF3_PTC-","SRSF3_PTC+","SRSF2_PTC-","SRSF2_PTC+","SRSF6_PTC+")),
             Rel.Norm.Levels,
             fill=factor(Sampletype, levels = c('siNC', 'si PYM1', 'si eIF4A3', 'si UPF1'))))+
  geom_hline(yintercept = 1, size = 0.5, linetype =3, color = "grey")+
  stat_summary(fun = mean, geom = "bar", position = "dodge", width= 0.6)+
  stat_summary(fun.data = mean_se, geom = "errorbar",width = 0.3, position = position_dodge(width = 0.6))+
  scale_fill_manual(values = c("grey60","#DE7D57","red3"))+
  #geom_point(position = position_dodge(width = 0.6), size = 0.5)+
  theme_classic()+
  labs(#title="Normalized RNA levels - NMD factor depletion",
       #subtitle = "3 Bio replicates, normalized to GAPDH and siNC",
       x = "Targets",y="Relative RNA levels")+
  guides(fill=guide_legend(title="Sample"))+
  theme(panel.grid.major = element_line(linetype = "blank"), 
        panel.grid.minor = element_line(linetype = "blank"),
        axis.text.x = element_text(size = 7,angle = 60,hjust=1),
        legend.direction = "horizontal",legend.position = "top",
        legend.justification = c(0, 1),
        legend.title=element_text(size=9), 
        legend.text=element_text(size=7),
        legend.key.size =unit(3, 'mm'),
        axis.line = element_line(linewidth  = 0.25),
        aspect.ratio = 0.6)
b


##Patchwork
library(patchwork)
layout <- c(
  area(1,1,2,2),area(1, 3, 2,4),
  area(3,1,4,2),area(3,3,4,4),area(3,6,4,6),
  area(5,1,8,3),area(5,4,8,6))
plot(layout)
a+b+c+d+plot_layout(design = layout)+
  plot_annotation(title = 'FIGURE 5', tag_levels = "A")  &
  theme(plot.tag = element_text(face = 'bold')) 
ggsave("S5patch.pdf",
       a+b+c+d+plot_layout(design = layout)+
         plot_annotation(title = 'FIGURE 5', tag_levels = "A")  &
         theme(plot.tag = element_text(face = 'bold')),
       width = 8.5, 
       height = 11, 
       units = "in")


##PLot of high exon numbers
PYMkdvsControlDESeq<-read_tsv("DESeqResults/PYMkdvsControlDESeq_tpm0.1.tsv")
PYMkdvsControlDESeq$ensembl_transcript_id<-str_split(PYMkdvsControlDESeq$ensembl_transcript_id, "\\.", simplify = T)[,1]
ExonNumbers <- getBM(attributes = c("ensembl_gene_id","ensembl_transcript_id",
                                 "rank","tmhmm","signalp"),
                  filters = "ensembl_transcript_id",
                  values = PYMkdvsControlDESeq$ensembl_transcript_id,
                  mart = ensembl100,
                  useCache = FALSE)
table(ExonNumbers$tmhmm)
table(ExonNumbers$signalp)

ExonNumbers<-ExonNumbers%>%arrange(ensembl_transcript_id,desc(rank))%>%group_by(ensembl_transcript_id)%>%dplyr::slice(1)
ExonNumbers%>%ggplot(aes(rank))+
  geom_bar()+coord_cartesian(xlim = c(0,25))
summary(ExonNumbers$rank)

names(ExonNumbers)[2]<-"ensembl_transcript_id"

PYMkdvsControlDESeq<-inner_join(PYMkdvsControlDESeq,ExonNumbers)
PYMkdvsControlDESeq<-PYMkdvsControlDESeq%>%filter(tmhmm != "TMhelix")
PYMkdvsControlDESeq<-PYMkdvsControlDESeq%>%filter(signalp != "SignalP-noTM")
PYMkdvsControlDESeq<-PYMkdvsControlDESeq%>%filter(signalp != "SignalP-TM")

table(PYMkdvsControlDESeq$signalp)
summary(PYMkdvsControlDESeq$rank)

PYMkdvsControlDESeq<-PYMkdvsControlDESeq%>%mutate(ExonNumber=case_when(rank>13 ~ "High",
                                                                       rank<6 ~ "Low",
                                                                       TRUE ~ "Others"))
table(PYMkdvsControlDESeq$ExonNumber)
PYMkdvsControlDESeq<-PYMkdvsControlDESeq%>%mutate(ExonNumber2=case_when(rank>12 ~ "High",
                                                                       #rank<4 ~ "Low",
                                                                       TRUE ~ "Others"))

table(PYMkdvsControlDESeq$ExonNumber2)



p_values <- pairwise.wilcox.test(
  PYMkdvsControlDESeq$log2FoldChange, 
  PYMkdvsControlDESeq$ExonNumber2, 
  p.adjust.method = "bonferroni"
)$p.value

# Convert p-values to a readable format for display
# Extract and format p-values correctly based on row and column names
pval_text <- paste("p-value (Wilcoxon test):",
                   paste(format(p_values, digits = 2)),  # "c vs l"
                   sep = "\n")

# Calculate the number of genes in each group
group_counts <- PYMkdvsControlDESeq %>%
  group_by(ExonNumber2) %>%
  dplyr::summarize(count = n(), .groups = 'drop')

# Create formatted labels with counts
labels_with_counts <- paste(group_counts$ExonNumber2, 
                            " (", group_counts$count, ")", sep = "")


s2<-PYMkdvsControlDESeq%>%ggplot(aes(log2FoldChange, color = ExonNumber2)) +
  geom_vline(aes(xintercept = 0), linewidth = 0.02)+
  geom_hline(aes(yintercept = 0.5), linewidth = 0.02)+
  stat_ecdf(geom = "line", linewidth = 0.4)+
  coord_cartesian(xlim = c(-2,2))+
  scale_color_manual(values = c("Others" = "grey45", "High" = "#0f7ba2"), 
                     labels = labels_with_counts) +  # Change legend text+
  theme_linedraw() +
  theme(axis.text = element_text(size=4),axis.title = element_text(size=6),
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    aspect.ratio = 1,
    legend.position = c(0.98, 0.02),         # Move legend to bottom right
    legend.justification = c(1, 0),    # Ensure correct alignment at bottom right
    legend.title=element_text(size=6), 
    legend.text=element_text(size=4),
    legend.key=element_blank(),
    legend.box.background = element_blank(),
    legend.background = element_blank(),
    legend.key.size =unit(2, 'mm')
  ) +
  labs(y = "Cumulative frequency", x = "log2FC PYM1 kd/Control") +
  annotate("text", x = -1.9, y = 0.92, label = pval_text, hjust = 0, size = 2, color = "black")  # Add p-values on top left

# Convert p-values to a readable format for display
# Extract and format p-values correctly based on row and column names
p_values <- pairwise.wilcox.test(
  PYMkdvsControlDESeq$log2FoldChange, 
  PYMkdvsControlDESeq$ExonNumber, 
  p.adjust.method = "bonferroni"
)$p.value

# Format p-value text
pval_text <- paste("p values (Wilcoxon test):",
                   paste("High vs Low =", format(p_values["Low", "High"], scientific = TRUE, digits = 2)),
                   paste("High vs Others =", format(p_values["Others", "High"], scientific = TRUE, digits = 2)),
                   paste("Others vs Low =", format(p_values["Others", "Low"], scientific = TRUE, digits = 2)),
                   sep = "\n")


# Calculate the number of genes in each group
group_counts <- PYMkdvsControlDESeq %>%
  group_by(ExonNumber) %>%
  dplyr::summarize(count = n(), .groups = 'drop')

# Create formatted labels with counts
labels_with_counts <- paste(group_counts$ExonNumber, 
                            " (", group_counts$count, ")", sep = "")
Egypt = list(c("#dd5129", "#0f7ba2", "#43b284", "#fab255"), c(1, 2, 3, 4), colorblind=TRUE)

s1<-PYMkdvsControlDESeq%>%ggplot(aes(log2FoldChange, color = ExonNumber)) +
  geom_vline(aes(xintercept = 0), linewidth = 0.02)+
  geom_hline(aes(yintercept = 0.5), linewidth = 0.02)+
  stat_ecdf(geom = "line", linewidth = 0.4)+
  coord_cartesian(xlim = c(-2,2))+
  scale_color_manual(values = c("Others" = "grey70", "High" = "#0f7ba2", "Low" = "#43b284"), 
                     labels = labels_with_counts) +  # Change legend text+
  theme_linedraw() +
  theme(axis.text = element_text(size=4),axis.title = element_text(size=6),
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    aspect.ratio = 1,
    legend.position = c(0.98, 0.02),   # Move legend to bottom right
    legend.justification = c(1, 0),    # Ensure correct alignment at bottom right
    legend.title=element_text(size=6), 
    legend.text=element_text(size=4),
    legend.key=element_blank(),
    legend.box.background = element_blank(),
    legend.background = element_blank(),
    legend.key.size =unit(2, 'mm')
  ) +
  labs(y = "Cumulative frequency", x = "log2FC PYM1 kd/Control") +
  annotate("text", x = -1.9, y = 0.92, label = pval_text, hjust = 0, size = 2, color = "black")  # Add p-values on top left

s1+s2
