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


#Plot:PYMkdvsControl_allPTC####
# calculate pairwise p-values (for example, Wilcoxon test)
p_values <- pairwise.wilcox.test(
  PYMkdvsControlDESeq$log2FoldChange, 
  PYMkdvsControlDESeq$PTC.Status, 
  p.adjust.method = "bonferroni"
)$p.value

wilcox.test(log2FoldChange ~ PTC.Status, data = PYMkdvsControlDESeq,
                           exact = FALSE)

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
s1<-PYMkdvsControlDESeq %>%
  ggplot(aes(log2FoldChange, color = PTC.Status)) +
  geom_vline(aes(xintercept = 0), linewidth = 0.02)+
  geom_hline(aes(yintercept = 0.5), linewidth = 0.02)+
  stat_ecdf(geom = "line", linewidth = 0.8) +
  scale_color_manual(values = c("FALSE" = "grey45", "TRUE" = "red3"), 
                     labels = labels_with_counts) +  # Change legend text
  coord_cartesian(xlim = c(-2, 2)) +
  theme_linedraw() +
  theme(axis.title = element_text(size = 5),axis.text = element_text(size = 3),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        aspect.ratio = 1,
        legend.position = c(0.98, 0.02),         # Move legend to bottom right
        legend.justification = c(1, 0),    # Ensure correct alignment at bottom right
        legend.title=element_text(size=7), 
        legend.text=element_text(size=5),
        legend.key.size =unit(2, 'mm')
  )+
  labs(#title = "PYM1 kd vs Control",
       y = "Cumulative frequency", x = "log2FC PYM1kd/Control") +
  annotate("text", x = -1.9, y = 0.92, label = pval_text, hjust = 0, size = 2, color = "black")  # Add p-values on top left

s1+ labs(title = "PYM1 kd vs Control")


##Comparing with CASC3 and SMGko
CASC3deseq<-read_tsv("DESeqResults/CASC3KOdeseqresults.tsv")
SMG67deseq<-read_tsv("DESeqResults/SMG67kdvsControlDESeq.tsv")
SMG67deseq <- SMG67deseq[, c(2:7, 1)]
PYMdeseq<-read_tsv("DESeqResults/PYMkdvsControlDESeq_tpm0.1.tsv")
PTClist <- read.delim("~/OneDrive - The Ohio State University/Singhlab/Bioinfo/CASC3/ENST_PTC-EPI-TFG.txt")

names(PYMkdvsControlDESeq)[7] <- "ENST.ID"
PYMkdvsControlDESeq$ENST.ID <- sub("\\..*", "", PYMkdvsControlDESeq$ENST.ID)
length(intersect(PYMkdvsControlDESeq$ENST.ID,PTClist$ENST.ID))

PYMkdvsControlDESeq<-inner_join(PYMkdvsControlDESeq, PTClist, by = "ENST.ID")

Foldchangedf <- data.frame()  # Initialize empty result data frame

for (name in c("PYMdeseq", "CASC3deseq", "SMG67deseq")) {
  tmp <- get(name)  # Retrieve the data frame by its name as a string
  names(tmp)[7] <- "ENST.ID"
  tmp$ENST.ID <- sub("\\..*", "", tmp$ENST.ID)
  tmp <- inner_join(tmp, PTClist, by = "ENST.ID")
  tmp$comparison <- str_remove(name, "deseq")
  Foldchangedf <- rbind(Foldchangedf, tmp)}
table(Foldchangedf$comparison)
#plot:
Foldchangedf<-Foldchangedf%>%
  mutate(PTC_status = case_when(PTC.Status == T ~ "PTC+",
                                 PTC.Status == F ~ "PTC-"))
Foldchangedf$comparison<-factor(Foldchangedf$comparison, levels=c("PYM","CASC3","SMG67"))
s2<-Foldchangedf%>%
  ggplot(aes(PTC_status, log2FoldChange, fill = PTC_status))+
  geom_hline(yintercept = 0,size = 0.05)+
  geom_boxplot(outliers = F)+
  facet_wrap(~comparison)+
  scale_fill_manual(values = c("PTC-" = "grey45", "PTC+" = "red3")) +
  theme_linedraw() +
  coord_cartesian(ylim = c(-3.2,4.5))+
  theme( axis.title = element_text(size = 5),axis.text = element_text(size = 3),
         panel.grid.major = element_blank(),
         panel.grid.minor = element_blank(),
         aspect.ratio = 3.2,
         legend.position = "none",
         strip.text = element_text(color="black",size =4, angle = 0, hjust =0),
         strip.background = element_blank()) +
  labs( y = "log2FC",
        x = "PTC status" ) +
  stat_compare_means(comparisons = list(c("PTC-", "PTC+")),
                     method = "wilcox.test",
                     label = "p.format",
                     label.y = 3.58, size = 1,tip.length = 0.005) # Adjust as needed for positioning
s2
Foldchangedf %>%
  +     dplyr::count(comparison, PTC_status)

# Initialize results dataframe
comparison_results <- data.frame(
  comparison = character(),
  test = character(),
  p_value = numeric(),
  statistic = numeric(),
  n_PTCpos = integer(),
  n_PTCneg = integer(),
  stringsAsFactors = FALSE
)

# Loop over each unique comparison group
for (comp in unique(Foldchangedf$comparison)) {
  print(comp)
  subdf <- Foldchangedf %>% filter(comparison == comp)
  # Grouped values
  ptc_plus  <- subdf %>% filter(PTC_status == "PTC+") %>% pull(log2FoldChange)
  ptc_minus <- subdf %>% filter(PTC_status == "PTC-") %>% pull(log2FoldChange)
  # Perform Wilcoxon test (non-parametric)
  test <- wilcox.test(ptc_plus, ptc_minus, alternative = "two.sided")
  # Append results
  comparison_results <- rbind(comparison_results, data.frame(
    comparison = comp,
    test = "Wilcoxon",
    p_value = test$p.value,
    statistic = unname(test$statistic),
    n_PTCpos = length(ptc_plus),
    n_PTCneg = length(ptc_minus)
  ))
}

# Apply Bonferroni correction
comparison_results$p_adj <- p.adjust(comparison_results$p_value, method = "bonferroni")

# View results
print(comparison_results)



#Second comparison: OE vs control#####

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
wilcox.test(log2FoldChange ~ PTC.Status, data = PYMOEvsControlDESeq,
            exact = FALSE)

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
s5<-PYMOEvsControlDESeq %>%
  ggplot(aes(log2FoldChange, color = PTC.Status)) +
  geom_vline(aes(xintercept = 0), linewidth = 0.02)+
  geom_hline(aes(yintercept = 0.5), linewidth = 0.02)+
  stat_ecdf(geom = "line", linewidth = 0.8) +
  scale_color_manual(values = c("FALSE" = "grey45", "TRUE" = "red3"), 
                     labels = labels_with_counts) +  # Change legend text
  coord_cartesian(xlim = c(-2, 2)) +
  theme_linedraw() +
  theme(axis.title = element_text(size = 5),axis.text = element_text(size = 3),
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    aspect.ratio = 1,
    legend.position = c(0.98, 0.02),         # Move legend to bottom right
    legend.justification = c(1, 0),    # Ensure correct alignment at bottom right
    legend.title=element_text(size=7), 
    legend.text=element_text(size=5),
    legend.key.size =unit(2, 'mm')
  ) +
  labs(#title = "PYM1 OE vs Control",
       y = "Cumulative frequency", x = "log2FC PYM1OE/Control") +
  annotate("text", x = -1.9, y = 0.92, label = pval_text, hjust = 0, size = 2, color = "black")  # Add p-values on top left
s5

#PYMOE filtered cdf####

SMG67kdvsControlDESeq<-read_tsv("~/OneDrive - The Ohio State University/Singhlab/Bioinfo/CASC3_RNPS1_SMG67_comparison/SMG67kdvsControlDESeq.tsv")

SMG67kdvsControlDESeq<-inner_join(SMG67kdvsControlDESeq,PTClist, by = "ENST.ID")
SMG67kdvsControlDESeq<-na.omit(SMG67kdvsControlDESeq)
table(SMG67kdvsControlDESeq$PTC.Status)

#Filtering for upregulated PTC+ and unchanged PTC-
PTCgenes_filtered <- SMG67kdvsControlDESeq %>%
  filter((PTC.Status == "FALSE") & log2FoldChange >= -0.58 & log2FoldChange <= 0.58 |
           (PTC.Status == "TRUE" & log2FoldChange > 0.58))
table(PTCgenes_filtered$PTC.Status)
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

wilcox.test(log2FoldChange ~ PTC.Status, data = FilteredPYMOEvsControlDESeq,
            exact = FALSE)

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
s3<-FilteredPYMOEvsControlDESeq %>%
  ggplot(aes(log2FoldChange, color = PTC.Status)) +
  geom_vline(aes(xintercept = 0), linewidth = 0.02)+
  geom_hline(aes(yintercept = 0.5), linewidth = 0.02)+
  stat_ecdf(geom = "line", linewidth = 0.8) +
  scale_color_manual(values = c("FALSE" = "grey45", "TRUE" = "red3"), 
                     labels = labels_with_counts) +  # Change legend text
  coord_cartesian(xlim = c(-2, 2)) +
  theme_linedraw() +
  theme(axis.title = element_text(size = 5),axis.text = element_text(size = 3),
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
s3+labs(title = "PYM1 OE vs Control",
       caption = "SMG6/7 validated")

#Third comparison: DelN vs control#####
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
wilcox.test(log2FoldChange ~ PTC.Status, data = PYMdelNOEvsControlDESeq,
            exact = FALSE)
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
s6<-PYMdelNOEvsControlDESeq %>%
  ggplot(aes(log2FoldChange, color = PTC.Status)) +
  geom_vline(aes(xintercept = 0), linewidth = 0.02)+
  geom_hline(aes(yintercept = 0.5), linewidth = 0.02)+
  stat_ecdf(geom = "line", linewidth = 0.8) +
  scale_color_manual(values = c("FALSE" = "grey45", "TRUE" = "red3"), 
                     labels = labels_with_counts) +  # Change legend text
  coord_cartesian(xlim = c(-2, 2)) +
  theme_linedraw() +
  theme(axis.title = element_text(size = 5),axis.text = element_text(size = 3),
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    aspect.ratio = 1,
    legend.position = c(0.98, 0.02),         # Move legend to bottom right
    legend.justification = c(1, 0),    # Ensure correct alignment at bottom right
    legend.title=element_text(size=7), 
    legend.text=element_text(size=5),
    legend.key.size =unit(2, 'mm')
  ) +
  labs(#title = "PYM1 DelN Overexpression vs Control",
       y = "Cumulative frequency", x = "log2FC PYM1 DelN OE/Control") +
  annotate("text", x = -1.9, y = 0.92, label = pval_text, hjust = 0, size = 2, color = "black")  # Add p-values on top left
res <- wilcox.test(log2FoldChange ~ PTC.Status, data = PYMdelNOEvsControlDESeq,
                   exact = FALSE)
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
wilcox.test(log2FoldChange ~ PTC.Status, data = FilteredPYMdelNOEvsControlDESeq,
            exact = FALSE)
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
s4<-FilteredPYMdelNOEvsControlDESeq %>%
  ggplot(aes(log2FoldChange, color = PTC.Status)) +
  geom_vline(aes(xintercept = 0), linewidth = 0.02)+
  geom_hline(aes(yintercept = 0.5), linewidth = 0.02)+
  stat_ecdf(geom = "line", linewidth = 0.8) +
  scale_color_manual(values = c("FALSE" = "grey45", "TRUE" = "red3"), 
                     labels = labels_with_counts) +  # Change legend text
  coord_cartesian(xlim = c(-2, 2)) +
  theme_linedraw() +
  theme(axis.title = element_text(size = 5),axis.text = element_text(size = 3),
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

s4+labs(title = "PYM1 DelN Overexpression vs Control", caption = "SMG6/7 validated")

#patchwork
library(patchwork)
layout <- c(
  area(1,1,2,2),area(1, 3, 2,4),
  area(3,1,4,2),area(3,3,4,4),area(3,5,4,6),
  area(3,7,4,8),area(5,4,8,6))
plot(layout)

ggsave("S5patch2.pdf",
       s1+s2+s3+s4+s5+s6+
         plot_layout(design = layout)+
         plot_annotation(title = 'FIGURE S10', tag_levels = "a")  &
         theme(plot.tag = element_text(face = 'bold')),
       width = 8.5, 
       height = 11, 
       units = "in")

  