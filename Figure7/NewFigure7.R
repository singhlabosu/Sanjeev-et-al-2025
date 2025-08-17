###Figure7 plots

library(tidyverse)
library(ggrepel)
library(patchwork)
library(tximport)
library(DESeq2)
library(ggsignif)

setwd("~/OneDrive - The Ohio State University/Singhlab/Bioinfo/PYMpaperFiguresUpdate/Figure7/")


#First: Comparing with SMG6/7kd#####

SMG67kdvsControlDESeq<-read_tsv("~/OneDrive - The Ohio State University/Singhlab/Bioinfo/CASC3_RNPS1_SMG67_comparison/SMG67kdvsControlDESeq.tsv")
PTClist <- read.delim("~/OneDrive - The Ohio State University/Singhlab/Bioinfo/CASC3/ENST_PTC-EPI-TFG.txt")
SMG67kdvsControlDESeq<-inner_join(SMG67kdvsControlDESeq,PTClist, by = "ENST.ID")
SMG67kdvsControlDESeq<-na.omit(SMG67kdvsControlDESeq)
table(SMG67kdvsControlDESeq$PTC.Status)

##Filtering for upregulated PTC+ and unchanged PTC-
PTCgenes_filtered <- SMG67kdvsControlDESeq %>%
  filter((PTC.Status == "FALSE") & log2FoldChange >= -0.58 & log2FoldChange <= 0.58 |
           (PTC.Status == "TRUE" & log2FoldChange > 0.58))
table(PTCgenes_filtered$PTC.Status)

#Getting transcript info 
library(biomaRt)
listMarts(host='https://apr2020.archive.ensembl.org')
ensembl100=useMart(host='https://apr2020.archive.ensembl.org', biomart='ENSEMBL_MART_ENSEMBL', dataset='hsapiens_gene_ensembl')
attributes <- listAttributes(ensembl100)
filters <- listFilters(ensembl100)

PTCgenes <- getBM(attributes = c("ensembl_gene_id","ensembl_transcript_id",
                                 "transcript_biotype"),
                  filters = "ensembl_transcript_id",
                  values = PTCgenes_filtered$ENST.ID,
                  mart = ensembl100,
                  useCache = FALSE)

names(PTCgenes)[2]<-"ENST.ID"
table(PTCgenes$transcript_biotype)
PTCgenes_filtered<-inner_join(PTCgenes_filtered,PTCgenes,by = "ENST.ID")



##Plotting each DESeq result in loop

library(patchwork)

layout <- c(
  area(1,1,2,2),area(1, 3, 2,4),
  area(1,5,2,6),area(3,1,4,2),area(3,3,4,4),
  area(5,1,8,3),area(5,4,9,7))
plot(layout)

# Initialize stat test results dataframe
comparison_results <- data.frame(
  Group = character(),
  Comparison = character(),
  n_set1 = integer(),
  n_set2 = integer(),
  test = character(),
  p_value = numeric(),
  W_statistic = numeric(),
  Alternate_Hypothesis = character(),
  stringsAsFactors = FALSE)


#NMD plots in a loop
files<-list.files("Tx_DESeq/")
n=1
for (id in files) {
  print(id)
df<-read_tsv(paste0("Tx_DESeq/", id),show_col_types = FALSE)
name<-str_split(id,"vs", simplify = T)[1]
print(name)
names(df)[7]<-"ENST.ID"
df<-inner_join(df,PTCgenes_filtered[,c(1,8:10)])
df <- df %>%
  dplyr::group_by(ensembl_gene_id) %>%                      # Group by gene ID
  filter(any(PTC.Status == TRUE) & any(PTC.Status == FALSE)) %>%  # Keep only genes with both TRUE & FALSE PTC statuses
  ungroup()
# calculate pairwise p-values (for example, Wilcoxon test)
test<- wilcox.test(log2FoldChange ~ PTC.Status, 
                   data = df,
  alternative =  "less")
p_values<-test$p.value

# Convert p-values to a readable format for display
# Extract and format p-values correctly based on row and column names
pval_text <- paste("p values (Wilcoxon test):",
                   paste(format(p_values, scientific = TRUE, digits = 2)),  # "PTC+ vs PTC-"
                   sep = "\n")

# Calculate the number of genes in each group
group_counts <- df %>%
  group_by(PTC.Status) %>%
  dplyr::summarize(count = n(), .groups = 'drop')

# Create formatted labels with counts
labels_with_counts <- paste(group_counts$PTC.Status, 
                            " (", group_counts$count, ")", sep = "")
labels_with_counts <-str_replace(labels_with_counts,"TRUE","PTC+")
labels_with_counts <-str_replace(labels_with_counts,"FALSE","PTC-")

# Define x-axis limits and annotation position for specific files
if (id %in% c("DcapvsNlucDESeq.tsv", "WcapvsNlucDESeq.tsv")) {
  x_limits <- c(-2, 2)   # Wider range for specific files
  annotate_x <- -1.9 }      # Adjust annotation position accordingly
  else {x_limits <- c(-4, 4)    # Default range
        annotate_x <- -3.9     # Default annotation position
  control<-"/Control"}     

if (id == "DcapvsNlucDESeq.tsv") {  
  name <- "DENV Capsid" 
  control<-"/Nluc"} 
if (id == "WcapvsNlucDESeq.tsv") {  
  name <- "WNV Capsid" 
  control<-"/Nluc"} 
if (id == "ZIKV(MR766)vsControlDESeq.tsv") {  
  x_limits <- c(-3, 3)   # Wider range for specific files
  annotate_x <- -2.9 }
print(name)
print(control)
# plotting as before
assign(paste0("n",n) , df %>%
  ggplot(aes(log2FoldChange, color = PTC.Status)) +
  geom_vline(aes(xintercept = 0), linewidth = 0.02)+
  geom_hline(aes(yintercept = 0.5), linewidth = 0.02)+
  stat_ecdf(geom = "line", linewidth = 0.4) +
  scale_color_manual(values = c("FALSE" = "grey45", "TRUE" = "red3"), 
                     labels = labels_with_counts) +  # Change legend text
  coord_cartesian(xlim = x_limits) +
  theme_linedraw() +
  theme(axis.text = element_text(size = 5),axis.title = element_text(size = 7),
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    aspect.ratio = 1,
    legend.position = c(0.95, 0.05),         # Move legend to bottom right
    legend.justification = c(1, 0),     # Ensure correct alignment at bottom right
    legend.title=element_text(size=6), 
    legend.text=element_text(size=4),
    legend.key=element_blank(),
    legend.box.background = element_blank(),
    legend.background = element_blank(),
    legend.key.size =unit(2, 'mm'),
    plot.background = element_blank(),
    plot.margin = unit(c(0, 0, 0, 0), "cm")) +
  labs(#title = paste0( name," vs Control"),
       y = "Cumulative frequency", x = paste0("log2FC ",name,control)) +
  annotate("text", x = annotate_x, y = 0.92, label = pval_text, hjust = 0, size = 1.5, color = "black")  # Add p-values on top left
)

comparison_results <- rbind(comparison_results, data.frame(
  Group = str_split(id, "\\.", simplify = T)[1],
  Comparison = "PTC+ vs PTC-",
  n_set1 = table(df$PTC.Status)[[2]],
  n_set2 = table(df$PTC.Status)[[1]],
  test = "Wilcoxon rank sum test with continuity correction",
  p_value = test$p.value,
  W_statistic = unname(test$statistic),
  Alternate_Hypothesis = "true location shift is greater than 0"
))
n=n+1}

print(comparison_results)

ggsave("NMDplots.pdf",
       n1+n2+n3+n4+n5+plot_layout(design = layout)+
         plot_annotation(title = 'NMDplots', tag_levels = "a")  &
         theme(plot.tag = element_text(face = 'bold')),
       width = 8.5, 
       height = 11, 
       units = "in")


##PlottingERTG plots similarly
sets_list<-readRDS("../sets_list.rds")
finaltestresults<-data.frame()
files<-list.files("MANE_DESeq/")
n=1
for (id in files) {
  print(id)
  df<-read_tsv(paste0("MANE_DESeq/", id),show_col_types = FALSE)
  name<-str_split(id,"vs", simplify = T)[1]
  print(name)
  df<-na.omit(df)
  
  ##NewClassification
  df <- df %>%
    mutate(Class1 = case_when(
      ensembl_gene_id %in% sets_list$Horste.ER ~ "ER+",
      ensembl_gene_id %in% sets_list$Ensembl.Mem  ~ "Mem",
      ensembl_gene_id %in% sets_list$Horste.TG ~ "TG+",
      TRUE ~ "Others" ))
  
  df<-df%>%filter(Class1 != "Mem")
  table(df$Class1)
  #plot:
  # Calculate pairwise p-values (e.g., Wilcoxon test)
  p_values <- pairwise.wilcox.test(
    df$log2FoldChange, 
    df$Class1, 
    p.adjust.method = "bonferroni"
  )$p.value
  
  # Format p-value text
  pval_text <- paste("p values (Wilcoxon test):",
                     paste("ER+ vs TG+ =", format(p_values["TG+", "ER+"], scientific = TRUE, digits = 2)),
                     paste("ER+ vs Others =", format(p_values["Others", "ER+"], scientific = TRUE, digits = 2)),
                     paste("TG+ vs Others =", format(p_values["TG+","Others"], scientific = TRUE, digits = 2)),
                     sep = "\n")
  
  
  # Calculate the number of entries in each class
  class_counts <- df %>%
    group_by(Class1) %>%
    dplyr::summarize(count = n(), .groups = 'drop')
  
  
  # Create formatted labels with counts
  labels_with_counts <- paste(class_counts$Class1, 
                              " (", class_counts$count, ")", sep = "")
  
  if (id %in% c("DcapvsNluc_MANEDESeq.tsv", "WcapvsNluc_MANEDESeq.tsv")) {
    x_limits <- c(-0.75, 0.75)   # Wider range for specific files
    annotate_x <- -0.7 # Adjust annotation position accordingly
    }      
  else {x_limits <- c(-3, 3)    # Default range
  annotate_x <- -2.9 # Default annotation position
  control<-"/Control"}     

  if (id == "DcapvsNluc_MANEDESeq.tsv") {  
    name <- "DENV Capsid" 
    control<-"/Nluc"} 
  if (id == "WcapvsNluc_MANEDESeq.tsv") {  
    name <- "WNV Capsid" 
    control<-"/Nluc"} 
  if (id == "ZIKV(MR766)vsControlMANEDEseq.tsv") {  
    x_limits <- c(-1, 1)   # Wider range for specific files
    annotate_x <- -0.9 
    control<-"/Control"} 
  
  # Plotting with the updated dataframe and columns
  assign(paste0("t",n) , df %>%
           ggplot(aes(log2FoldChange, color = Class1)) +
           geom_vline(aes(xintercept = 0), linewidth = 0.02)+
           geom_hline(aes(yintercept = 0.5), linewidth = 0.02)+
           coord_cartesian(xlim = x_limits)+
           stat_ecdf(geom = "line", linewidth = 0.4) +
           scale_color_manual(values = c("ER+" = "#9D3F97", "TG+" = "#39B54A", "Others" = "grey45"), 
                              labels = labels_with_counts) +  # Change legend text
           theme_linedraw() +
           theme(axis.text = element_text(size = 5),axis.title = element_text(size = 7),
                 panel.grid.major = element_blank(), 
                 panel.grid.minor = element_blank(),
                 aspect.ratio = 1,
                 legend.position = c(0.98, 0.02),         # Move legend to bottom right
                 legend.justification = c(1, 0),     # Ensure correct alignment at bottom right
                 legend.title=element_text(size=6), 
                 legend.text=element_text(size=5),
                 legend.key=element_blank(),
                 legend.box.background = element_blank(),
                 legend.background = element_blank(),
                 legend.key.size =unit(2, 'mm'),
                 plot.background = element_blank(),
                 plot.margin = unit(c(0, 0, 0, 0), "cm")) +
           labs(y = "Cumulative frequency", x = paste0("log2FC ",name, control),#title = paste0( name," vs Control"),
                color = "Compartment") +
           annotate("text", x = annotate_x, y = 0.85, label = pval_text, hjust = 0, size = 1.5, color = "black") )  # Add p-values on top left
  
  # Get all pairwise combinations
  group_levels <- unique(df[["Class1"]])
  pairs <- combn(group_levels, 2, simplify = FALSE)
  
  # Perform tests and store results
  results <- lapply(pairs, function(pair) {
    group1 <- pair[1]
    group2 <- pair[2]
    
    x <- df[["log2FoldChange"]][df[["Class1"]] == group1]
    y <- df[["log2FoldChange"]][df[["Class1"]] == group2]
    
    test <- wilcox.test(x, y, exact = FALSE)
    
    data.frame(Group = str_split(id, "\\.", simplify = T)[1],
               Comparison = paste0(group1," vs ", group2),
               n_set1 = length(x),
               n_set2 = length(y),
               test = "Wilcoxon rank sum test with continuity correction",
               p_value = test$p.value,
               W_statistic = test$statistic,
               Alternate_Hypothesis = "true location shift is not equal to 0",
               stringsAsFactors = FALSE
    )
  })
  
  results_df <- bind_rows(results)
  
  # Adjust p-values (Bonferroni, same as in your test)
  results_df <- results_df %>%
    mutate(p_adjusted_bonferroni = p.adjust(p_value, method = "bonferroni"))
  finaltestresults<-bind_rows(finaltestresults, results_df)
           
  n=n+1}

print(finaltestresults)
ggsave("ERTGplots.pdf",
       t1+t2+t3+t4+t5+plot_layout(design = layout)+
         plot_annotation(title = 'ERTGplots', tag_levels = "a")  &
         theme(plot.tag = element_text(face = 'bold')),
       width = 8.5, 
       height = 11, 
       units = "in")


##Plottinng final figures

layout <- c(
  area(1,3,2,5),area(1, 6, 2,8),
  area(3,1,4,2),area(3,3,4,4),area(3,5,4,6),area(3,7,4,8),
  area(5,1,8,3),area(5,4,9,7))
plot(layout)
layout <- c(
  area(1,5,2,6),area(1,7, 2,8),
  area(3,1,4,2),area(3,3,4,4),area(3,5,4,6),area(3,7,4,8),
  area(5,1,8,3),area(5,4,9,7))
plot(layout)
ggsave("Figure7patch.pdf",
       t3+n3+t2+n2+t5+n5+plot_layout(design = layout)+
         plot_annotation(title = 'Figure7', tag_levels = "a")  &
         theme(plot.tag = element_text(face = 'bold')),
       width = 8.5, 
       height = 11, 
       units = "in")

layout <- c(area(1,3, 2,4),area(1,5, 2,6),
            area(3,1,4,2),area(3,3,4,4),
  area(7,6,8,7),
  area(5,1,8,3),area(5,4,9,7))
plot(layout)

ggsave("FigureS8patch.pdf",
       t1+n1+t4+n4+plot_layout(design = layout)+
         plot_annotation(title = 'FigureS8', tag_levels = "a")  &
         theme(plot.tag = element_text(face = 'bold')),
       width = 8.5, 
       height = 11, 
       units = "in")



layout <- c(
  area(1,5,2,6),area(1,7, 2,8),
  area(3,1,4,2),area(3,3,4,4),area(3,5,4,6),area(3,7,4,8))
plot(layout)
ggsave("Figure7patch2.pdf",
       t3+n3+t2+n2+t5+n5+plot_layout(design = layout)+
         plot_annotation(tag_levels = "a")  &
         theme(plot.tag = element_text(face = 'bold')),
       width = 180, 
       height = 90, 
       units = "mm")
