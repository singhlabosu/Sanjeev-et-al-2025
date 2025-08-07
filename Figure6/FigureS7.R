setwd("~/OneDrive - The Ohio State University/Singhlab/Bioinfo/PYMpaperFiguresUpdate/Figure6/")
library(tidyverse)
#plotting signalP/transmembrane/mitocarta/pexdb######
PYMkdDESeqlong<-read_csv("wild_vs_kd_CHECKED_res.csv")
names(PYMkdDESeqlong)[1]<-"ensembl_gene_id"
PYMkdDESeqlong<-na.omit(PYMkdDESeqlong)

#Creating a new df
Membrane<-read_tsv("Membrane.tsv")
Membrane<-Membrane%>%filter(chromosome_name %in% c(1:22, "X","Y","MT"))
table(Membrane$chromosome_name)
PYMkdDESeqlong<-inner_join(PYMkdDESeqlong[,1:3],Membrane[,1:3])
names(PYMkdDESeqlong)
#Creating new columns in YES/No format
PYMkdDESeqlong<-PYMkdDESeqlong%>%mutate(
  SignalPeptide = case_when(
    tmhmm=="TMhelix"  ~ "Yes",
    TRUE ~"No"))
PYMkdDESeqlong<-PYMkdDESeqlong%>%mutate(
  TransMembraneDomain = case_when(
    signalp %in% c("SignalP-noTM","SignalP-TM") ~ "Yes",
    TRUE ~"No"))

PYMkdDESeqlong<-PYMkdDESeqlong[,c(1,3,6,7)]


# Read the FASTA file
pexdb_sequences <- readAAStringSet("PexDB_Homo_sapiens.fas")

# Extract headers from the FASTA file
names(pexdb_sequences)
length(pexdb_sequences) # only 100

#Checking ensembl GO annotation
# genes_peroxisome <- getBM(
#   attributes = c("external_gene_name","ensembl_gene_id", "go_id", "name_1006"), # Gene name, GO ID, GO term name, 
#   filters = "go", 
#   values = "GO:0005777", # GO term for peroxisome
#   mart = ensembl100)
# write_tsv(genes_peroxisome,"genes_peroxisome.tsv")
genes_peroxisome<-read_tsv("genes_peroxisome.tsv")
unique(genes_peroxisome$external_gene_name)#121 gene names
unique(genes_peroxisome$ensembl_gene_id)#129 gene ids


#Reading MitoCarta#####
MitoCarta<- read_excel("Human.MitoCarta3.0.xls", 
                       sheet = "A Human MitoCarta3.0")
names(MitoCarta)
MitoCarta<-MitoCarta[,c(1,3,4,15)]
MitoCarta$ensembl_gene_id<-str_split(MitoCarta$EnsemblGeneID_mapping_version_20200130, "\\|",simplify = T)[,1]


PYMkdDESeqlong<-PYMkdDESeqlong%>%
  mutate(MitocartaAnnotation=
           case_when(ensembl_gene_id %in% MitoCarta$ensembl_gene_id~ "Yes",
                     TRUE ~ "No"))
table(PYMkdDESeqlong$MitocartaAnnotation)

PYMkdDESeqlong<-PYMkdDESeqlong%>%
  mutate(PeroxisomeAnnotation=
           case_when(ensembl_gene_id %in% genes_peroxisome$ensembl_gene_id ~ "Yes",
                     TRUE ~ "No"))
table(PYMkdDESeqlong$PeroxisomeAnnotation)

sets_list <- readRDS("sets_list.rds")
PYMkdDESeqlong<-PYMkdDESeqlong%>%
  mutate(MembraneTranslation=
           case_when(ensembl_gene_id %in% sets_list$Jan.ER ~ "Yes",
                     TRUE ~ "No"))


#Converting to long format
PYMkdDESeqlong <- PYMkdDESeqlong %>% 
  pivot_longer(cols = !c(ensembl_gene_id, log2FoldChange), # which columns do we want to "pivot"
               names_to = "Classification", # where should the column names go
               values_to = "Annotation")

PYMkdDESeqlong$Classification <- factor(PYMkdDESeqlong$Classification,
                                        levels = c("MembraneTranslation",
                                                   "TransMembraneDomain",
                                                   "SignalPeptide",
                                                   "MitocartaAnnotation", 
                                                   "PeroxisomeAnnotation"))

table(PYMkdDESeqlong$Classification)
library(ggpubr)
# Calculate medians for each PTC.Status group
median_values <- PYMkdDESeqlong %>%
  group_by(Annotation, Classification) %>%
  dplyr::summarize(median_log2FC = median(log2FoldChange),
                   count = n(), .groups = 'drop')

c<-PYMkdDESeqlong %>%
  ggplot(aes(x = Annotation, y = log2FoldChange, fill = Annotation)) +
  geom_hline(yintercept = 0, size = 0.05) +
  #geom_point(aes(color = Group), position = position_jitterdodge(), alpha = 0.5, size=0.5) +
  #geom_boxplot(position = position_dodge(), varwidth = T, alpha = 0.1, outlier.size = 0.8) +
  geom_boxplot(position = position_dodge(), varwidth = T,outlier.shape = NA)+
  scale_fill_manual(values = c("No" = "grey45", "Yes" = "red3")) +
  scale_color_manual(values = c("No" = "grey45", "Yes" = "red3")) +
  theme_linedraw() +
  coord_cartesian(ylim = c(-1.5,1.5))+
  theme( panel.grid.major = element_blank(),
         panel.grid.minor = element_blank(),
         aspect.ratio = 3,
         legend.position = "none",
         strip.background = element_blank(),
         strip.text = element_text(color= "black",size =6, angle = 0, hjust =0)) +
  labs( y = "log2FC PYM1kd/Control",
        x = "Annotation" ) +
  facet_wrap(~Classification, nrow = 1) +
  stat_compare_means(comparisons = list(c("Yes", "No")),
                     method = "wilcox.test",
                     label = "p.format",
                     label.y = 0.8,size = 3,tip.length = 0.003) + # Adjust as needed for positioning
          # Add medians below each boxplot
  geom_text(data = median_values,
          aes(x = Annotation, y = -1.1, label = sprintf("%.2f", median_log2FC), color = Annotation),
          size = 3, vjust = 1.5, inherit.aes = FALSE) +
  geom_text(data = median_values,
            aes(x = Annotation, y = -1.3, label = sprintf("%s", count)),
            size = 3, vjust = 1.5, inherit.aes = FALSE, color = "black")
c+labs(title = "PYM1kd vs Control foldchanges",
       subtitle = "Grouped by different membrane annotations")

##Plotting exon number cdf#####
#read PYM1kd
library(readr)
PYMkdvsControlDESeq <- read.csv("~/Library/CloudStorage/OneDrive-TheOhioStateUniversity/Singhlab/Bioinfo/PYM/PYMKDvsControlTPMDESeq.txt", sep="")
PYMkdvsControlDESeq<-na.omit(PYMkdvsControlDESeq)
PYMkdvsControlDESeq$ensembl_transcript_id<-str_split(rownames(PYMkdvsControlDESeq), "\\.",simplify = T)[,1]
#PYMkdvsControlDESeq<-PYMkdvsControlDESeq%>%filter(baseMean > 10)
#get isoforminfo
# library(biomaRt)
# listMarts(host='https://apr2020.archive.ensembl.org')
# ensembl100=useMart(host='https://apr2020.archive.ensembl.org', biomart='ENSEMBL_MART_ENSEMBL', dataset='hsapiens_gene_ensembl')
# attributes <- listAttributes(ensembl100)
# filters <- listFilters(ensembl100)
# 
# isoforminfo<-getBM(attributes = c("ensembl_gene_id","ensembl_transcript_id",
#                                   "ensembl_exon_id", "rank","transcript_biotype"),
#                    filters = "ensembl_transcript_id",
#                    values = PYMkdvsControlDESeq$ensembl_transcript_id,
#                    mart = ensembl100,
#                    useCache = FALSE)
# table(isoforminfo$transcript_biotype)
# 
# # #filter for rank
# isoforminfo<-isoforminfo%>%
#   arrange(ensembl_transcript_id,desc(rank))%>%
#   group_by(ensembl_transcript_id)%>%dplyr::slice(1)
# write_tsv(isoforminfo, "isoforminfo.tsv")
isoforminfo<-read_tsv("isoforminfo.tsv")
names(isoforminfo)
PYMkdvsControlDESeq<-inner_join(PYMkdvsControlDESeq,isoforminfo)
summary(PYMkdvsControlDESeq$rank)

PYMkdvsControlDESeq<-PYMkdvsControlDESeq%>%filter(transcript_biotype=="protein_coding")
PYMkdvsControlDESeq<-PYMkdvsControlDESeq%>%filter(!ensembl_gene_id%in%sets_list$Ensembl.Mem)
#PYMkdvsControlDESeq<-PYMkdvsControlDESeq%>%filter(baseMean>20)


#Calculate the third quartile (75th percentile) of the "rank" column
n <- quantile(PYMkdvsControlDESeq$rank, 0.75, na.rm = TRUE)
m <- quantile(PYMkdvsControlDESeq$rank, 0.25, na.rm = TRUE)

# Add a new column with the value of n included in the string
PYMkdvsControlDESeq <- PYMkdvsControlDESeq %>%
  mutate(exon_classification = case_when(rank > n ~ paste0(n, "+", " exons"), 
                                         rank < m ~ paste0(m, "-", " exons"), 
                                         TRUE ~ "others"))


table(PYMkdvsControlDESeq$exon_classification)
PYMkdvsControlDESeq<-PYMkdvsControlDESeq%>%filter(exon_classification != "others")
#plot cdf
#plot:
# Calculate pairwise p-values (e.g., Wilcoxon test)
p_values <- pairwise.wilcox.test(
  PYMkdvsControlDESeq$log2FoldChange, 
  PYMkdvsControlDESeq$exon_classification
)$p.value

# Format p-value text
pval_text <- paste("p-values (Wilcoxon test):",
                   paste("High vs low  =", 
                         format(p_values["5- exons","14+ exons"],
                                scientific = TRUE, digits = 2)),
                   sep = "\n")


# Calculate the number of entries in each class
library(dplyr)
class_counts <- PYMkdvsControlDESeq %>%
  group_by(exon_classification) %>%
  dplyr::summarize(count = dplyr::n(), .groups = 'drop')


# Create formatted labels with counts
labels_with_counts <- paste(class_counts$exon_classification, 
                            " (", class_counts$count, ")", sep = "")
labels_with_counts[1]<-paste0("High, ",labels_with_counts[1])
labels_with_counts[2]<-paste0("Low, ",labels_with_counts[2])

# Plotting with the updated dataframe and columns
d<-PYMkdvsControlDESeq %>%
  ggplot(aes(log2FoldChange, color = exon_classification)) +
  geom_vline(aes(xintercept = 0), linewidth = 0.02)+
  geom_hline(aes(yintercept = 0.5), linewidth = 0.02)+
  stat_ecdf(geom = "line", linewidth = 0.8) +
  scale_color_manual(values = c("14+ exons" = "#127CA3", "5- exons" = "#43B285"), 
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
  labs(y = "Cumulative frequency", x = "log2FC PYM1kd/Control",
       color = "Exon count") +
  annotate("text", x = -1.9, y = 0.92, label = pval_text, 
           hjust = 0, size = 4, color = "black")  # Add p-values on top left

d+labs(title = "PYM1 kd  vs Control",
       subtitle = "Grouped based on Exon number")
res <- wilcox.test(log2FoldChange ~ exon_classification, data = PYMkdvsControlDESeq,
                   exact = FALSE)
res
#plotting Single vs multi exon######
PYMkdDESeq<-read_csv("wild_vs_kd_CHECKED_res.csv")
names(PYMkdDESeq)[1]<-"ensembl_gene_id"
PYMkdDESeq<-na.omit(PYMkdDESeq)
#MANE isoforms
MANEisoforms<-read_tsv("../MANEexons.tsv")
names(MANEisoforms)
PYMkdDESeq<-inner_join(PYMkdDESeq,MANEisoforms)
summary(PYMkdDESeq$rank)


# Add a new column with the value of n included in the string
PYMkdDESeq <- PYMkdDESeq %>%
  mutate(Exon_count = case_when(rank > 1 ~ "Multi", 
                                rank == 1 ~ " Single"))


table(PYMkdDESeq$Exon_count)
#plot cdf
#plot:
# Calculate pairwise p-values (e.g., Wilcoxon test)
p_values <- pairwise.wilcox.test(
  PYMkdDESeq$log2FoldChange, 
  PYMkdDESeq$Exon_count
)$p.value

# Format p-value text
pval_text <- paste("p-value (Wilcoxon test)",
                   paste("Single exon vs Multi exon:\n", 
                         format(p_values[1],
                                scientific = TRUE, digits = 2)),
                   sep = "\n")

# Calculate the number of entries in each class
class_counts <- PYMkdDESeq %>%
  group_by(Exon_count) %>%
  dplyr::summarize(count = dplyr::n(), .groups = 'drop')


# Create formatted labels with counts
labels_with_counts <- paste(class_counts$Exon_count, 
                            " (", class_counts$count, ")", sep = "")

# Plotting with the updated dataframe and columns
b<-PYMkdDESeq %>%
  ggplot(aes(log2FoldChange, color = Exon_count)) +
  geom_vline(aes(xintercept = 0), linewidth = 0.02)+
  geom_hline(aes(yintercept = 0.5), linewidth = 0.02)+
  stat_ecdf(geom = "line", linewidth = 0.8) +
  scale_color_manual(values = c("#37469D", "#747474"), 
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
  labs(y = "Cumulative frequency", x = "log2FC PYM1kd/Control",
       color = "Exon count") +
  annotate("text", x = -1.9, y = 0.92, label = pval_text, 
           hjust = 0, size = 4, color = "black")  # Add p-values on top left

b+labs(title = "PYM1 kd vs Control",
       subtitle = "Grouped based on Exon number")
res <- wilcox.test(log2FoldChange ~ Exon_count, data = PYMkdDESeq,
                   exact = FALSE)
res
#Comparing with SMG6/7dkd####

SMG67DESeq<-read_tsv("SMG67_MANE_DESeqResult.tsv")
SMG67DESeq<-na.omit(SMG67DESeq)
sets_list <- readRDS("sets_list.rds")

SMG67DESeq<-SMG67DESeq%>%mutate(Class1 = case_when(
  gene_id %in% sets_list$Horste.ER ~ "ER+",
  gene_id %in% sets_list$Ensembl.Mem  ~ "Mem",
  gene_id %in% sets_list$Horste.TG ~ "TG+",
  TRUE ~ "Others" ))
SMG67DESeq<-SMG67DESeq%>%filter(Class1 != "Mem")
table(SMG67DESeq$Class1)

#names(SMG67DESeq)[7]<-"ensembl_gene_id"
# Calculate pairwise p-values (e.g., Wilcoxon test)
p_values <- pairwise.wilcox.test(
  SMG67DESeq$log2FoldChange, 
  SMG67DESeq$Class1, 
  p.adjust.method = "bonferroni"
)$p.value

# Format p-value text
pval_text <- paste("p-values (Wilcoxon test):",
                   paste("ER+ vs TG+ =", format(p_values["TG+", "ER+"], scientific = TRUE, digits = 2)),
                   paste("ER+ vs Others =", format(p_values["Others", "ER+"], scientific = TRUE, digits = 2)),
                   paste("TG+ vs Others =", format(p_values["TG+","Others"], scientific = TRUE, digits = 2)),
                   sep = "\n")


# Calculate the number of entries in each class
library(tidyverse)
class_counts <- SMG67DESeq %>%
  dplyr::group_by(Class1) %>%
  dplyr::summarize(count = dplyr::n(), .groups = 'drop')




# Create formatted labels with counts
labels_with_counts <- paste(class_counts$Class1, 
                            " (", class_counts$count, ")", sep = "")

# Plotting with the updated dataframe and columns
e<-SMG67DESeq %>%
  ggplot(aes(log2FoldChange, color = Class1)) +
  geom_vline(aes(xintercept = 0), linewidth = 0.02)+
  geom_hline(aes(yintercept = 0.5), linewidth = 0.02)+
  stat_ecdf(geom = "line", linewidth = 0.8) +
  scale_color_manual(values = c("ER+" = "#9D3F97", "TG+" = "#39B54A", "Others" = "grey60"), 
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
  labs(y = "Cumulative frequency", x = "log2FC SMG7KO+SMG6kd/Control",
       color = "Compartment") +
  annotate("text", x = -1.9, y = 0.92, label = pval_text, hjust = 0, size = 4, color = "black")  # Add p-values on top left

e+labs(title = "SMG6/7kd fold changes",
       subtitle = "Grouped by Compartments")


#Comparing with thapsigargin treatment####

DESeqFritz<-read_tsv("DESEQFritz.tsv")

names(DESeqFritz)[7]<-"ensembl_gene_id"

DESeqFritz<-DESeqFritz%>%mutate(Class1 = case_when(
  ensembl_gene_id %in% sets_list$Horste.ER ~ "ER+",
  ensembl_gene_id %in% sets_list$Ensembl.Mem  ~ "Mem",
  ensembl_gene_id %in% sets_list$Horste.TG ~ "TG+",
  TRUE ~ "Others" ))
DESeqFritz<-DESeqFritz%>%filter(Class1 != "Mem")


# Calculate pairwise p-values (e.g., Wilcoxon test)
p_values <- pairwise.wilcox.test(
  DESeqFritz$log2FoldChange, 
  DESeqFritz$Class1, 
  p.adjust.method = "bonferroni"
)$p.value

# Format p-value text
pval_text <- paste("p-values (Wilcoxon test):",
                   paste("ER+ vs TG+ =", format(p_values["TG+", "ER+"], scientific = TRUE, digits = 2)),
                   paste("ER+ vs Others =", format(p_values["Others", "ER+"], scientific = TRUE, digits = 2)),
                   paste("TG+ vs Others =", format(p_values["TG+","Others"], scientific = TRUE, digits = 2)),
                   sep = "\n")


# Calculate the number of entries in each class
library(dplyr)
class_counts <- DESeqFritz %>%
  group_by(Class1) %>%
  dplyr::summarize(count = n(), .groups = 'drop')



# Create formatted labels with counts
labels_with_counts <- paste(class_counts$Class1, 
                            " (", class_counts$count, ")", sep = "")

# Plotting with the updated dataframe and columns
h<-DESeqFritz %>%
  ggplot(aes(log2FoldChange, color = Class1)) +
  geom_vline(aes(xintercept = 0), linewidth = 0.02)+
  geom_hline(aes(yintercept = 0.5), linewidth = 0.02)+
  coord_cartesian(xlim = c(-2,2))+
  stat_ecdf(geom = "line", linewidth = 0.8) +
  scale_color_manual(values = c("ER+" = "#9D3F97", "TG+" = "#39B54A", "Others" = "grey60"), 
                     labels = labels_with_counts) +  # Change legend text
  theme_linedraw() +
  theme(
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    aspect.ratio = 1,
    legend.position = c(0.95, 0.05),         # Move legend to bottom right
    legend.justification = c(1, 0)     # Ensure correct alignment at bottom right
  ) +
  labs(y = "Cumulative frequency", x = "Log2FC Thapsigargin/DMSO",
       color = "Compartment") +
  annotate("text", x = -1.9, y = 0.92, label = pval_text, hjust = 0, size = 4, color = "black")  # Add p-values on top left
h+labs(title = "Thapsigargin treatment fold changes",
       subtitle = "Grouped by Compartments")

#CDF: UPF1LL upgenes in PYMKD#####
UPF1LLkdDEseq <- read_delim("UPF1LLkdDEseq.txt", 
                            delim = "\t", escape_double = FALSE, 
                            trim_ws = TRUE)

UPF1LLkdDEseq<-UPF1LLkdDEseq%>%filter(`siUPF1LL vs NT log2(fold change)`>0.5&`False discovery rate (FDR)`<0.05)

PYMkdDESeq<-read_csv("wild_vs_kd_CHECKED_res.csv")
names(PYMkdDESeq)[1]<-"ensembl_gene_id"
PYMkdDESeq<-na.omit(PYMkdDESeq)
PYMkdDESeq<-PYMkdDESeq%>%mutate(UPF1LLkdStatus=
                                  case_when(ensembl_gene_id %in% UPF1LLkdDEseq$`ENSEMBL Gene ID` ~ "Upregulated",
                                            TRUE ~ "Others"))
table(PYMkdDESeq$UPF1LLkdStatus)
# Calculate pairwise p-values (e.g., Wilcoxon test)
p_values <- wilcox.test(log2FoldChange ~ UPF1LLkdStatus, data = PYMkdDESeq,
                        exact = FALSE)$p.value

# Format p-value text
pval_text <- paste("p-values (Wilcoxon test):","UPF1LL targets",
                   paste( "vs Others =", format(p_values, scientific = TRUE, digits = 2)),
                   sep = "\n")


# Calculate the number of entries in each class
class_counts <- PYMkdDESeq %>%
  dplyr::group_by(UPF1LLkdStatus) %>%
  dplyr::summarize(count = dplyr::n(), .groups = 'drop')

# Create formatted labels with counts
labels_with_counts <- paste(class_counts$UPF1LLkdStatus, 
                            " (", class_counts$count, ")", sep = "")

g<-PYMkdDESeq%>%ggplot(aes(log2FoldChange, colour = UPF1LLkdStatus)) +
  geom_vline(aes(xintercept = 0), linewidth = 0.02)+
  geom_hline(aes(yintercept = 0.5), linewidth = 0.02)+
  stat_ecdf(geom = "step",size=1.0)+
  coord_cartesian(xlim = c(-2, 2)) +
  theme_linedraw() +
  scale_color_discrete(labels = labels_with_counts,
                     name = "UPF1LL kd status\n(Fritz et al., (2022))")+
  theme(
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    aspect.ratio = 1,
    legend.background = element_blank(),        # Removes the white box background
    legend.key = element_blank() ,
    legend.position = c(0.95, 0.05),         # Move legend to bottom right
    legend.justification = c(1, 0)     # Ensure correct alignment at bottom right
  ) +
  labs(x = "log2 FC PYM1kd/Control", y = "Cumulative frequency")+
  annotate("text", x = -1.9, y = 0.92, label = pval_text, hjust = 0, size = 4, color = "black") 

g+labs(title = "CDF plot, log2FC PYM1kd vs Control",
       subtitle = "Grouped based on UPF1 LL kd status, Fritz et al 2022")

res <- wilcox.test(log2FoldChange ~ UPF1LLkdStatus, data = PYMkdDESeq,
                   exact = FALSE)
res




#CDF: NBAS upgenes in PYMKD#####

NBASkdDEseq<-read_delim("LongmanTableS1.txt", 
                        delim = "\t", escape_double = FALSE, 
                        trim_ws = TRUE)
NBASkdDEseq<-NBASkdDEseq%>%filter(NBAS_padj < 0.05 & NBAS_log2FC > 0)

PYMkdDESeq<-PYMkdDESeq%>%mutate(NBASkdStatus=
                                case_when(ensembl_gene_id %in% NBASkdDEseq$gene_id ~ "Upregulated",
                                            TRUE ~ "Others"))

table(PYMkdDESeq$NBASkdStatus)
# Calculate pairwise p-values (e.g., Wilcoxon test)
p_values <- wilcox.test(log2FoldChange ~ NBASkdStatus, data = PYMkdDESeq,
                        exact = FALSE)$p.value

# Format p-value text
pval_text <- paste("p-values (Wilcoxon test):","NBAS targets",
                   paste( "vs Others =", format(p_values, scientific = TRUE, digits = 2)),
                   sep = "\n")


# Calculate the number of entries in each class
class_counts <- PYMkdDESeq %>%
  dplyr::group_by(NBASkdStatus) %>%
  dplyr::summarize(count = dplyr::n(), .groups = 'drop')

# Create formatted labels with counts
labels_with_counts <- paste(class_counts$NBASkdStatus, 
                            " (", class_counts$count, ")", sep = "")

f<-PYMkdDESeq%>%ggplot(aes(log2FoldChange, colour = NBASkdStatus)) +
  geom_vline(aes(xintercept = 0), linewidth = 0.02)+
  geom_hline(aes(yintercept = 0.5), linewidth = 0.02)+
  stat_ecdf(geom = "step",size=1.0)+
  coord_cartesian(xlim = c(-2, 2)) +
  theme_linedraw() +
  scale_color_discrete(labels = labels_with_counts,
                       name = "NBAS kd status\n(Longman et al., (2020))")+
  theme(
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    aspect.ratio = 1,
    legend.background = element_blank(),        # Removes the white box background
    legend.key = element_blank() ,
    legend.position = c(0.95, 0.05),         # Move legend to bottom right
    legend.justification = c(1, 0)     # Ensure correct alignment at bottom right
  ) +
  labs(x = "log2 FC PYM1kd/Control", y = "Cumulative frequency")+
  annotate("text", x = -1.9, y = 0.92, label = pval_text, hjust = 0, size = 4, color = "black") 

f+labs(title = "CDF plot, L2FC PYMkd vs Control",
       subtitle = "Grouped based on NBAS kd status, Longman et al 2020")

res <- wilcox.test(log2FoldChange ~ NBASkdStatus, data = PYMkdDESeq,
                   exact = FALSE, alternative = "greater")
res





#Making overlap Venn#####

###makinglists
sets_list2 <- list( NBAS_Up= na.omit(NBASkdDEseq$gene_id[NBASkdDEseq$NBAS_padj < 0.05 & NBASkdDEseq$NBAS_log2FC > 0]),
                   UPF1LL_Up= na.omit(UPF1LLkdDEseq$`ENSEMBL Gene ID`[UPF1LLkdDEseq$`siUPF1LL vs NT log2(fold change)`>0.5&UPF1LLkdDEseq$`False discovery rate (FDR)`<0.05]),
                   PYM_Down=na.omit(PYMkdDESeq$ensembl_gene_id[PYMkdDESeq$log2FoldChange>0&PYMkdDESeq$padj<0.05]),
                   PYM_Up=na.omit(PYMkdDESeq$ensembl_gene_id[PYMkdDESeq$log2FoldChange<0&PYMkdDESeq$padj<0.05]))

# Fit the euler diagram
fit <- euler(sets_list2, shape= "ellipse")
#fit <- euler(sets_list)
# Plot
i<-plot(fit, quantities = T)
i

#Patchwork######
library(patchwork)
layout <- c(area(1, 9, 3,11),
            area(4,1,6,5),area(4, 6, 6,8),area(4, 9, 6,11),
 area(7,4,14,6))
plot(layout)




c+d+e+f+g+h+i+plot_layout(design = layout)+
  plot_annotation(title = 'FIGURE S7', tag_levels = "a")  &
  theme(plot.tag = element_text(face = 'bold')) 
ggsave("S7patch2.pdf",
       b+c+d+e+plot_layout(design = layout)+
         plot_annotation(title = 'FIGURE 6', tag_levels = "a")  &
         theme(plot.tag = element_text(face = 'bold')),
       width = 17, 
       height = 22, 
       units = "in")

layout <- c(
  area(4,1,6,3),area(4,4,6,6),area(4,7,6,9),area(4,10,6,11),
  area(5,4,8,6))
ggsave("S7patch3.pdf",
       f+g+h+i+plot_layout(design = layout)+
         plot_annotation(title = 'FIGURE 6', tag_levels = "a")  &
         theme(plot.tag = element_text(face = 'bold')),
       width = 17, 
       height = 22, 
       units = "in")


