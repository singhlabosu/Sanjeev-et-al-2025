#Figure6

setwd("~/OneDrive - The Ohio State University/Singhlab/Bioinfo/PYMpaperFiguresUpdate/Figure6/")

#Comparing overlaps of gene classifications
library(tidyverse)
##Reading Jan et. al table & Classification#####
ERdata<-read_tsv("NIHMS646671-supplement-Table_S6.txt",col_names = T,skip = 23)
names(ERdata)
ERdata<-ERdata[,c(1:3,14:18)]

ERdata%>%ggplot(aes(log2.enrichment))+geom_density()+geom_vline(xintercept = 0.58, size = 0.5, linetype = 2)

ERdata<-ERdata %>% filter(!is.na(log2.enrichment))
ERdata<-ERdata%>%mutate(
  Class=case_when(log2.enrichment > 0.58 ~ "ER enriched",
                  TRUE ~ "Others"))
table(ERdata$Class)
glimpse(ERdata)

##Getting geneids for ucsc id
UCSCtable<-read_tsv("UCSCtoEnsembl.tsv")
names(UCSCtable)
UCSCtable<-UCSCtable[,c(1,12)]
UCSCtable<-UCSCtable%>%mutate(alignID=str_split(alignID, "\\.", simplify = T)[,1])
ERdata<-ERdata%>%mutate(alignID=str_split(UCSC.gene.id, "\\.", simplify = T)[,1])
dim(left_join(ERdata,UCSCtable))
dim(inner_join(ERdata,UCSCtable))
4804-4478
#326 rows with missing ensembl transcript id
ERdata<-inner_join(ERdata,UCSCtable)




###Retreiving transcript information from biomart
# library(biomaRt)
# listMarts(host='https://apr2020.archive.ensembl.org')
# ensembl100=useMart(host='https://apr2020.archive.ensembl.org', biomart='ENSEMBL_MART_ENSEMBL', dataset='hsapiens_gene_ensembl')
# attributes <- listAttributes(ensembl100)
# filters <- listFilters(ensembl100)
# 
# names(ERdata)[2]<-"uniprot_gn_id"
# names(ERdata)[11]<-"ensembl_transcript_id"
# ERdata<-ERdata%>%mutate(ensembl_transcript_id=str_split(ensembl_transcript_id, "\\.", simplify = T)[,1])
# 
# Geneids <- getBM(attributes = c("ensembl_gene_id",
#                                 "ensembl_transcript_id"),
#                  filters = "ensembl_transcript_id",
#                  values = ERdata$ensembl_transcript_id,
#                  mart = ensembl100,
#                  useCache = FALSE)
# write_tsv(Geneids,"Geneids.tsv")
Geneids<-read_tsv("Geneids.tsv")
names(ERdata)[11]<-"ensembl_transcript_id"
ERdata$ensembl_transcript_id<-str_split(ERdata$ensembl_transcript_id, "\\.", simplify = T)[,1]
ERdata<-inner_join(ERdata,Geneids)

#Reading Horste et. al table#####
library(readxl)
Subcell <- read_excel("1-s2.0-S109727652300970X-mmc2.xlsx", na = "")

Subcell[Subcell == "."] <- NA # Replace all dots with NA
Subcell<-Subcell%>%mutate(Compartment=case_when(
  `TMD classification` %in% c(11,1) ~ "TG+",
  `TMD classification` %in% c(31,3) ~ "ER+",
  `TMD classification` %in% c(41,4) ~ "CY+",
  `TMD classification` %in% c(51,5) ~ "Unbiased"))
table(Subcell$Compartment)

1476*100/(1476+226+212+226)

names(Subcell)
Subcell<-Subcell[,c(1:5,33)]
Subcell<-Subcell%>%distinct()
# 
# GeneNames <- getBM(attributes = c("ensembl_gene_id",
#                                 "external_gene_name","chromosome_name"),
#                  filters = "external_gene_name",
#                  values = Subcell$Gene_name,
#                  mart = ensembl100,
#                  useCache = FALSE)
# write_tsv(GeneNames,"GeneNames.tsv")
GeneNames<-read_tsv("GeneNames.tsv")
GeneNames %>% group_by(external_gene_name) %>% filter(n() > 1)
#1448 duplicated rows
9155+1448

GeneNames<-GeneNames%>%filter(chromosome_name %in% c(1:22, "X", "Y", "MT"))
GeneNames %>% group_by(external_gene_name) %>% filter(n() > 1)
#20 duplicated rows->10 extra rows

names(Subcell)[1]<-"external_gene_name"
Subcell<-left_join(Subcell, GeneNames)

#Membrane ann. Ensembl####

#Getting membrane annotation from ensembl
# Membrane<-getBM(attributes = c("ensembl_gene_id","signalp","tmhmm",
#                                "gene_biotype","chromosome_name"),
#                 filters = "biotype",
#                 values = "protein_coding",
#                 mart = ensembl100,
#                 useCache = FALSE)#28463
# write_tsv(Membrane,"Membrane.tsv")
Membrane<-read_tsv("Membrane.tsv")
Membrane<-Membrane%>%filter(chromosome_name %in% c(1:22, "X","Y","MT"))
table(Membrane$chromosome_name)
Membrane<-Membrane%>%mutate(Annotation=case_when(tmhmm == "TMhelix" ~ "Membrane",
                                                 signalp == "SignalP-noTM" ~ "Membrane",
                                                 signalp == "SignalP-TM" ~ "Membrane",
                                                 TRUE ~ "Other"))
table(Membrane$Annotation)
table(Membrane$tmhmm)
table(Membrane$signalp)

#############Overlap between Classifications############

###makinglists

sets_list <- list( Jan.ER= unique(ERdata$ensembl_gene_id[ERdata$Class == "ER enriched"]),
                   Horste.ER= na.omit(unique(Subcell$ensembl_gene_id[Subcell$Compartment == "ER+"])),
                   Horste.CY= na.omit(unique(Subcell$ensembl_gene_id[Subcell$Compartment == "CY+"])),
                   Horste.TG= na.omit(unique(Subcell$ensembl_gene_id[Subcell$Compartment == "TG+"])),
                   Horste.UB= na.omit(unique(Subcell$ensembl_gene_id[Subcell$Compartment == "Unbiased"])),
                   Ensembl.Mem= unique(Membrane$ensembl_gene_id[Membrane$Annotation == "Membrane"]))
saveRDS(sets_list, file = "sets_list.rds")
sets_list <- readRDS("sets_list.rds")
#sets_list <- na.omit(sets_list)            
intersect(sets_list$Horste.ER,sets_list$Horste.TG)
library(eulerr)
# Fit the euler diagram
fit <- euler(sets_list, shape= "ellipse")

# Plot
plot(fit, quantities = T)
#Reading PexDB#####
library(Biostrings)

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
#plotting ER_TG_cdf####

PYMkdDESeq<-read_csv("wild_vs_kd_CHECKED_res.csv")
names(PYMkdDESeq)[1]<-"ensembl_gene_id"
PYMkdDESeq<-na.omit(PYMkdDESeq)

##NewClassification
PYMkdDESeq <- PYMkdDESeq %>%
  mutate(Class1 = case_when(
    ensembl_gene_id %in% sets_list$Horste.ER ~ "ER+",
    ensembl_gene_id %in% sets_list$Ensembl.Mem  ~ "Mem",
    ensembl_gene_id %in% sets_list$Horste.TG ~ "TG+",
    TRUE ~ "Others" ))
Subcell<-inner_join(Subcell,PYMkdDESeq)

PYMkdDESeq<-PYMkdDESeq%>%filter(Class1 != "Mem")
table(PYMkdDESeq$Class1)
#plot:
# Calculate pairwise p-values (e.g., Wilcoxon test)
p_values <- pairwise.wilcox.test(
  PYMkdDESeq$log2FoldChange, 
  PYMkdDESeq$Class1, 
  p.adjust.method = "bonferroni"
)$p.value

# Format p-value text
pval_text <- paste("p values (Wilcoxon test):",
                   paste("ER+ vs TG+ =", format(p_values["TG+", "ER+"], scientific = TRUE, digits = 2)),
                   paste("ER+ vs Others =", format(p_values["Others", "ER+"], scientific = TRUE, digits = 2)),
                   paste("TG+ vs Others =", format(p_values["TG+","Others"], scientific = TRUE, digits = 2)),
                   sep = "\n")


# Calculate the number of entries in each class
class_counts <- PYMkdDESeq %>%
  group_by(Class1) %>%
  dplyr::summarize(count = n(), .groups = 'drop')


# Create formatted labels with counts
labels_with_counts <- paste(class_counts$Class1, 
                            " (", class_counts$count, ")", sep = "")

# Plotting with the updated dataframe and columns
a<-PYMkdDESeq %>%
  ggplot(aes(log2FoldChange, color = Class1)) +
  geom_vline(aes(xintercept = 0), linewidth = 0.02)+
  geom_hline(aes(yintercept = 0.5), linewidth = 0.02)+
  coord_cartesian(xlim = c(-2,2))+
  stat_ecdf(geom = "line", linewidth = 0.4) +
  scale_color_manual(values = c("ER+" = "#9D3F97", "TG+" = "#39B54A", "Others" = "grey45"), 
                     labels = labels_with_counts) +  # Change legend text
  theme_linedraw() +
  theme(
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    aspect.ratio = 1,
    legend.position = c(0.98, 0.02),         # Move legend to bottom right
    legend.justification = c(1, 0),     # Ensure correct alignment at bottom right
    legend.key.size =unit(2, 'mm'),
    legend.title=element_text(size=6), 
    legend.text=element_text(size=5),
    axis.title=element_text(size=7), 
    axis.text=element_text(size=5),
    legend.key=element_blank(),
    legend.box.background = element_blank(),
    panel.background = element_blank(),
    plot.background = element_blank(),
    plot.margin = unit(c(0, 0, 0, 0), "cm")) +
  labs(y = "Cumulative frequency", x = "log2FC PYM1kd/Control",
       color = "Compartment") +
  annotate("text", x = -2, y = 0.85, label = pval_text, hjust = 0, size = 1.5, color = "black")  # Add p-values on top left
a


#plotting ER_TGqPCR#####
#ER targets first
exp1<-read_csv("~/OneDrive - The Ohio State University/Singhlab/qPCR/20230316/MS_2023-03-16 19-48-37_PYMERexp1.csv",
               skip = 19)
exp2<-read_csv("~/OneDrive - The Ohio State University/Singhlab/qPCR/20230316/MS_2023-03-16 21-39-52_PYMERexp2.csv",
               skip = 19)
exp3<-read_csv("~/OneDrive - The Ohio State University/Singhlab/qPCR/20230316/MS_2023-03-21 12-41-06_PYMERexp3.csv",
               skip = 19)
combinedCqs<-bind_rows(exp1,exp2,exp3)
table(combinedCqs$Target)
combinedCqs <- combinedCqs %>%
  filter(!Target %in% c("-ve", "C1QBP_A", "C5ORF15_A", "RPS9-", "RPS9+", "SRSF3-", "SRSF3+"))
combinedCqs<-combinedCqs%>%
  mutate(Cq=as.numeric(as.character(Cq)))
combinedCqs<-na.omit(combinedCqs)
library(Rmisc)
CqSummary <- summarySE(combinedCqs, measurevar="Cq", groupvars=c("Sample","Target"))
CqSummary<-CqSummary[,1:5]
delcqdf<-tibble()
for(i in unique(CqSummary$Sample)){
  print(i)
  filteredCqsummary<-CqSummary%>%filter(Sample==i)%>%mutate(
    Delcq = Cq - filter(CqSummary,Sample==i&Target=="ACTB")[,4] )
  delcqdf<-rbind(delcqdf,filteredCqsummary)
}
delcqdf<-mutate(delcqdf, Rel.levels = 2^-Delcq)
library(stringr)
delcqdf$Treatment<-str_split(delcqdf$Sample, "\\_", simplify = T)[,1]
delcqdf$Experiment<-str_split(delcqdf$Sample, "\\_", simplify = T)[,2]

finaldf<-tibble()
filteredCqsummary<-tibble()

for(i in unique(delcqdf$Target)){
  print(i)
  for(j in unique(delcqdf$Experiment)){
    print(j)
    filteredCqsummary<-delcqdf%>%filter(Target==i&Experiment==j)
    filteredCqsummary<-filteredCqsummary%>%mutate(
      Rel.level.tarnorm = Rel.levels/filter(filteredCqsummary,Treatment=="siNC"&Target==i)[,7] )
    finaldf<-rbind(finaldf,filteredCqsummary)}}

finaldf%>%filter(Treatment %in% c("siNC", "siPYM"))%>%
  ggplot(aes(Target, Rel.level.tarnorm, 
             fill = factor(Treatment,
                           levels=c("siNC","siPYM","siE","siPE"))))+
  scale_fill_brewer(palette="Dark2")+
  stat_summary(fun = mean, geom = "bar", position = "dodge", width=0.6, alpha =0.7)+
  stat_summary(fun.data = mean_se, 
               geom = "errorbar",width = 0.2, position = position_dodge(width = 0.6))+
  geom_point(shape = 1,size=0.75, position = position_dodge(width = 0.6))+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 60,hjust=1))+
  labs(title="Relative levels : 3 Bioreps", y="Relative Level (Normalized to Actin/siNC)")+
  guides(fill=guide_legend(title="Treatment"))

#TG transcripts
exp1<-read.csv("~/OneDrive - The Ohio State University/Singhlab/qPCR/20240518_PYM_4A3_ERTG/MS_2024-05-17 12-29-59_PYM1_4A3dkd_ERTG.csv",header = T)
exp2<-read.csv("~/OneDrive - The Ohio State University/Singhlab/qPCR/20240518_PYM_4A3_ERTG/MS_2024-05-17 14-59-34_PYM1_4A3dkd_ERTG.csv",header = T)
exp3<-read.csv("~/OneDrive - The Ohio State University/Singhlab/qPCR/20240518_PYM_4A3_ERTG/MS_2024-05-17 16-47-38_PYM1_4A3dkd_ERTG.csv",header = T)
exp4<-read.csv("~/OneDrive - The Ohio State University/Singhlab/qPCR/20240518_PYM_4A3_ERTG/MS_2024-05-17 18-23-21_PYM1_4A3dkd_ERTG.csv",header = T)

qPCR<-rbind(exp1,exp2,exp3,exp4)
table(qPCR$Target)
qPCR <- qPCR %>%
  mutate_at(vars(Target), ~ gsub("Z394", "ZNF394", .))
CqSummary <- summarySE(qPCR, measurevar="Cq", groupvars=c("Sample","Target"))
table(CqSummary$Sample)
CqSummary <- CqSummary[2:117,1:5]%>%filter(Target != "ZFP36L1")
delcqdf<-tibble()
for(i in unique(CqSummary$Sample)){
  print(i)
  filteredCqsummary<-CqSummary%>%filter(Sample==i)%>%mutate(
    Delcq = Cq - filter(CqSummary,Sample==i&Target=="ACTB")[,4] )
  delcqdf<-rbind(delcqdf,filteredCqsummary)
}
delcqdf<-mutate(delcqdf, Rel.levels = 2^-Delcq)

delcqdf$Treatment<-str_split(delcqdf$Sample, "\\_", simplify = T)[,1]
delcqdf$Experiment<-str_split(delcqdf$Sample, "\\_", simplify = T)[,2]

delcqdf<-delcqdf%>%
  filter(Treatment %in% c("siNC","siP"))
filteredCqsummary<-tibble()

for(i in unique(delcqdf$Target)){
  print(i)
  for(j in unique(delcqdf$Experiment)){
    print(j)
    filteredCqsummary<-delcqdf%>%filter(Target==i&Experiment==j)
    filteredCqsummary<-filteredCqsummary%>%mutate(
      Rel.level.tarnorm = Rel.levels/filter(filteredCqsummary,Treatment=="siNC"&Target==i)[,7] )
    finaldf<-rbind(finaldf,filteredCqsummary)}}


names(finaldf)

###modified plot
library(wesanderson)
table(finaldf$Target)
finaldf$Target2<-str_split(finaldf$Target, "\\_",simplify = T)[,1]
finaldf%>%filter(Target2=="SCD")
finaldf<-finaldf%>%filter(Target!="SCD_B")
finaldf%>%filter(Target2=="SCD")
finaldf<-finaldf%>% mutate(Treatment = case_when(Target == "SCD" & Treatment == "siP" ~ "siPYM",
                                                 TRUE ~ Treatment))


b<-finaldf%>%filter(Treatment %in% c("siNC", "siPYM", "siP"))%>%
  ggplot(aes(factor(Target2,
                    levels=c("ACTB","PYM1","TMEM109","SFXN4","REEP5","SCD",
                             "C5ORF15","C1QBP","MSMO1","ZNF394","TIS11B","LIPT1",
                             "ARC","THAP1","SERTAD1","ATF3")), 
             Rel.level.tarnorm, 
             fill = factor(Treatment,
                           levels=c("siNC","siPYM","siP"))))+
  scale_fill_manual(values = c("siPYM" = "#9D3F97", "siP" = "#39B54A", "siNC" = "grey60"))+
  stat_summary(fun = mean, geom = "bar", position = "dodge", width=0.6)+
  stat_summary(fun.data = mean_se, 
               geom = "errorbar",width = 0.2, position = position_dodge(width = 0.6))+
  geom_point(aes(group = Treatment),
             shape = 1, size = 0.3, alpha=0.5,
             position = position_dodge(width = 0.6)) +
  theme_classic()+
  theme(panel.grid.major = element_line(linetype = "blank"), 
        panel.grid.minor = element_line(linetype = "blank"),
        axis.text.x = element_text(size = 5,angle = 60,hjust=1),
        legend.direction = "horizontal", 
        legend.position = "top",
        axis.line = element_line(linewidth  = 0.25),
        aspect.ratio = 0.4,
        legend.key.size =unit(2, 'mm'),
        legend.title=element_text(size=6), 
        legend.text=element_text(size=5),
        axis.title=element_text(size=7), 
        axis.text.y=element_text(size=5),
        legend.key=element_blank(),
        legend.box.background = element_blank(),
        panel.background = element_blank(),
        plot.background = element_blank(),
        plot.margin = unit(c(0, 0, 0, 0), "cm"))+
  labs(y = "Relative Levels", x = "Targets", fill = "Treatment:")
b


# Make sure RepNo is consistent and numeric for correct pairing
finaldf <- finaldf %>%
  mutate(RepNo2 = as.numeric(Experiment),
         log2fc = log(Rel.level.tarnorm, base = 2))  # extract number if needed

# Initialize empty results dataframe
test_results <- data.frame(
  Target = character(),
  Comparison = character(),
  t_statistic = numeric(),
  df = numeric(),
  p_value = numeric(),
  stringsAsFactors = FALSE
)

# List of unique targets
targets <- unique(finaldf$Target2)
targets <- targets[2:16]
# Define comparisons
comparisons <- list(
  c("siNC", "siPYM"),
  c("siNC", "siP")
)

# Loop through each target
for (tgt in targets) {
  print(tgt)
  for (cmp in comparisons) {
    comp1 <- cmp[1]
    comp2 <- cmp[2]
    print(cmp[2])
    # Subset and reshape data for paired comparison
    sub_df <- finaldf %>%
      filter(Target2 == tgt, Treatment %in% c(comp1, comp2)) %>%
      dplyr::select(RepNo2, Treatment, log2fc) %>%
      pivot_wider(names_from = Treatment, values_from = log2fc)
    
    # Skip if both conditions not present or not enough data
    if (!all(c(comp1, comp2) %in% names(sub_df))) next
    sub_df <- drop_na(sub_df)
    if (nrow(sub_df) < 2) next  # need at least 2 paired values
    
    # Perform one-tailed paired t-test: treatment > control (siX > siNC)
    test <- t.test(sub_df[[comp2]], sub_df[[comp1]],
                   paired = TRUE,
                   alternative = "two.sided")
    
    # Append results
    test_results <- rbind(test_results, data.frame(
      Target = tgt,
      Comparison = paste(comp2, "vs", comp1),
      t_statistic = unname(test$statistic),
      df = unname(test$parameter),
      p_value = unname(test$p.value),
      stringsAsFactors = FALSE
    ))
  }
}

# Apply Bonferroni correction across all tests
test_results$p_adjusted <- p.adjust(test_results$p_value, method = "bonferroni")

# View results
print(test_results)
#write_tsv(finaldf,"ERTG_qPCR_finalresult.tsv")

#Comparing Stability of TG, ER####

#Reading in stability data

StabilityTable<-read_tsv("all_.99_stability.filtered.mx.txt")
StabilityTable<-StabilityTable%>%
  mutate(AveControl = C1+C2+C3)%>%
  mutate(AveKD = P1 + P2+P3)%>%
  mutate(RelativeStability=AveKD-AveControl)%>%
  mutate(Compartment = case_when(
    GeneID %in% sets_list$Horste.ER ~ "ER+",
    GeneID %in% sets_list$Ensembl.Mem  ~ "Mem",
    GeneID %in% sets_list$Horste.TG ~ "TG+",
    TRUE ~ "Others" ))
StabilityTable<-StabilityTable%>%filter(Compartment != "Mem")

table(StabilityTable$Compartment)

# Calculate medians for each group
medians <- StabilityTable %>%
  group_by(Compartment) %>%
  dplyr::summarize(median_value = median(RelativeStability, na.rm = TRUE),
                   n = n())%>%
  mutate(count = paste0("(", n,")"))

# Define colors to match the fill colors in the boxplot
median_colors <- c("Others" = "grey30", "ER+" = "#9D3F97", "TG+" = "#39B54A")

# Plot with medians
library(ggsignif)
e<-StabilityTable %>%
  ggplot(aes(x = factor(Compartment, levels = c('Others', 'ER+', 'TG+')), 
             y = RelativeStability)) +
  geom_boxplot(aes(fill = Compartment), varwidth = TRUE, outlier.shape = 1,outlier.size = 0.3, outlier.alpha = 0.5) +
  theme_linedraw() +
  geom_signif(comparisons = list(c("Others", "ER+")), map_signif_level = FALSE,
              y_position = 2.9, tip_length = 0.005, size = 0.15,
              textsize = 1.8) +
  geom_signif(comparisons = list(c("Others", "TG+")), map_signif_level = FALSE,
              y_position = 3.55, tip_length = 0.005, size = 0.15,
              textsize = 1.8) +
  coord_cartesian(ylim = c(-3.6, 4.1)) +
  scale_fill_manual(values = c("#9D3F97", "grey40","#39B54A")) +
  labs( x = "Compartments", y= "Relative mRNA stability") +
  theme(panel.grid.major = element_line(linetype = "blank"), 
        panel.grid.minor = element_line(linetype = "blank"),
        legend.position = "none",aspect.ratio = 2,
        axis.title=element_text(size=7), 
        axis.text.y=element_text(size=5),
        panel.background = element_blank(),
        plot.background = element_blank(),
        plot.margin = unit(c(0, 0, 0, 0), "cm")
        ) +
  # Add median values below each boxplot
  geom_text(data = medians, aes(x = Compartment, y = -3.53, 
                                label = sprintf("%.2f", median_value), 
                                color = Compartment), 
            inherit.aes = FALSE, size = 1.5, vjust = 1) +
  geom_text(data = medians, aes(x = Compartment, y = -3.75, label = sprintf("%s", count)),
           size = 1.5, vjust = 1, inherit.aes = FALSE, color = "black")+
  scale_color_manual(values = median_colors)
e



#plotting ER/TG foldchages of RIPiT#####

E117RvsWT<-read_csv("../DESeq2 for PYM samples/all_NT_res.csv")
E117RvsWT<-na.omit(E117RvsWT)
names(E117RvsWT)[1]<-"ensembl_gene_id"

E117RvsWT <- E117RvsWT %>%
  mutate(Class1 = case_when(
    ensembl_gene_id %in% sets_list$Horste.ER ~ "ER+",
    ensembl_gene_id %in% sets_list$Ensembl.Mem  ~ "Mem",
    ensembl_gene_id %in% sets_list$Horste.TG ~ "TG+",
    TRUE ~ "Others" ))
table(E117RvsWT$Class1)
E117RvsWT<-E117RvsWT%>%filter(Class1 != "Mem")

#plot:
# Calculate pairwise p-values (e.g., Wilcoxon test)
  p_values <- pairwise.wilcox.test(
    E117RvsWT$log2FoldChange, 
    E117RvsWT$Class1, 
    p.adjust.method = "bonferroni"
  )$p.value
  
  # Format p-value text
  pval_text <- paste("p values (Wilcoxon test):",
                     paste("ER+ vs TG+ =", format(p_values["TG+", "ER+"], scientific = TRUE, digits = 2)),
                     paste("ER+ vs Others =", format(p_values["Others", "ER+"], scientific = TRUE, digits = 2)),
                     paste("TG+ vs Others =", format(p_values["TG+","Others"], scientific = TRUE, digits = 2)),
                     sep = "\n")


# Calculate the number of entries in each class
library(dplyr)
class_counts <- E117RvsWT %>%
  group_by(Class1) %>%
  dplyr::summarize(count = n(), .groups = 'drop')

# Create formatted labels with counts
labels_with_counts <- paste(class_counts$Class1, 
                            " (", class_counts$count, ")", sep = "")

# Plotting with the updated dataframe and columns
c<-E117RvsWT %>%
  ggplot(aes(log2FoldChange, color = Class1)) +
  geom_vline(aes(xintercept = 0), linewidth = 0.02)+
  geom_hline(aes(yintercept = 0.5), linewidth = 0.02)+
  stat_ecdf(geom = "line", linewidth = 0.4) +
  scale_color_manual(values = c("ER+" = "#9D3F97", "TG+" = "#39B54A", "Others" = "grey60"), 
                     labels = labels_with_counts) +  # Change legend text
  coord_cartesian(xlim = c(-3, 3)) +
  theme_linedraw() +
  theme(
    axis.title = element_text(size=7),
    axis.text = element_text(size=5),
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    aspect.ratio = 1,
    legend.position = c(0.98, 0.02),         # Move legend to bottom right
    legend.justification = c(1, 0),
    legend.title=element_text(size=6), 
    legend.text=element_text(size=5),
    legend.key.size =unit(2, 'mm'))+     # Ensure correct alignment at bottom right
  labs(y = "Cumulative frequency", x = "log2FC E117R/WT",
       color = "Compartment") +
  annotate("text", x = -2, y = 0.85, label = pval_text, hjust = 0, size = 2, color = "black")+  # Add p-values on top left
  coord_cartesian(xlim = c(-2,2))
c


#slcn for ERTG######
#readslcn counts

E117RvsWTslcn<-read_csv("../SLCN_counts/slcn-trimmed-last-allcans-final_NT_res.csv")
names(E117RvsWTslcn)[1]<-"ensembl_gene_id"
E117RvsWTslcn<-E117RvsWTslcn%>%filter(baseMean != 0)
E117RvsWTslcn<-na.omit(E117RvsWTslcn)
E117RvsWTslcn$Geneid <- str_split(E117RvsWTslcn$ensembl_gene_id, "-", simplify = T)[,1]
E117RvsWTslcn$Region <- str_split(E117RvsWTslcn$ensembl_gene_id, "-", simplify = T)[,2]
table(E117RvsWTslcn$Region)
#Filter into groups based on ERTG
E117RvsWTslcn <- E117RvsWTslcn %>%
  mutate(Class1 = case_when(
    Geneid %in% sets_list$Horste.ER ~ "ER+",
    Geneid %in% sets_list$Ensembl.Mem  ~ "Mem",
    Geneid %in% sets_list$Horste.TG ~ "TG+",
    TRUE ~ "Others" ))

E117RvsWTslcn<-E117RvsWTslcn%>%filter(Class1 != "Mem")
E117RvsWTslcn<-E117RvsWTslcn%>%filter(Class1 != "Others")
E117RvsWTslcn<-E117RvsWTslcn%>%filter(Region != "s")
E117RvsWTslcn <- E117RvsWTslcn %>%
  mutate(Region = recode(Region,
                         "c" = "Canonical",
                         "n" = "Non-canonical",
                         "l" = "Last exon"))
UTR3<-read_tsv("TIS_vs_others_3primeUTR_50_all_genes.txt")


#Boxplots
median_values <- E117RvsWTslcn %>%
  group_by(Class1, Region) %>%
  dplyr::summarize(median_log2FC = median(log2FoldChange), .groups = 'drop',
                   n = n())%>%
  mutate(count = paste0("(", n,")"))
  
library(ggpubr)
d<-E117RvsWTslcn %>%
  ggplot(aes(x = Class1, y = log2FoldChange, fill = Class1)) +
  geom_hline(yintercept = 0, size = 0.05) +
  #geom_point(aes(color = Group), position = position_jitterdodge(), alpha = 0.5, size=0.5) +
  #geom_boxplot(position = position_dodge(), varwidth = T, alpha = 0.1, outlier.size = 0.8) +
  geom_boxplot(position = position_dodge(), varwidth = T, outlier.shape = 1,outlier.size = 0.3, outlier.alpha = 0.5)+
  scale_fill_manual(values = c("ER+" = "#9D3F97", "TG+" = "#39B54A")) +
  theme_linedraw() +
  coord_cartesian(ylim = c(-1.3,1.7))+
  theme( 
    axis.title = element_text(size=7),
    axis.text = element_text(size=5),
    panel.grid.major = element_blank(),
         panel.grid.minor = element_blank(),
         aspect.ratio = 2.1,
         legend.position = "none",
         strip.text = element_text(color="black",size =6, angle = 0, hjust =0),
         strip.background = element_blank(),
         strip.clip = "off",
         strip.placement = "inside") +
  labs( y = "log2FC E117R vs WT",
        x = "Region" ) +
  facet_wrap(~Region, nrow = 1) +
  stat_compare_means(
    comparisons = list(c("ER+", "TG+")),
    method = "wilcox.test",
    label = "p.format",
    label.y = 1.45,
    size = 2,
    tip.length = 0.01)+ 
  geom_text(data = median_values,
            aes(x = Class1, y = -1.1, label = sprintf("%.2f", median_log2FC), color = Class1),
            size = 1.5, vjust = 1.5, inherit.aes = FALSE) +
  geom_text(data = median_values,
            aes(x = Class1, y = -1.3, label = sprintf("%s", count)),
            size = 1.5, vjust = 1.5, inherit.aes = FALSE, color = "black")+
  scale_color_manual(values = c("ER+" = "#9D3F97", "TG+" = "#39B54A"))
  
d


UTR3<-read_tsv("TIS_vs_others_3primeUTR_50_all_genes.txt")
UTR3%>%ggplot(aes(log2foldchange, colour = group))+
  stat_ecdf()


library(patchwork)
layout <- c(
  area(1,1,2,2),area(1, 3, 2,6),
  area(3,1,4,2),area(3,3,4,5),area(3,6,4,6))
plot(layout)
a+b+c+d+e+plot_layout(design = layout)+
  plot_annotation(title = 'FIGURE 6', tag_levels = "a")  &
  theme(plot.tag = element_text(face = 'bold')) 
ggsave("F6patch1_3.pdf",
       a+b+c+d+e+
         theme(
           axis.text = element_text(size = 5),
           axis.title = element_text(size = 7),
           legend.title=element_text(size=6), 
           legend.text=element_text(size=5),
           legend.key=element_blank(),
           legend.box.background = element_blank())+
         plot_layout(design = layout)+
         plot_annotation(tag_levels = "a")  &
         theme(plot.tag = element_text(face = 'bold'),
               plot.background = element_blank(),
               plot.margin = unit(c(0, 0, 0, 0), "cm")),
       width = 180, 
       height = 180, 
       units = "mm")





 #Machine Learning plots######
#Kd vs Control

MLresults<-read_csv("../FinalML/permutation_model_adam_34_256_350_2e-06_240_plotting_info_testing_3split_KD.txt",col_names = F)
names(MLresults)<-c("Value", "Error","MLCorrelation","Correlation")
MLresults$Term<-rownames(MLresults)

# Create a vector with the feature names
MLresults$features <- as.factor(c("Gene length", "Exon number", "Exon length geometric mean", 
                                  "Intron length geometric mean", "Length of largest exon", 
                                  "Length of largest intron", "Total exon length", 
                                  "Total intron length", "Length of shortest exon", 
                                  "Length of shortest intron", "Mean exon length", 
                                  "Mean intron length", "Length of first intron", 
                                  "Length of largest exon excluding the last", 
                                  "Total exon length excluding the last", 
                                  "Exon length geometric mean excluding the last", 
                                  "Length of shortest exon excluding the last", 
                                  "Mean exon length excluding the last", 
                                  "Percentage GC content", "Number of transmembrane domains", 
                                  "Presence of signal peptide", "Stop codon:TAA", "Stop codon:TGA", "Stop codon:TAG", 
                                  "Length of 3\'UTR", "Length of 5\'UTR", "Number of exons starting in-frame", 
                                  "Number of exons starting in +1 frame", "Number of exons starting in -1 frame", 
                                  "Distance from AUG to exon junction", "Number of 5\'UTR junctions", 
                                  "Distance from exon junction to stop codon", 
                                  "Number of 3\'UTR junctions", "Average CDS exon length"))

correlation<-read_tsv("../FinalML/KD_total_correlation_coeffs.txt",col_names = F)
names(correlation)<-c("features","Correlation")
correlation<-correlation[1:34,]
correlation$features <- MLresults$features

# First plot
g <- MLresults %>%
  ggplot(aes(Value, reorder(features, Value))) +
  geom_vline(xintercept = 0, size = 0.5, linetype = 3, color = "grey") +
  geom_point() + theme_light() +
  geom_errorbar(aes(xmin = Value - Error, xmax = Value + Error), width = 0.2) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), aspect.ratio = 2,
        axis.text.y = element_text(size = 7),axis.text.x = element_text(size = 7))+#,
#        axis.title.x = element_text(hjust=-0.5)) +
  labs(y=NULL,title = "Features", x= "Relative importance:\nPYM1kd vs Control")
g

f<-correlation %>%
  ggplot(aes(y = reorder(features, Correlation), x = Correlation, fill = Correlation)) +
  geom_bar(stat = "identity", color = "black", size = 0.2) +
  theme_light() +
  scale_fill_gradient2(low = "red3", mid = "grey90", high = "blue3", midpoint = 0, limits = c(-0.3, 0.3)) +  # Set color gradient limits
  #scale_y_continuous(limits = c(-0.3, 0.3)) +  # Set y-axis range
  theme(
    panel.grid.major = element_blank(), panel.grid.minor = element_blank(), aspect.ratio =2,
    legend.position = c(0.05, 0.95), legend.justification = c(0, 1),
    legend.key.size = unit(0.5, "lines"), 
    legend.text = element_text(size = 6),
    legend.title=element_text(size=7),  # Adjust legend size
    axis.text.y = element_text(size = 7),axis.text.x = element_text(size = 7)
  ) +
  labs(y=NULL, x = "Correlation:\nPYM1kd vs Control", title = "Features",
       fill = "Correlation\nSign") #+
  #theme(axis.text.y = element_text(size = 10))
f




##Patchwork####




x<-free(f)+free(g)

plot(layout)
f+g
ggsave("S6patch2_3.pdf",
       f+g+
         plot_annotation(tag_levels = "A")  &
         theme(plot.tag = element_text(face = 'bold')),
       width = 8.5, 
       height = 11, 
       units = "in")



#plotting ERTG for PYM1 overexpression####

# Load required libraries
library(DESeq2)
library(tidyverse)
counts <- read_tsv("PYM_MANE_counts.tsv", col_names = TRUE, skip = 1)
counts <- counts %>%
  dplyr::select(Geneid, 7:18) %>% 
  column_to_rownames("Geneid")
colnames(counts) <- c("control1", "control2", "control3", 
                      "PYMkd1", "PYMkd2", "PYMkd3", 
                      "PYMOE1", "PYMOE2", "PYMOE3", 
                      "PYMN1", "PYMN2", "PYMN3")
sample_info <- data.frame(
  SampleID = colnames(counts),
  Condition = rep(c("Control", "PYMkd", "PYMOE", "PYMN"), each = 3)
)
sample_info$Condition <- factor(sample_info$Condition, levels = c("Control", "PYMkd", "PYMOE", "PYMN"))
keep <- rowSums(counts >= 10) >= 2
counts <- counts[keep, ]
dds <- DESeqDataSetFromMatrix(countData = counts,
                              colData = sample_info,
                              design = ~ Condition)
dds <- DESeq(dds)

results_pyoe_vs_control <- results(dds, contrast = c("Condition", "PYMOE", "Control"))

results_pyoe_vs_control_df <- as.data.frame(results_pyoe_vs_control)
results_pyoe_vs_control_df$GeneID <- rownames(results_pyoe_vs_control_df)

plotMA(results_pyoe_vs_control, ylim = c(-5, 5), main = "MA Plot: PYMOE vs Control")

names(results_pyoe_vs_control_df)[7]<-"ensembl_gene_id"

results_pyoe_vs_control_df<-results_pyoe_vs_control_df%>%mutate(Class1 = case_when(ensembl_gene_id %in% sets_list$Horste.ER ~ "ER+",
                                                                                   ensembl_gene_id %in% sets_list$Horste.TG ~ "TG+",
                                                                                   TRUE ~ "Others"))
table(results_pyoe_vs_control_df$Class1)

results_pyoe_vs_control_df <- results_pyoe_vs_control_df %>%
  mutate(Class1 = case_when(
    ensembl_gene_id %in% sets_list$Horste.ER ~ "ER+",
    ensembl_gene_id %in% sets_list$Ensembl.Mem  ~ "Mem",
    ensembl_gene_id %in% sets_list$Horste.TG ~ "TG+",
    TRUE ~ "Others" ))
table(results_pyoe_vs_control_df$Class1)
results_pyoe_vs_control_df<-results_pyoe_vs_control_df%>%filter(Class1 != "Mem")

#plot:
# Calculate pairwise p-values (e.g., Wilcoxon test)
p_values <- pairwise.wilcox.test(
  results_pyoe_vs_control_df$log2FoldChange, 
  results_pyoe_vs_control_df$Class1, 
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
class_counts <- results_pyoe_vs_control_df %>%
  group_by(Class1) %>%
  dplyr::summarize(count = dplyr::n(), .groups = 'drop')

# Create formatted labels with counts
labels_with_counts <- paste(class_counts$Class1, 
                            " (", class_counts$count, ")", sep = "")

# Plotting with the updated dataframe and columns
results_pyoe_vs_control_df %>%
  ggplot(aes(log2FoldChange, color = Class1)) +
  geom_vline(aes(xintercept = 0), linewidth = 0.02)+
  geom_hline(aes(yintercept = 0.5), linewidth = 0.02)+
  stat_ecdf(geom = "line", linewidth = 0.8) +
  scale_color_manual(values = c("ER+" = "red2", "TG+" = "blue2", "Others" = "grey75"), 
                     labels = labels_with_counts) +  # Change legend text
  coord_cartesian(xlim = c(-2, 2)) +
  theme_linedraw() +
  theme(
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    aspect.ratio = 1,
    legend.position = c(0.95, 0.05),         # Move legend to bottom right
    legend.justification = c(1, 0))+     # Ensure correct alignment at bottom right
  labs(title = "PYM1 OE vs Control fold changes",
       subtitle = "Grouped by Compartments",
       y = "Cumulative frequency", x = "Log2FC PYM1 OE/Control",
       color = "Compartment") +
  annotate("text", x = -1.9, y = 0.92, label = pval_text, hjust = 0, size = 4, color = "black")  # Add p-values on top left

results_pyoe_vs_control_df %>%
  ggplot(aes(log2FoldChange, color = Class1)) +
  geom_vline(aes(xintercept = 0), linewidth = 0.02)+
  geom_hline(aes(yintercept = 0.5), linewidth = 0.02)+
  stat_ecdf(geom = "line", linewidth = 0.8) +
  scale_color_manual(values = c("ER+" = "red2", "TG+" = "blue2", "Others" = "grey75"), 
                     labels = labels_with_counts) +  # Change legend text
  coord_cartesian(xlim = c(-2, 2)) +
  theme_linedraw() +
  theme(
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    aspect.ratio = 1,
    legend.position = c(0.95, 0.05),         # Move legend to bottom right
    legend.justification = c(1, 0))+     # Ensure correct alignment at bottom right
  labs(title = "PYM1 DelN OE vs Control fold changes",
       subtitle = "Grouped by Compartments",
       y = "Cumulative frequency", x = "Log2FC PYM1 DelN OE/Control",
       color = "Compartment") +
  annotate("text", x = -1.9, y = 0.92, label = pval_text, hjust = 0, size = 4, color = "black")  # Add p-values on top left

results_pyoe_vs_control_df %>%
  ggplot(aes(log2FoldChange, color = Class1)) +
  geom_vline(aes(xintercept = 0), linewidth = 0.02)+
  geom_hline(aes(yintercept = 0.5), linewidth = 0.02)+
  stat_ecdf(geom = "line", linewidth = 0.8) +
  scale_color_manual(values = c("ER+" = "red2", "TG+" = "blue2", "Others" = "grey75"), 
                     labels = labels_with_counts) +  # Change legend text
  coord_cartesian(xlim = c(-2, 2)) +
  theme_linedraw() +
  theme(
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    aspect.ratio = 1,
    legend.position = c(0.95, 0.05),         # Move legend to bottom right
    legend.justification = c(1, 0))+     # Ensure correct alignment at bottom right
  labs(title = "PYM1 vs DelN fold changes",
       subtitle = "Grouped by Compartments",
       y = "Cumulative frequency", x = "Log2FC PYM1 /PYM1delN",
       color = "Compartment") +
  annotate("text", x = -1.9, y = 0.92, label = pval_text, hjust = 0, size = 4, color = "black")  # Add p-values on top left

