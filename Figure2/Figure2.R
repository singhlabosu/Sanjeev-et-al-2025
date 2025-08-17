setwd("~/OneDrive - The Ohio State University/Singhlab/Bioinfo/PYMpaperFiguresUpdate/Figure2/")
library(tidyverse)

CHXvsE117R<-read_csv("../DESeq2 for PYM samples/all_mutant_res.csv")
CHXvsE117R<-na.omit(CHXvsE117R)
names(CHXvsE117R)[1]<-"V1"
CHXvsE117R<-CHXvsE117R%>%mutate(log2FoldChange=-1*log2FoldChange)

CHXvsWT<-read_csv("../DESeq2 for PYM samples/all_WT_res.csv")
CHXvsWT<-na.omit(CHXvsWT)
names(CHXvsWT)[1]<-"V1"
CHXvsWT<-CHXvsWT%>%mutate(log2FoldChange=-1*log2FoldChange)

FoldchangeComp<-inner_join(CHXvsWT,CHXvsE117R, by = "V1")

library(ggpmisc)
library(ggpubr)
library(ggpointdensity)
reg<-lm(formula = log2FoldChange.y ~ log2FoldChange.x, 
        data=FoldchangeComp)                       

#get intercept and slope value 
coeff<-coefficients(reg)           
intercept<-coeff[1] 
slope<- coeff[2] 

a<-FoldchangeComp%>%filter(baseMean.x>100)%>%
  ggplot(aes(log2FoldChange.x,log2FoldChange.y))+
  geom_pointdensity(shape = 1, alpha = 0.8)+
  scale_colour_gradient(low = "grey70", high = "black", na.value = NA)+
  theme_linedraw()+
  theme(axis.text = element_text(size = 5),axis.title = element_text(size = 7),
        aspect.ratio = 1,panel.grid.major = element_line(linetype = "blank"),
        panel.grid.minor = element_line(linetype = "blank"))+
  coord_cartesian(xlim = c(-5,5), ylim = c(-3.5,3.5))+
  labs(#title = "Comparison of foldchanges upon CHX treatment", subtitle = "WT vs E117R", caption = "baseMean > 100",
       x = "log2FC WT , +/- CHX",  y = "log2FC E117R , +/- CHX")+
  stat_cor(aes(label = after_stat(rr.label)),label.x = -4.5, label.y = 3.2, size = 2)+
  geom_abline(intercept = intercept, slope = slope, color = "black", linewidth = 0.2, linetype="dashed")+ 
  theme(legend.position = "none")
a
FoldchangeComp%>%filter(log2FoldChange.x>2.5 & log2FoldChange.y < -1)%>%dplyr::select(V1)

E117RvsWT<-read_csv("../DESeq2 for PYM samples/all_NT_res.csv")
names(E117RvsWT)[1]<-"V1"
E117RvsWT<-na.omit(E117RvsWT)
CombinedScatterdf<-inner_join(CHXvsWT,E117RvsWT, by = "V1")

reg<-lm(formula = log2FoldChange.y ~ log2FoldChange.x, 
        data=CombinedScatterdf)                       

#get intercept and slope value 
coeff<-coefficients(reg)           
intercept<-coeff[1] 
slope<- coeff[2]

b<-CombinedScatterdf%>%ggplot(aes(log2FoldChange.x,log2FoldChange.y))+
  geom_pointdensity(shape = 1, alpha = 0.8)+
  scale_colour_gradient(low = "grey70", high = "black", na.value = NA)+
  theme_linedraw()+
  theme(axis.text = element_text(size = 5),axis.title = element_text(size = 7),
        aspect.ratio = 1,panel.grid.major = element_line(linetype = "blank"),
        panel.grid.minor = element_line(linetype = "blank"))+
  coord_cartesian(xlim = c(-5,5), ylim = c(-5,5))+
  labs(#title = "Comparison of foldchanges", subtitle = "CHX vs MAGOHE117R",
       x = "log2FC CHX vs WT",  y = "log2FC E117R vs WT")+
  stat_cor(aes(label = after_stat(rr.label)),label.x = -4.5, label.y = 4.8, size=2)+
  geom_abline(intercept = intercept, slope = slope, color = "black", linewidth = 0.2, linetype="dashed")+ 
  theme(legend.position = "none")
b


library(biomaRt)
listMarts(host='https://apr2020.archive.ensembl.org')
ensembl100=useMart(host='https://apr2020.archive.ensembl.org', biomart='ENSEMBL_MART_ENSEMBL', dataset='hsapiens_gene_ensembl')
attributes <- listAttributes(ensembl100)
filters <- listFilters(ensembl100)

lncRNA<- getBM(attributes = c("ensembl_gene_id","gene_biotype"),
               mart = ensembl100,
               useCache = FALSE)
table(lncRNA$gene_biotype)
lncRNA<-lncRNA%>%filter(gene_biotype=="lncRNA")


##plotting cdfs for RP mRNAs and lncRNAs
names(E117RvsWT)[1]<-"ensembl_gene_id"
E117RvsWT<-E117RvsWT%>%mutate(Type = case_when(ensembl_gene_id %in% lncRNA$ensembl_gene_id ~ "lncRNA",
                                               TRUE ~ "Protein_coding"))

##PLotting lncRNA #####
# Perform the Wilcoxon rank-sum test for 'Type' column
wilcox_test <- wilcox.test(log2FoldChange ~ Type, data = E117RvsWT,
                           exact = FALSE, alternative = "greater")

# Extract p-value and gene counts based on 'Type'
p_value <- formatC(wilcox_test$p.value, format = "e", digits = 2)  # Format p-value in scientific notation
count_type1 <- sum(E117RvsWT$Type == unique(E117RvsWT$Type)[1])
count_type2 <- sum(E117RvsWT$Type == unique(E117RvsWT$Type)[2])

# Modify legend text to include gene counts
legend_labels <- c(
  paste0(unique(E117RvsWT$Type)[1], "\n(n=", count_type1, ")"),
  paste0(unique(E117RvsWT$Type)[2], "\n(n=", count_type2, ")"))

legend_labels<-str_replace(legend_labels, "\\_", "\n")
# Ensure that 'Type' levels are in the intended order
E117RvsWT$Type <- factor(E117RvsWT$Type, levels = unique(E117RvsWT$Type))

# Plot with p-value annotation and modified legend text based on 'Type'
c<-E117RvsWT %>%
  ggplot(aes(log2FoldChange, color = Type)) +
  stat_ecdf(geom = "line", linewidth = 0.5) +
  coord_cartesian(xlim = c(-2, 2)) +
  scale_color_manual(values = c("grey40","firebrick"), labels = legend_labels) +  # Update legend labels
  theme_linedraw() +
  geom_vline(aes(xintercept = 0), linewidth = 0.02) +
  geom_hline(aes(yintercept = 0.5), linewidth = 0.02) +
  theme(axis.text = element_text(size = 5),axis.title = element_text(size = 7),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        aspect.ratio = 1,
        legend.position = c(0.99, 0.01),
        legend.justification = c(1, 0),
        legend.spacing.y = unit(5, 'mm'),
        legend.title=element_text(size=6), 
        legend.text=element_text(size=5),
        legend.key=element_blank(),
        legend.box.background = element_blank(),
        legend.background = element_blank(),
        legend.key.size =unit(2, 'mm')) +
  labs(#title = "E117R vs WT foldchanges - lncRNAs vs others",
       y = "Cumulative frequency",
       x = "log2FC E117R/WT",
       color = "Type" ) +  # Label for the legend
  annotate("text",x = -1.6, y = 0.95,label = paste("p =", p_value),hjust = 0,size = 2)
c
compare_means(log2FoldChange ~ Type, data = E117RvsWT, method = "wilcox.test")

#####



###Reading NucCytoRanks
NucCyto<-read_tsv("nuclear fractions and percentiles.txt")
names(NucCyto)[1]<-"ensembl_gene_id"


lncRNA_NucCyto<-inner_join(lncRNA,NucCyto)
glimpse(lncRNA_NucCyto)
lncRNA_NucCyto$Fraction<-as.numeric(lncRNA_NucCyto$`Nuclear Fraction`)
summary(lncRNA_NucCyto$Fraction)

lncRNA_NucCyto<-lncRNA_NucCyto%>%mutate(NuclearlncRNA = case_when(Fraction > 0.9580 ~ "Yes",
                                                                  TRUE ~ "No"))
table(lncRNA_NucCyto$NuclearlncRNA)
names(lncRNA_NucCyto)

names(E117RvsWT)[1]<-"ensembl_gene_id"

lncRNA_NucCyto_E117R<-inner_join(lncRNA_NucCyto,E117RvsWT)
table(lncRNA_NucCyto_E117R$NuclearlncRNA)

# Perform the Wilcoxon rank-sum test for 'Type' column
wilcox_test <- wilcox.test(log2FoldChange ~ NuclearlncRNA, data = lncRNA_NucCyto_E117R,
                           exact = FALSE, alternative = "less")

# Extract p-value and gene counts based on 'Type'
p_value <- formatC(wilcox_test$p.value, format = "e", digits = 2)  # Format p-value in scientific notation
count_type1 <- sum(lncRNA_NucCyto_E117R$NuclearlncRNA == unique(lncRNA_NucCyto_E117R$NuclearlncRNA)[1])
count_type2 <- sum(lncRNA_NucCyto_E117R$NuclearlncRNA == unique(lncRNA_NucCyto_E117R$NuclearlncRNA)[2])

# Modify legend text to include gene counts
legend_labels <- c(
  paste0(unique(lncRNA_NucCyto_E117R$NuclearlncRNA)[1], " (n=", count_type1, ")"),
  paste0(unique(lncRNA_NucCyto_E117R$NuclearlncRNA)[2], " (n=", count_type2, ")"))

# Ensure that 'Type' levels are in the intended order
lncRNA_NucCyto_E117R$NuclearlncRNA <- factor(lncRNA_NucCyto_E117R$NuclearlncRNA, levels = unique(lncRNA_NucCyto_E117R$NuclearlncRNA))

d<-lncRNA_NucCyto_E117R%>%
  ggplot(aes(log2FoldChange, color = NuclearlncRNA))+
  stat_ecdf(geom = "line", linewidth = 0.5)+
  coord_cartesian(xlim = c(-2,2))+
  scale_color_manual(values=c("red", "red4"), labels = legend_labels) +
  theme_linedraw()+
  geom_vline(aes(xintercept = 0), linewidth = 0.02)+
  geom_hline(aes(yintercept = 0.5), linewidth = 0.02)+
  theme(axis.text = element_text(size = 5),axis.title = element_text(size = 7),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        aspect.ratio = 1,
        legend.position = c(0.99, 0.01),
        legend.justification = c(1, 0),
        legend.title=element_text(size=6), 
        legend.text=element_text(size=5),
        legend.key=element_blank(),
        legend.box.background = element_blank(),
        legend.background = element_blank(),
        legend.key.size =unit(2, 'mm'))+
  annotate("text",x = -1.6, y = 0.95,label = paste("p =", p_value),hjust = 0,size = 2)+
  labs(#title = "E117R vs WT foldchanges - lncRNA by nuclear fraction",
       y = "Cumulative frequency", x = "log2FC E117R/WT", color ="Nuclear\nlncRNA")

d
wilcox_test
compare_means(log2FoldChange ~ NuclearlncRNA, data = lncRNA_NucCyto_E117R, method = "wilcox.test")
library(patchwork)
(b+a)/(c+d)

layout <- c(
  area(1,1,2,2),area(1,3,2,4),
  area(3,1,4,2),area(3,3,4,4))
plot(layout)

ggsave("Figure2patch.pdf",
       b+a+c+d+plot_layout(design = layout)+
         plot_annotation(title = 'Figure2', tag_levels = "a")  &
         theme(plot.tag = element_text(face = 'bold')),
       width = 8.5, 
       height = 11, 
       units = "in")

ggsave("Figure2patch_final.pdf",
       b+a+c+d+plot_layout(design = layout)+
         plot_annotation(tag_levels = "a")  &
         theme(plot.tag = element_text(face = 'bold'),
               panel.background = element_blank(),   
               plot.background = element_blank(),
               plot.margin = unit(c(0, 0, 0, 0), "cm")),
       width = 88, 
       height = 90, 
       units = "mm")
