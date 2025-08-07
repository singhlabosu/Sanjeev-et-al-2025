###Figure4 plots

library(tidyverse)
library(ggrepel)
library(patchwork)
library(tximport)
library(DESeq2)
library(ggsignif)
  
setwd("~/OneDrive - The Ohio State University/Singhlab/Bioinfo/PYMpaperFiguresUpdate/Figure4/")
#Exon. numbers plot####
###reading files in 
E117RvsWT<-read_csv("../DESeq2 for PYM samples/all_NT_res.csv")
E117RvsWT<-na.omit(E117RvsWT)
names(E117RvsWT)[1]<-"ensembl_gene_id"
SEgenes_MANE_noInt<-read_tsv("Kovalak_etal/SEgenes_MANE_noInt.tsv")

###Getting exon numbers from ensembl
Exons<-read_tsv("../MANEexons.tsv")

E117RvsWT<-inner_join(E117RvsWT, Exons[,c(1,4)])
E117RvsWT<-E117RvsWT%>%mutate(ExonGroup = case_when(rank > 1 ~ "Multi",
                                                    rank == 1~ "Single"))
E117RvsWT$ExonGroup
table(E117RvsWT$ExonGroup)
E117RvsWT <- E117RvsWT %>%
  filter(!(ExonGroup == "Single" & !ensembl_gene_id %in% SEgenes_MANE_noInt$ensembl_gene_id))


# calculate pairwise p-values (for example, Wilcoxon test)
p_values <- pairwise.wilcox.test(
  E117RvsWT$log2FoldChange, 
  E117RvsWT$ExonGroup, 
  p.adjust.method = "bonferroni"
)$p.value

test<-wilcox.test(log2FoldChange ~ ExonGroup, data = E117RvsWT,
            exact = FALSE, alternative = "two.sided")
test$p.value

# Convert p-values to a readable format for display
# Extract and format p-values correctly based on row and column names
pval_text <- paste("p-value","(Wilcoxon test):",
                   paste(format(p_values, scientific = TRUE, digits = 2)),  # "c vs l"
                   sep = "\n")

# Calculate the number of genes in each group
group_counts <- E117RvsWT %>%
  group_by(ExonGroup) %>%
  dplyr::summarize(count = n(), .groups = 'drop')



# Create formatted labels with counts
labels_with_counts <- paste(group_counts$ExonGroup, 
                            " (", group_counts$count, ")", sep = "")


# plotting as before
a<-E117RvsWT %>%
  #filter(Region != "s") %>%
  ggplot(aes(log2FoldChange, color = ExonGroup)) +
  geom_vline(aes(xintercept = 0), linewidth = 0.02)+
  geom_hline(aes(yintercept = 0.5), linewidth = 0.02)+
  stat_ecdf(geom = "line", linewidth = 0.6) +
  scale_color_manual(values = c("Multi" = "grey45", "Single" = "#27AAE1"), 
                     labels = labels_with_counts) +  # Change legend text
  coord_cartesian(xlim = c(-2, 2)) +
  theme_linedraw() +
  theme(
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    aspect.ratio = 1,
    legend.position = c(0.95, 0.05),         # Move legend to bottom right
    legend.justification = c(1, 0),          # Ensure correct alignment at bottom right
    axis.text = element_text(size = 5),
    axis.title = element_text(size = 7),
    legend.title=element_text(size=6), 
    legend.text=element_text(size=5),
    legend.key=element_blank(),
    legend.key.size =unit(2, 'mm'),
    legend.box.background = element_blank()
  ) +
  labs(y = "Cumulative frequency", x = "log2FC E117R/WT") +
  annotate("text", x = 0.9, y = 0.32, label = pval_text, hjust = 0, size = 2, color = "black")  # Add p-values on top left
a+labs(title = "E117R vs WT fold changes")

##Volcano plot with SEgenes####
E117RvsWT<-E117RvsWT%>%mutate(log10padj=log10(padj))
E117RvsWT%>%
  ggplot(aes(log2FoldChange,-log10padj))+
  geom_point(data=subset(E117RvsWT, ExonGroup == "Multi"),
             color="gray60",size = 3.0, alpha= 0.2)+
  geom_point(data=subset(E117RvsWT, ExonGroup == "Single"),
             color="blue3",alpha = 0.3, size = 3.0)+
  coord_cartesian(xlim =c(-5,5), ylim=c(0,13))+
  theme_linedraw()+
  theme(aspect.ratio = 6/14)+
  labs(title = "E117R vs WT") + theme(panel.grid.major = element_line(linetype = "blank"),
                                      panel.grid.minor = element_line(linetype = "blank"),
                                      axis.title = element_text(size = 15),
                                      axis.text = element_text(size = 15),
                                      plot.title = element_text(size = 11))

##slcn plot####

# Reading in files
E117RvsWTslcn <- read_csv("../SLCN_counts/slcn-trimmed-last-allcans-final_NT_res.csv")
names(E117RvsWTslcn)[1] <- "ensembl_gene_id"
E117RvsWTslcn <- E117RvsWTslcn %>% filter(baseMean > 10)

# Extracting Gene ID and Region
E117RvsWTslcn$Geneid <- str_split(E117RvsWTslcn$ensembl_gene_id, "-", simplify = T)[, 1]
E117RvsWTslcn$Region <- str_split(E117RvsWTslcn$ensembl_gene_id, "-", simplify = T)[, 2]

# Check counts per region
table(E117RvsWTslcn$Region)
E117RvsWTslcn <- E117RvsWTslcn %>%
  filter(!(Region == "s" & !Geneid %in% SEgenes_MANE_noInt$ensembl_gene_id))

# E117RvsWTslcn <- E117RvsWTslcn %>%
#   filter(!(Region == "l" & !Geneid %in% NoInt3UTR$ensembl_gene_id))
# table(E117RvsWTslcn$Region)

# Pairwise Wilcoxon test (include all regions, including "s")
p_values <- pairwise.wilcox.test(
  E117RvsWTslcn$log2FoldChange, 
  E117RvsWTslcn$Region, 
  p.adjust.method = "bonferroni"
)$p.value

# Convert p-values to readable format for all comparisons
pval_text <- paste("p-values (Wilcoxon test):",
                   paste("c vs l =", format(p_values["l", "c"], scientific = TRUE, digits = 2)),  # c vs l
                   paste("c vs n =", format(p_values["n", "c"], scientific = TRUE, digits = 2)),  # c vs n
                   paste("l vs n =", format(p_values["n", "l"], scientific = TRUE, digits = 2)),  # l vs n
                   paste("c vs s =", format(p_values["s", "c"], scientific = TRUE, digits = 2)),  # c vs s
                   paste("l vs s =", format(p_values["s", "l"], scientific = TRUE, digits = 2)),  # l vs s
                   paste("n vs s =", format(p_values["s", "n"], scientific = TRUE, digits = 2)),  # n vs s
                   sep = "\n")

# Map region codes to descriptive labels
E117RvsWTslcn <- E117RvsWTslcn %>%
  mutate(Regions = recode(Region,
                          "c" = "canonical",
                          "l" = "last exon",
                          "n" = "non-canonical",
                          "s" = "single exon"))

# Calculate the number of entries in each region
region_counts <- E117RvsWTslcn %>%
  group_by(Regions) %>%
  dplyr::summarize(count = n(), .groups = 'drop')

# Create formatted labels with counts
labels_with_counts <- paste(region_counts$Regions, 
                            " (", region_counts$count, ")", sep = "")

library(ggtext)

region_colors <- c(
  c = "#A97C50",  # Canonical (from Set1)
  l = "#F15A29",  # Last exon (from Set1)
  n = "#FBB040",  # Non-canonical (from Set1)
  s = "#27AAE1"   # Single exon (from Set1)
)

# Create pval_text with colored region symbols (c, l, n, s)
pval_text <- paste(
  "p-values (Wilcoxon test):",
  paste0("<span style='color:", region_colors["c"], "'>c</span> vs ",
         "<span style='color:", region_colors["l"], "'>l</span> = ",
         format(p_values["l", "c"], scientific = TRUE, digits = 2)),
  paste0("<span style='color:", region_colors["c"], "'>c</span> vs ",
         "<span style='color:", region_colors["n"], "'>n</span> = ",
         format(p_values["n", "c"], scientific = TRUE, digits = 2)),
  paste0("<span style='color:", region_colors["l"], "'>l</span> vs ",
         "<span style='color:", region_colors["n"], "'>n</span> = ",
         format(p_values["n", "l"], scientific = TRUE, digits = 2)),
  paste0("<span style='color:", region_colors["c"], "'>c</span> vs ",
         "<span style='color:", region_colors["s"], "'>s</span> = ",
         format(p_values["s", "c"], scientific = TRUE, digits = 2)),
  paste0("<span style='color:", region_colors["l"], "'>l</span> vs ",
         "<span style='color:", region_colors["s"], "'>s</span> = ",
         format(p_values["s", "l"], scientific = TRUE, digits = 2)),
  paste0("<span style='color:", region_colors["n"], "'>n</span> vs ",
         "<span style='color:", region_colors["s"], "'>s</span> = ",
         format(p_values["s", "n"], scientific = TRUE, digits = 2)),
  sep = "<br>"
)

# Add caption to the plot
b<-E117RvsWTslcn %>%
  ggplot(aes(log2FoldChange, color = Region)) +
  stat_ecdf(geom = "line", linewidth = 0.6) +
  scale_color_manual(
    values = c(c = "#A97C50",  
               l = "#F15A29", 
               n = "#FBB040", 
               s = "#27AAE1"), 
    labels = labels_with_counts
  ) +
  coord_cartesian(xlim = c(-2, 2)) +
  theme_linedraw() +
  theme(plot.tag = element_markdown(size = 10),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    aspect.ratio = 1,
    legend.position = c(0.99, 0.01),
    legend.justification = c(1, 0),
    axis.text = element_text(size = 5),
    axis.title = element_text(size = 7),
    legend.title=element_text(size=6), 
    legend.text=element_text(size=5),
    legend.key=element_blank(),
    legend.key.size =unit(2, 'mm'),
    legend.box.background = element_blank(),
    plot.caption = element_markdown(size = 5)  # Enable markdown for caption
  ) +
  labs(y = "Cumulative frequency", x = "log2FC E117R/WT",
    caption = pval_text  # Add styled p-value text as caption
  )
b+labs(title = "E117R vs WT fold changes",
  subtitle = "Grouped by EJC regions")

# Prepare variables
data <- E117RvsWTslcn
group_var <- "Region"
value_var <- "log2FoldChange"

# Get all pairwise combinations
group_levels <- unique(data[[group_var]])
pairs <- combn(group_levels, 2, simplify = FALSE)

# Perform tests and store results
results <- lapply(pairs, function(pair) {
  group1 <- pair[1]
  group2 <- pair[2]
  
  x <- data[[value_var]][data[[group_var]] == group1]
  y <- data[[value_var]][data[[group_var]] == group2]
  
  test <- wilcox.test(x, y, exact = FALSE)
  
  data.frame(
    group1 = group1,
    group2 = group2,
    p_value = test$p.value,
    W_statistic = test$statistic,
    alternative = test$alternative,
    method = test$method,
    stringsAsFactors = FALSE
  )
})

# Combine into a dataframe
results_df <- bind_rows(results)

# Adjust p-values (Bonferroni, same as in test)
results_df <- results_df %>%
  mutate(p_adjusted = p.adjust(p_value, method = "bonferroni"),
         hypothesis = paste(group1, "vs", group2))

# Reorder columns
results_df <- results_df %>%
  dplyr::select(hypothesis, group1, group2, W_statistic, p_value, p_adjusted, alternative, method)

# pval_df <- as.data.frame(as.table(stat$p.value)) %>%
#   filter(!is.na(Freq)) %>%
#   dplyr::rename(group1 = Var1, group2 = Var2, p_adjusted = Freq)

##Making plot for SE genes

######
#RNAseq vs RIPiT

PYMkdCounts<-read_tsv("Counts/PYMMANEcounts.txt",skip = 1)
HauerClipCounts<-read_tsv("Counts/HauerMANEcounts.txt",skip = 1)
WT_NT1<-read_tsv("Counts/WT_all_NT_1", skip = 1)
WT_NT2<-read_tsv("Counts/WT_all_NT_2", skip = 1)

CombinedCounts<-inner_join(PYMkdCounts[,c(1,7:9)], HauerClipCounts[,c(1,10:12)])
CombinedCounts<-inner_join(CombinedCounts, WT_NT1[,c(1,7)])
CombinedCounts<-inner_join(CombinedCounts, WT_NT2[,c(1,7)])
names(CombinedCounts)
CombinedCounts$RNA<-rowMeans(CombinedCounts[,2:4])
CombinedCounts$Clip<-rowMeans(CombinedCounts[,5:7])
CombinedCounts$RIPiT<-rowMeans(CombinedCounts[,8:9])
CombinedCounts <- CombinedCounts[rowSums(CombinedCounts[, 2:9]) != 0, ]


#Plotting RNA vs RIPiT
#getting SEgenes
SEgeneslist<-Exons%>%filter(rank == 1)%>%dplyr::select(ensembl_gene_id)
#ENSG00000286522(H3C2),ENSG00000187837(H1-2)
CombinedCounts<-CombinedCounts%>%mutate(Name= case_when(Geneid == "ENSG00000286522" ~ "H3C2",
                                                        Geneid == "ENSG00000187837" ~ "H1-2",
                                                        TRUE ~ ""))

SEgenesScatter<-CombinedCounts%>%filter(RNA != 0 & RIPiT != 0)%>%filter(Geneid %in% SEgeneslist$ensembl_gene_id)
library(scales)
reg<-lm(formula = RIPiT ~ RNA, 
        data=SEgenesScatter)                       

library(ggpubr)
library(ggrepel)
c<-SEgenesScatter%>%
  ggplot(aes(RNA,RIPiT, label = Name))+
  geom_point(shape = 1, alpha =0.8, color="black", size=0.8)+geom_label_repel(size = 2)+
  geom_point(data = SEgenesScatter[SEgenesScatter$Name != "",], color = "red")+
  theme_linedraw()+
  labs(x= "RNA-seq",y="RIPiT-seq")+
  theme(aspect.ratio = 1,
        axis.text = element_text(size = 5),
        axis.title = element_text(size = 7),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = c(0.99, 0.01),
        legend.justification = c(1, 0),
        legend.title=element_text(size=6), 
        legend.text=element_text(size=5),
        legend.key=element_blank(),
        legend.box.background = element_blank())+
  stat_cor(aes(label = after_stat(rr.label)),label.x = 0, label.y = 2.7)+
  geom_smooth(method = lm, color = "grey", linewidth = 0.4, linetype="dashed", se = F)+
  scale_x_log10()+scale_y_log10()
c+labs(title = "SE genes scatter: RNA vs RIPiT")

SEgeneswithin2X<-read_tsv("SEgeneswithin2X.tsv")
d<-SEgenesScatter%>%filter(Clip != 0)%>%filter(Geneid %in% SEgeneswithin2X$ensembl_gene_id)%>%
  ggplot(aes(Clip,RIPiT, label = Name))+
  geom_point(shape = 1, alpha =0.8, color="black",size=0.8)+geom_label_repel(size = 2)+
  geom_point(data = SEgenesScatter[SEgenesScatter$Name != "",], color = "red")+
  theme_linedraw()+
  labs(x= "CLIP-seq",y="RIPiT-seq")+
  theme(aspect.ratio = 1,
        axis.text = element_text(size = 5),
        axis.title = element_text(size = 7),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = c(0.99, 0.01),
        legend.justification = c(1, 0),
        legend.title=element_text(size=6), 
        legend.text=element_text(size=5),
        legend.key=element_blank(),
        legend.box.background = element_blank())+
  stat_cor(aes(label = after_stat(rr.label)),label.x = 0, label.y = 2.7)+
  geom_smooth(method = lm, color = "grey", linewidth = 0.4, linetype="dashed", se = F)+
  scale_x_log10()+scale_y_log10()
d+ labs(title = "SE genes scatter: Clip vs RIPiT")

SEgenesScatter%>%filter(Clip > 10 & RIPiT >10)%>%filter(Geneid %in% SEgeneswithin2X$ensembl_gene_id)%>%
  ggplot(aes(Clip,RIPiT, label = Name))+
  geom_point(shape = 1, alpha =0.8, color="black",size=0.8)+geom_label_repel(size = 2)+
  geom_point(data = SEgenesScatter[SEgenesScatter$Name != "",], color = "red")+
  theme_linedraw()+
  labs(x= "CLIP-seq",y="RIPiT-seq")+
  theme(panel.grid.major = element_line(linetype = "blank"),
        panel.grid.minor = element_line(linetype = "blank"),
        aspect.ratio = 1)+
  stat_cor(aes(label = after_stat(rr.label)),label.x = 1, label.y = 2.7)+
  geom_smooth(method = lm, color = "grey", linewidth = 0.4, linetype="dashed", se = F)+
  scale_x_log10()+scale_y_log10()

library(patchwork)
b+c

layout <- c(
  area(1,1,2,2),area(3,1,4,2),area(1, 3, 2,4),
  area(1,5,2,6),area(1,1))
plot(layout)
a+b+c+d+plot_layout(design = layout)+
  plot_annotation(title = 'FIGURE 4', tag_levels = "a")  

ggsave("test.pdf",
       a+b+c+d+plot_layout(design = layout)+
         plot_annotation(tag_levels = "a") &
         theme(panel.background = element_blank(),   
               plot.background = element_blank(),
               plot.margin = unit(c(0, 0, 0, 0), "cm")),
       width = 180,
       height = 120,
       units = "mm")

final_plot <- (a + b + c + d) +
  plot_layout(design = layout) +
  plot_annotation(
    tag_levels = "a",
    theme = theme(
      plot.tag = element_text(face = 'bold'),
      panel.background = element_blank(),
      plot.background = element_blank(),
      plot.margin = unit(c(0, 0, 0, 0), "cm")
    )
  )

ggsave("F4patch3.pdf",
       plot = final_plot,
       width = 180,
       height = 120,
       units = "mm")



###reading files in 

E117RvsWT<-E117RvsWT%>%mutate(ExonGroup = case_when(rank > 9 ~ "10+",
                                                    TRUE ~ as.character(rank)))
E117RvsWT$ExonGroup <- factor(E117RvsWT$ExonGroup, levels = c(as.character(1:9), "10+"))

library(RColorBrewer)  

E117RvsWT %>%
  ggplot(aes(log2FoldChange, color = ExonGroup)) +
  stat_ecdf(geom = "line", linewidth = 0.8) +
  coord_cartesian(xlim = c(-2, 2)) +
  scale_color_manual(values = c("#FFFFCC",  # Light Yellow
                                "#FFEDA0",  # Pale Yellow
                                "#FED976",  # Light Orange-Yellow
                                "#FEB24C",  # Soft Orange
                                "#FD8D3C",  # Medium Orange
                                "#FC4E2A",  # Strong Orange-Red
                                "#E31A1C",  # Bright Red
                                "#BD0026",  # Dark Red
                                "#800026",  # Deep Burgundy Red
                                "#67000D"))+# Darkest Red
  theme_linedraw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        aspect.ratio = 1) +
  labs(y = "Cumulative frequency", x = "log2FC E117R/WT")

E117RvsWT %>%
  ggplot(aes(log2FoldChange, color = ExonGroup)) +
  stat_ecdf(geom = "line", linewidth = 0.8) +
  coord_cartesian(xlim = c(-2, 2)) +
  scale_x_continuous(n.breaks = 3)+
  scale_y_continuous(n.breaks = 3)+
  scale_color_manual(values = c("#27AAE1","#FED976", "#FEBB50", "#FD9E3B", 
                                "#FC7F32","#E65A29", "#CC3722", "#AA1E1B",
                                "#870912","#67000D"))+
  theme_linedraw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        aspect.ratio = 1,
        legend.position = c(0.99, 0.01),
        legend.justification = c(1, 0),
        legend.title=element_text(size=5), 
        legend.text=element_text(size=4),
        legend.key.size =unit(2, 'mm'),
        legend.box.background = element_blank()) +
  labs(y = "Cumulative frequency", x = "log2FC E117R/WT")

E117RvsWT %>%
  ggplot(aes(log2FoldChange, color = ExonGroup)) +
  stat_ecdf(geom = "line", linewidth = 0.8) +
  coord_cartesian(xlim = c(-2, 2)) +
  scale_color_manual(values = c("#27AAE1","#FED976", "#FEBB50", "#FD9E3B", 
                                "#FC7F32","#E65A29", "#CC3722", "#AA1E1B",
                                "#870912","#67000D"))+
  theme_linedraw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        aspect.ratio = 1,
        legend.position = c(0.99, 0.01),
        legend.justification = c(1, 0),
        legend.title=element_text(size=5), 
        legend.text=element_text(size=4),
        legend.key.size =unit(2, 'mm'),
        legend.box.background = element_blank()) +
  labs(y = "Cumulative frequency", x = "log2FC E117R/WT")

library(MetBrewer)
colorblind_palettes
E117RvsWT %>%
  ggplot(aes(log2FoldChange, color = ExonGroup)) +
  stat_ecdf(geom = "line", linewidth = 0.8) +
  coord_cartesian(xlim = c(-2, 2)) +
  scale_color_manual(values=met.brewer("Hiroshige"))+
  theme_linedraw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        aspect.ratio = 1) +
  labs(y = "Cumulative frequency", x = "log2FC E117R/WT")

E117RvsWT %>%
  ggplot(aes(log2FoldChange, color = ExonGroup)) +
  stat_ecdf(geom = "line", linewidth = 0.8) +
  coord_cartesian(xlim = c(-2, 2)) +
  scale_color_manual(values=rev(met.brewer("Benedictus",n=10)))+
  theme_linedraw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        aspect.ratio = 1) +
  labs(y = "Cumulative frequency", x = "log2FC E117R/WT")

library(scales)

# Define the start and end colors
color_gradient <- colorRampPalette(c("#f9e0e0", "#f28aaa","#9a133d","#3b0000"))(9)
color_gradient <-c("#27AAE1",color_gradient)


color_gradient <- colorRampPalette(c("#D3B843", "#788f33", "#165d43","#002d1b"))(9)
color_gradient <-c("#27AAE1",color_gradient)

library(ggplot2)
library(scales)  # for pretty_breaks, if needed

ecdf_plot <- E117RvsWT %>%
  ggplot(aes(log2FoldChange, color = ExonGroup)) +
  stat_ecdf(geom = "line", linewidth = 0.2) +
  coord_cartesian(xlim = c(-2, 2), ylim = c(0, 1)) +
  scale_color_manual(values = color_gradient) +
  scale_x_continuous(breaks = c(-2, 0, 2)) +
  scale_y_continuous(breaks = c(0, 0.5, 1)) +
  theme_linedraw(base_size = 5) +  # sets all text elements to 5pt
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    aspect.ratio = 1,
    legend.position = c(0.99, 0.01),
    legend.justification = c(1, 0),
    legend.title = element_text(size = 5),
    legend.text = element_text(size = 5),
    legend.key.size = unit(1.5, 'mm'),
    legend.box.background = element_blank(),
    legend.background = element_blank(), 
    axis.title = element_text(size = 5),
    axis.text = element_text(size = 5),
    plot.background = element_blank(),
    panel.background = element_blank(),
    plot.margin = unit(c(0, 0, 0, 0), "cm")
  ) +
  labs(y = "Cumulative frequency", x = "log2FC E117R/WT")
ggsave("ecdf_E117RvsWT.pdf", ecdf_plot,
       width = 20, height = 20, units = "mm", useDingbats = FALSE)
 



