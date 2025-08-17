setwd("~/OneDrive - The Ohio State University/Singhlab/Bioinfo/PYMpaperFiguresUpdate/Figure3/")
library(tidyverse)
#Exon. numbers plot####
###reading files in 
E117RvsWT<-read_csv("../DESeq2 for PYM samples/all_NT_res.csv")
E117RvsWT<-na.omit(E117RvsWT)
names(E117RvsWT)[1]<-"ensembl_gene_id"

###Getting exon numbers from ensembl
Exons<-read_tsv("../MANEexons.tsv")
SEgenes_MANE_noInt<-read_tsv("../Figure4/Kovalak_etal/SEgenes_MANE_noInt.tsv")

E117RvsWT<-inner_join(E117RvsWT, Exons[,c(1,4)])
E117RvsWT<-E117RvsWT%>%mutate(ExonGroup = case_when(rank > 1 ~ "Multi",
                                                    rank == 1~ "Single"))
E117RvsWT$ExonGroup
table(E117RvsWT$ExonGroup)
E117RvsWT <- E117RvsWT %>%
  filter(!(ExonGroup == "Single" & !ensembl_gene_id %in% SEgenes_MANE_noInt$ensembl_gene_id))


E117RvsWT<-E117RvsWT%>%mutate(ExonGroup = case_when(rank > 9 ~ "10+",
                                                    TRUE ~ as.character(rank)))
E117RvsWT$ExonGroup <- factor(E117RvsWT$ExonGroup, levels = c(as.character(1:9), "10+"))


color_gradient <- colorRampPalette(c("#D3B843", "#788f33", "#165d43","#002d1b"))(9)
color_gradient <-c("#27AAE1",color_gradient)

E117RvsWT %>%
  ggplot(aes(log2FoldChange, color = ExonGroup)) +
  stat_ecdf(geom = "line", linewidth = 0.8) +
  coord_cartesian(xlim = c(-2, 2)) +
  scale_color_manual(values = color_gradient)+
  theme_linedraw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        aspect.ratio = 1,
        legend.position = c(0.99, 0.01),
        legend.justification = c(1, 0),
        legend.title=element_text(size=9), 
        legend.text=element_text(size=7),
        legend.key.size =unit(4, 'mm'),
        legend.box.background = element_blank()) +
  labs(y = "Cumulative frequency", x = "log2FC E117R/WT")


##For fig S3

##SLCN plot for CHX treatment#####

# Reading in files
CHXvsNTslcn <- read_csv("../SLCN_counts/slcn-trimmed-last-allcans-final_WT_res.csv")
#CHXvsNTslcn <- read_csv("../Figure2/SLCN_counts/slcn-trimmed-last-allcans-final_NT_res.csv")
#plot#####
names(CHXvsNTslcn)[1] <- "ensembl_gene_id"
CHXvsNTslcn <- CHXvsNTslcn %>% filter(baseMean > 10)

# Extracting Gene ID and Region
CHXvsNTslcn$Geneid <- str_split(CHXvsNTslcn$ensembl_gene_id, "-", simplify = T)[, 1]
CHXvsNTslcn$Region <- str_split(CHXvsNTslcn$ensembl_gene_id, "-", simplify = T)[, 2]

# Check counts per region
table(CHXvsNTslcn$Region)
CHXvsNTslcn <- CHXvsNTslcn %>%
  filter(!(Region == "s" & !Geneid %in% SEgenes_MANE_noInt$ensembl_gene_id))

# CHXvsNTslcn <- CHXvsNTslcn %>%
#   filter(!(Region == "l" & !Geneid %in% NoInt3UTR$ensembl_gene_id))
# table(CHXvsNTslcn$Region)

# Pairwise Wilcoxon test (include all regions, including "s")
p_values <- pairwise.wilcox.test(
  CHXvsNTslcn$log2FoldChange, 
  CHXvsNTslcn$Region, 
  p.adjust.method = "bonferroni"
)$p.value
stat<-pairwise.wilcox.test(
  CHXvsNTslcn$log2FoldChange, 
  CHXvsNTslcn$Region, 
  p.adjust.method = "bonferroni")
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
CHXvsNTslcn <- CHXvsNTslcn %>%
  mutate(Regions = recode(Region,
                          "c" = "canonical",
                          "l" = "last exon",
                          "n" = "non-canonical",
                          "s" = "single exon"))

# Calculate the number of entries in each region
region_counts <- CHXvsNTslcn %>%
  group_by(Regions) %>%
  summarize(count = n(), .groups = 'drop')

# Create formatted labels with counts
labels_with_counts <- paste(region_counts$Regions, 
                            " (", region_counts$count, ")", sep = "")
region_colors <- c(
  c = "#A97C50",  
  l = "#F15A29",  
  n = "#FBB040",  
  s = "#27AAE1")  


# Plotting
s1<-CHXvsNTslcn %>%
  ggplot(aes(log2FoldChange, color = Region)) +
  stat_ecdf(geom = "line", linewidth = 0.8) +
  scale_color_manual(values = region_colors, 
                     labels = labels_with_counts) +  # Change legend text
  coord_cartesian(xlim = c(-2, 2)) +
  theme_linedraw() +
  theme(
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    aspect.ratio = 1,
    legend.position = c(0.99, 0.01),         # Move legend to bottom right
    legend.justification = c(1, 0)          # Ensure correct alignment at bottom right
  ) +
  labs(title = "CHX vs NT fold changes",
       subtitle = "Grouped by EJC regions",
       y = "Cumulative frequency", x = "Log2FC CHX/NT") +
  annotate("text", x = -1.9, y = 0.82, label = pval_text, hjust = 0, size = 3, color = "black")  # Add p-values on top left
s1


#Plotting exon numbers for CHXvsWT
CHXvsNT<-read_csv("../DESeq2 for PYM samples/all_WT_res.csv")
CHXvsNT<-na.omit(CHXvsNT)
names(CHXvsNT)[1]<-"ensembl_gene_id"

###Getting exon numbers 

CHXvsNT<-inner_join(CHXvsNT, Exons[,c(1,4)])

CHXvsNT<-CHXvsNT%>%mutate(ExonGroup = case_when(rank > 9 ~ "10+",
                                                TRUE ~ as.character(rank)))
CHXvsNT$ExonGroup <- factor(CHXvsNT$ExonGroup, levels = c(as.character(1:9), "10+"))

color_gradient <- colorRampPalette(c("#D3B843", "#788f33", "#165d43","#002d1b"))(10)
s2<-CHXvsNT%>%ggplot(aes(log2FoldChange, color = ExonGroup))+
  stat_ecdf(geom = "line", linewidth = 0.8)+
  coord_cartesian(xlim = c(-2,2))+
  scale_color_manual(values = color_gradient)+
  theme_linedraw()+
  theme(panel.grid.major = element_line(linetype = "blank"),
        panel.grid.minor = element_line(linetype = "blank"),
        aspect.ratio = 1)+
  labs(title = "E117R vs WT foldchanges by exon number",
       y = "Cumulative frequency", x = "log2FC CHX/NT")
s2
s1+s2


#####Testing p-values
# Prepare variables
data <- CHXvsNTslcn
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

# Adjust p-values (Bonferroni, same as in your test)
results_df <- results_df %>%
  mutate(p_adjusted = p.adjust(p_value, method = "bonferroni"),
         hypothesis = paste(group1, "vs", group2))

# Reorder columns
results_df <- results_df %>%
  select(hypothesis, group1, group2, W_statistic, p_value, p_adjusted, alternative, method)

print(results_df)


pval_df <- as.data.frame(as.table(stat$p.value)) %>%
  filter(!is.na(Freq)) %>%
  rename(group1 = Var1, group2 = Var2, p_adjusted = Freq)