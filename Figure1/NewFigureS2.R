
setwd("~/OneDrive - The Ohio State University/Singhlab/Bioinfo/PYMpaperFiguresUpdate/Figure1/")

#Comparing overlaps of gene classifications
library(tidyverse)

###Making replicate correlation figure based on michael's new counts####

WT_NT1<-read_tsv("counts for PYM samples/WT_all_NT_1", skip = 1)
WT_NT2<-read_tsv("counts for PYM samples/WT_all_NT_2", skip = 1)
WT_CHX1<-read_tsv("counts for PYM samples/WT_all_CHX_1", skip = 1)
WT_CHX2<-read_tsv("counts for PYM samples/WT_all_CHX_2", skip = 1)
E117R_NT1<-read_tsv("counts for PYM samples/mutant_all_NT_1", skip = 1)
E117R_NT2<-read_tsv("counts for PYM samples/mutant_all_NT_2", skip = 1)
E117R_CHX1<-read_tsv("counts for PYM samples/mutant_all_CHX_1", skip = 1)
E117R_CHX2<-read_tsv("counts for PYM samples/mutant_all_CHX_2", skip = 1)


list<-c("WT_NT1","WT_NT2","WT_CHX1","WT_CHX2",
        "E117R_NT1","E117R_NT2","E117R_CHX1","E117R_CHX2")
counts_combined<-get(list[1])[,c(1,7)]
for (i in list[-1]) {
  print(i)
  df<-get(i)
  counts_combined<-cbind(counts_combined, df[,7])}

rownames(counts_combined)<-counts_combined$Geneid
counts_combined<-counts_combined[,-1]
str_extract(names(counts_combined),"MAGOH-[:graph:]+[1,2]")
list
names(counts_combined)<-list
names(counts_combined)
library(ggpointdensity)
library(scales)
counts_combined<-counts_combined[,1:8]
b<-counts_combined%>%ggplot(aes(WT_NT1,WT_NT2))+
  geom_pointdensity(shape = 1, size=0.5)+
  scale_colour_gradient(low = "grey80", high = "black", na.value = NA)+
  scale_x_log10(labels = scales::comma)+scale_y_log10(labels = scales::comma)+theme_linedraw()+
  theme(aspect.ratio = 1,panel.grid.major = element_line(linetype = "blank"),
        panel.grid.minor = element_line(linetype = "blank"),legend.position = "none",
        axis.text = element_text(size = 5),axis.title = element_text(size = 7))

c<-counts_combined%>%ggplot(aes(WT_CHX1,WT_CHX2))+
  geom_pointdensity(shape = 1, size=0.5)+
  scale_colour_gradient(low = "grey80", high = "black", na.value = NA)+
  scale_x_log10(labels = scales::comma)+scale_y_log10(labels = scales::comma)+theme_linedraw()+
  theme(aspect.ratio = 1,panel.grid.major = element_line(linetype = "blank"),
        panel.grid.minor = element_line(linetype = "blank"),legend.position = "none",
        axis.text = element_text(size = 5),axis.title = element_text(size = 7))


d<-counts_combined%>%ggplot(aes(E117R_NT1,E117R_NT2))+
  geom_pointdensity(shape = 1, size=0.5)+
  scale_colour_gradient(low = "grey80", high = "black", na.value = NA)+
  scale_x_log10(labels = scales::comma)+scale_y_log10(labels = scales::comma)+theme_linedraw()+
  theme(aspect.ratio = 1,panel.grid.major = element_line(linetype = "blank"),
        panel.grid.minor = element_line(linetype = "blank"),legend.position = "none",
        axis.text = element_text(size = 5),axis.title = element_text(size = 7))

e<-counts_combined%>%ggplot(aes(E117R_CHX1,E117R_CHX2))+
  geom_pointdensity(shape = 1, size=0.5)+
  scale_colour_gradient(low = "grey80", high = "black", na.value = NA)+
  scale_x_log10(labels = scales::comma)+scale_y_log10(labels = scales::comma)+theme_linedraw()+
  theme(aspect.ratio = 1,panel.grid.major = element_line(linetype = "blank"),
        panel.grid.minor = element_line(linetype = "blank"),legend.position = "none",
        axis.text = element_text(size = 5),axis.title = element_text(size = 7))

a+b+c+d


###DESeq2 for PCA plot####
#Creating a metadata table with information on each sample
sampleinfo<-as.data.frame(colnames(counts_combined))
sampleinfo$treatment<-as.factor(c("WT","WT","WT_CHX","WT_CHX",
                                  "E117R","E117R","E117R_CHX","E117R_CHX"))
names(sampleinfo)[1]<-"Sample"
sampleinfo

library(DESeq2)
ddscounts.table<-DESeqDataSetFromMatrix(counts_combined, sampleinfo, formula(~ treatment))
ddscounts <- DESeq(ddscounts.table)

#PCAplots
library(ggplot2)
library(ggrepel)
vsd <- vst(ddscounts, blind = FALSE)
a<-plotPCA(vsd,ntop = 3000, intgroup = "treatment")+
  #labs(title = "PCA plot of RIPiT samples")+
  theme_light()+
  theme(axis.text = element_text(size = 5),axis.title = element_text(size = 7),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        aspect.ratio = 1,
        legend.position = c(0.98, 0.3),         # Move legend to bottom right
        legend.justification = c(1, 0),     # Ensure correct alignment at bottom right
        legend.title=element_text(size=6), 
        legend.text=element_text(size=4),
        legend.key=element_blank(),
        legend.box.background = element_blank(),
        legend.background = element_rect(linetype = 1, size = 0.25, colour = 1),
        legend.key.size =unit(2, 'mm'))



##Plotting fold changes: density and scatter####

E117RvsWT<-read_csv("../DESeq2 for PYM samples/all_NT_res.csv")
E117RvsWT<-na.omit(E117RvsWT)
CHXvsWT<-read_csv("../DESeq2 for PYM samples/all_WT_res.csv")
CHXvsWT<-CHXvsWT%>%mutate(log2FoldChange = -1*log2FoldChange)#michaels email from 2/2 explains that the deseq was performed untreated/treated 
CHXvsWT<-na.omit(CHXvsWT)

f<-CHXvsWT%>%ggplot(aes(log2FoldChange))+
  geom_density(fill = "cyan", alpha= 0.4)+
  geom_density(data = E117RvsWT, fill = "orange", alpha= 0.4)+
  theme_linedraw()+
  theme(aspect.ratio = 1)+
  labs(#title = "Distribution of Fold changes",
       x="log2FC")+ 
  theme(panel.grid.major = element_line(linetype = "blank"),
        panel.grid.minor = element_line(linetype = "blank"),
        axis.title = element_text(size = 6),
        axis.text = element_text(size = 4))

library(patchwork)
layout <- c(
  area(1,1, 2,2),area(1,3,2,4),area(1,5,2,6),
  area(3,1,4,2),area(3,3,4,4),area(3,5,4,6),
  area(5,1,6,2),area(5,3,6,4),
  area(7,5,10,8))
plot(layout)

ggsave("FigureS2patch2.pdf",
       a+b+c+d+e+f+plot_layout(design = layout)+
         plot_annotation(title = 'FigureS2', tag_levels = "a")  &
         theme(plot.tag = element_text(face = 'bold')),
       width = 8.5, 
       height = 11, 
       units = "in")


setwd("~/OneDrive - The Ohio State University/Singhlab/Westerns/2025/")
library(tidyverse)
rawquant<-read_csv("FractionationAfterCHXQuant.csv")


normvalues<-tibble()
for(i in unique(rawquant$Replicate)){
  print(i)
  for(j in unique(rawquant$Protein)){
    print(j)
    for(k in unique(rawquant$Treatment)){
      filtereddf<-rawquant%>%filter(Replicate==i & Protein==j & Treatment == k)%>%mutate(
        NormLevels = Signal/sum(Signal))
      normvalues<-rbind(normvalues,filtereddf)}}}
localizations <- c("Cytoplasm", "Membrane", "Chromatin")
normvalues$Localization <- rep(localizations, length.out = nrow(normvalues))
normvalues$Treatment <- factor(normvalues$Treatment, levels = c("Control", "CHX"))

normvalues %>%
  ggplot(aes(Treatment, NormLevels, fill = Localization)) +
  stat_summary(fun = mean, geom = "bar", position = position_dodge(width = 0.6), width = 0.6) +
  stat_summary(fun.data = mean_cl_boot, geom = "errorbar", width = 0.3, position = position_dodge(width = 0.6)) +
  geom_point(aes(colour = Replicate, group = Localization), 
             shape = 1, 
             position = position_dodge(width = 0.6)) +
  scale_fill_brewer(palette = "Set2") +
  scale_color_brewer(palette = "Dark2")+
  theme_linedraw() +theme(aspect.ratio = 2,panel.grid.major = element_line(linetype = "blank"),
                          panel.grid.minor = element_line(linetype = "blank"))+
  labs(title = "EJC proteins localization",
       y ="Normalized Protein Levels" ) + 
  facet_wrap(~Protein)



rawquantIP<-read_csv("IPafterCHXQuant.csv")

rawquantIP<-rawquantIP%>%filter(Protein != "RBM8A")
normvalues<-tibble()
for(i in unique(rawquantIP$Replicate)){
  print(i)
  for(j in unique(rawquantIP$Treatment)){
    print(j)
    for(k in unique(rawquantIP$Sample)){
      filtereddf<-rawquantIP%>%filter(Replicate==i & Treatment==j & Sample == k)
      ref_value <- filtereddf$Signal[1]
      filtereddf<-filtereddf%>%mutate(
        NormLevels = Signal/ref_value)
      normvalues<-rbind(normvalues,filtereddf)}}}
normvalues$Treatment <- factor(normvalues$Treatment, levels = c("Control", "CHX"))
library(ggpubr)
normvalues<-normvalues %>% mutate(NormLevels = 1 / NormLevels)%>%
  filter(Protein == "MAGOH"& Sample=="IP")

a<-normvalues %>%
  ggplot(aes(Protein, NormLevels, fill = Treatment)) + 
  stat_summary(fun = mean, geom = "bar", position = position_dodge(width = 0.6), width = 0.6) +
  stat_summary(fun.data = mean_cl_boot, geom = "errorbar", width = 0.3, position = position_dodge(width = 0.6)) +
  geom_point(aes(colour = Replicate, group = Treatment), 
             shape = 1, 
             position = position_dodge(width = 0.6)) +
  stat_compare_means(
    aes(group = Treatment),
    method = "wilcox.test", 
    label = "p.format",
    label.y = 18, size = 3,tip.length = 0.005
  ) +
  scale_fill_brewer(palette = "Set2") +
  scale_color_brewer(palette = "Dark2") +
  theme_linedraw() +theme(aspect.ratio = 3)+
  labs(title = "IP : FLAG-PYM1", y = "Normalized Protein Levels") + 
  theme(panel.grid.major = element_line(linetype = "blank"),
        panel.grid.minor = element_line(linetype = "blank"))
a

compare_means(NormLevels ~ Treatment, data = normvalues, method = "wilcox.test")
wilcox.test(NormLevels ~ Treatment, data = normvalues)

rawquantIP<-read_csv("IPafterCHXQuant.csv")

rawquantIP<-rawquantIP%>%filter(Protein != "RBM8A")
rawquantIP<-rawquantIP%>%filter(Sample=="TotalExtract")
normvalues<-tibble()
for(i in unique(rawquantIP$Replicate)){
  print(i)
  for(j in unique(rawquantIP$Treatment)){
    print(j)
    filtereddf<-rawquantIP%>%filter(Replicate==i & Treatment==j)
    ref_value <- filtereddf %>% filter(Protein == "TUBULIN") %>% pull(Signal)
    filtereddf<-filtereddf%>%mutate(
      NormLevels = Signal/ref_value)
    normvalues<-rbind(normvalues,filtereddf)}}
normvalues$Treatment <- factor(normvalues$Treatment, levels = c("Control", "CHX"))
library(ggpubr)



b<-normvalues %>%
  filter(Protein !="TUBULIN") %>%mutate(NormLevels=NormLevels*10)%>%
  ggplot(aes(Protein, NormLevels, fill = Treatment)) + 
  stat_summary(fun = mean, geom = "bar", position = position_dodge(width = 0.6), width = 0.6) +
  stat_summary(fun.data = mean_cl_boot, geom = "errorbar", width = 0.3, position = position_dodge(width = 0.6)) +
  geom_point(aes(colour = Replicate, group = Treatment), 
             shape = 1, 
             position = position_dodge(width = 0.6)) +
  stat_compare_means(
    aes(group = Treatment),
    method = "wilcox.test", paired = TRUE,
    label = "p.format",
    label.y = 2.5, size = 3,tip.length = 0.005
  ) +
  scale_fill_brewer(palette = "Set2") +
  scale_color_brewer(palette = "Dark2") +
  theme_linedraw() +theme(aspect.ratio = 2,panel.grid.major = element_line(linetype = "blank"),
                          panel.grid.minor = element_line(linetype = "blank"))+scale_y_log10()+
  labs(title = "Total Extract", y = "Normalized Protein Levels")
b

compare_means(NormLevels ~ Treatment, group.by = "Protein", data = normvalues, method = "wilcox.test")
wilcox.test(NormLevels ~ Treatment, data = normvalues,subset = Protein %in% c("MAGOH"))
wilcox.test(NormLevels ~ Treatment, data = normvalues,subset = Protein %in% c("PYM1"))

library(patchwork)
(b + a) + plot_layout(guides = "collect")



####GOTERM ANALYSIS####
#CHXvsWT
setwd("~/Library/CloudStorage/OneDrive-TheOhioStateUniversity/Singhlab/Bioinfo/PYMpaperFiguresUpdate/Figure1/GOTermAnalysis")
# CHXvsWT%>%filter(log2FoldChange > 0, padj < 0.05)%>%select(V1)
# write_tsv(CHXvsWT%>%filter(log2FoldChange > 0, padj < 0.05)%>%select(V1),"CHXup.txt", col_names = F)
# write_tsv(CHXvsWT%>%filter(log2FoldChange < 0, padj < 0.05)%>%select(V1),"CHXdown.txt", col_names = F)
# write_tsv(CHXvsWT%>%select(V1),"CHXgenes.txt", col_names = F)
# 
# E117RvsWT
# summary(E117RvsWT$log2FoldChange)
# write_tsv(E117RvsWT%>%filter(log2FoldChange > 0.26395)%>%select(V1),"E117Rup.txt", col_names = F)
# write_tsv(E117RvsWT%>%filter(log2FoldChange < -0.270291)%>%select(V1),"E117Rdown.txt", col_names = F)
# write_tsv(E117RvsWT%>%select(V1),"E117Rgenes.txt", col_names = F)
# 
# #plot
# CHXdownGOresult<-read_tsv("CHXdownResults.txt")
# CHXupGOresult<-read_tsv("CHXupResults.txt")
# E117RdownGOresult<-read_tsv("E117RdownResults.txt")
# E117RupGOresult<-read_tsv("E117RupResults.txt")
# 
# CHXdownGOresult %>%
#   group_by(Category) %>%
#   arrange(Bonferroni) %>%
#   top_n(-5, Bonferroni) %>%  
#   ggplot(aes(x=Benjamini,y=Term, fill=Category))+
#   geom_col()+
#   scale_x_log10()
# 
# CHXupGOresult %>%
#   group_by(Category) %>%
#   arrange(Bonferroni) %>%
#   top_n(-5, Bonferroni) %>%  
#   ggplot(aes(x=Benjamini,y=Term, fill=Category))+
#   geom_col()+
#   scale_x_log10()
# 
# E117RupGOresult %>%
#   group_by(Category) %>%
#   arrange(Bonferroni) %>%
#   top_n(-5, Bonferroni) %>%  
#   ggplot(aes(x=Benjamini,y=Term, fill=Category))+
#   geom_col()+
#   scale_x_log10()
# 
# E117RdownGOresult %>%
#   group_by(Category) %>%
#   arrange(Bonferroni) %>%
#   top_n(-5, Bonferroni) %>%  
#   ggplot(aes(x=Benjamini,y=Term, fill=Category))+
#   geom_col()+
#   scale_x_log10()


CHXdownKW<-read_tsv("KW/CHXdownKW.txt")
CHXupKW<-read_tsv("KW/CHXupKW.txt")
E117RdownKW<-read_tsv("KW/E117RdownKW.txt")
E117RupKW<-read_tsv("KW/E117RupKW.txt")


CHXdownKW$Comparison<-"CHXdown"
CHXupKW$Comparison<-"CHXup"
E117RdownKW$Comparison<-"E117Rdown"
E117RupKW$Comparison<-"E117Rup"

CHXdownKW %>%
  arrange(desc(Bonferroni)) %>%
  top_n(-7, Bonferroni) %>%  
  ggplot(aes(x = Benjamini, y = Term, color = Category, size = Count)) +
  geom_col(size = 0.02, width = 0.02, fill = "grey")+
  geom_point(position = position_dodge(width = 0.5)) +  # Adjust dodge width as needed
  scale_x_log10()+
  theme_light() + theme(panel.grid.major = element_line(linetype = "blank"),
                        panel.grid.minor = element_line(linetype = "blank"))

CHXupKW %>%
  arrange(desc(Bonferroni)) %>%
  top_n(-7, Bonferroni) %>%  
  ggplot(aes(x = Benjamini, y = Term, color = Category, size = Count)) +
  geom_col(size = 0.02, width = 0.02, fill = "grey")+
  geom_point(position = position_dodge(width = 0.5)) +  # Adjust dodge width as needed
  scale_x_log10()+
  theme_light() + theme(panel.grid.major = element_line(linetype = "blank"),
                        panel.grid.minor = element_line(linetype = "blank"))


E117RdownKW %>%
  arrange(desc(Bonferroni)) %>%
  top_n(-7, Bonferroni) %>%  
  ggplot(aes(x = Benjamini, y = Term, color = Category, size = Count)) +
  geom_col(size = 0.02, width = 0.02, fill = "grey")+
  geom_point(position = position_dodge(width = 0.5)) +  # Adjust dodge width as needed
  scale_x_log10()+
  theme_light() + theme(panel.grid.major = element_line(linetype = "blank"),
                        panel.grid.minor = element_line(linetype = "blank"))

E117RupKW %>%
  arrange(desc(Bonferroni)) %>%
  top_n(-7, Bonferroni) %>%  
  ggplot(aes(x = Benjamini, y = Term, color = Category, size = Count)) +
  geom_col(size = 0.02, width = 0.02, fill = "grey")+
  geom_point(position = position_dodge(width = 0.5)) +  # Adjust dodge width as needed
  scale_x_log10()+
  theme_light() + theme(panel.grid.major = element_line(linetype = "blank"),
                        panel.grid.minor = element_line(linetype = "blank"))

CHXdownKW %>%
  arrange(Bonferroni) %>%
  top_n(-7, Bonferroni) %>% 
  mutate(Term = factor(Term, levels = rev(Term[order(Bonferroni)]))) %>%  # Reorder Term by Bonferroni
  ggplot(aes(x = Benjamini, y = Term, color = Category,size = Count)) +
  geom_col(size = 0.02, width = 0.02, fill = "grey")+
  geom_point(position = position_dodge(width = 0.5)) +  # Adjust dodge width as needed
  scale_x_log10()+
  theme_light() + theme(panel.grid.major = element_line(linetype = "blank"),
                        panel.grid.minor = element_line(linetype = "blank"))

KWcombined<-rbind(CHXdownKW,CHXupKW,E117RdownKW,E117RupKW)
KWcombined$Term<-str_split(KWcombined$Term, "~", simplify = T)[,2]
KWcombined$Category <- str_split(KWcombined$Category, "KW_", simplify = T)[,2]

a<-KWcombined %>% filter(Comparison == "CHXdown")%>%
  arrange(desc(Bonferroni)) %>%
  top_n(-7, Bonferroni) %>%  
  mutate(Term = factor(Term, levels = rev(Term[order(Bonferroni)]))) %>%  # Reorder Term by Bonferroni
  ggplot(aes(x = Benjamini, y = Term, color = Category, size = Count)) +
  geom_col(size = 0.02, width = 0.02, fill = "grey", color = NA)+
  geom_point(position = position_dodge(width = 0.5)) +  # Adjust dodge width as needed
  scale_x_log10()+labs(title = "CHXdown")+
  theme_light() + theme(panel.grid.major = element_line(linetype = "blank"),
                        panel.grid.minor = element_line(linetype = "blank"))


b<-KWcombined %>% filter(Comparison == "CHXup")%>%
  arrange(desc(Bonferroni)) %>%
  top_n(-7, Bonferroni) %>%  
  mutate(Term = factor(Term, levels = rev(Term[order(Bonferroni)]))) %>%  # Reorder Term by Bonferroni
  ggplot(aes(x = Benjamini, y = Term, color = Category, size = Count)) +
  geom_col(size = 0.02, width = 0.02, fill = "grey", color = NA)+
  geom_point(position = position_dodge(width = 0.5)) +  # Adjust dodge width as needed
  scale_x_log10()+labs(title = "CHXup")+
  theme_light() + theme(panel.grid.major = element_line(linetype = "blank"),
                        panel.grid.minor = element_line(linetype = "blank"))

c<-KWcombined %>% filter(Comparison == "E117Rdown")%>%
  arrange(desc(Bonferroni)) %>%
  top_n(-7, Bonferroni) %>%  
  mutate(Term = factor(Term, levels = rev(Term[order(Bonferroni)]))) %>%  # Reorder Term by Bonferroni
  ggplot(aes(x = Benjamini, y = Term, color = Category, size = Count)) +
  geom_col(size = 0.02, width = 0.02, fill = "grey", color = NA)+
  geom_point(position = position_dodge(width = 0.5)) +  # Adjust dodge width as needed
  scale_x_log10()+labs(title = "E117Rdown")+
  theme_light() + theme(panel.grid.major = element_line(linetype = "blank"),
                        panel.grid.minor = element_line(linetype = "blank"))


d<-KWcombined %>% filter(Comparison == "E117Rup")%>%
  arrange(desc(Bonferroni)) %>%
  top_n(-7, Bonferroni) %>%  
  mutate(Term = factor(Term, levels = rev(Term[order(Bonferroni)]))) %>%  # Reorder Term by Bonferroni
  ggplot(aes(x = Benjamini, y = Term, color = Category, size = Count)) +
  geom_col(size = 0.02, width = 0.02, fill = "grey", color = NA)+
  geom_point(position = position_dodge(width = 0.5)) +  # Adjust dodge width as needed
  scale_x_log10()+labs(title = "E117Rup")+
  theme_light() + theme(panel.grid.major = element_line(linetype = "blank"),
                        panel.grid.minor = element_line(linetype = "blank"))

a+b+c+d
table(KWcombined$Comparison)

filtered_data <- KWcombined %>%
  filter(Comparison %in% c("CHXdown", "CHXup", "E117Rdown", "E117Rup")) %>%
  group_by(Comparison) %>%
  distinct(Term, .keep_all = TRUE) %>%
  arrange(Benjamini) %>%
  slice_min(Benjamini, n = 7) %>%
  ungroup() %>%
  mutate(Term = factor(paste(Comparison, Term, sep = "_"), levels = rev(unique(paste(Comparison, Term, sep = "_")[order(Benjamini)]))))

# Plot with faceting
ggplot(filtered_data, aes(x = Benjamini, y = Term, color = Category, size = Count)) +
  geom_col(size = 0.02, width = 0.02, fill = "grey", color = NA) +
  geom_point() +
  scale_x_log10(n.breaks = 4) +
  scale_size_continuous(range = c(3, 10)) +  # Adjust point size
  scale_color_manual(values = c("BIOLOGICAL_PROCESS" = "#dd5129", 
                                "CELLULAR_COMPONENT" = "#43b284",
                                "MOLECULAR_FUNCTION" = "#fab255")) +  # Customize colors
  scale_y_discrete(labels = function(x) sub(".*_", "", x)) +  # Remove appended comparison from labels
  labs(title = "Comparison Plots", x = "Benjamini", y = "Term") +
  theme_light() +
  theme(aspect.ratio = 1.2,
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "bottom",   # Move legend to the bottom
        legend.box = "horizontal",    # Arrange legend horizontally
        legend.title = element_blank()  # Optional: remove legend title
  ) +
  facet_wrap(~Comparison, scales = "free", ncol = 2)  # Facet by Comparison



rawquant<-read_csv("FractionationAfterCHXQuant.csv")


normvalues<-tibble()
for(i in unique(rawquant$Replicate)){
  print(i)
  for(j in unique(rawquant$Protein)){
    print(j)
    for(k in unique(rawquant$Treatment)){
      filtereddf<-rawquant%>%filter(Replicate==i & Protein==j & Treatment == k)%>%mutate(
        NormLevels = Signal/sum(Signal))
      normvalues<-rbind(normvalues,filtereddf)}}}
localizations <- c("Cytoplasm", "Membrane", "Chromatin")
normvalues$Localization <- rep(localizations, length.out = nrow(normvalues))
normvalues$Treatment <- factor(normvalues$Treatment, levels = c("Control", "CHX"))

normvalues %>%
  ggplot(aes(Treatment, NormLevels, fill = Localization)) +
  stat_summary(fun = mean, geom = "bar", position = position_dodge(width = 0.6), width = 0.6) +
  stat_summary(fun.data = mean_cl_boot, geom = "errorbar", width = 0.3, position = position_dodge(width = 0.6)) +
  geom_point(aes(colour = Replicate, group = Localization), 
             shape = 1, 
             position = position_dodge(width = 0.6)) +
  scale_fill_brewer(palette = "Set2") +
  scale_color_brewer(palette = "Dark2")+
  theme_linedraw() +theme(aspect.ratio = 2,panel.grid.major = element_line(linetype = "blank"),
                          panel.grid.minor = element_line(linetype = "blank"))+
  labs(title = "EJC proteins localization",
       y ="Normalized Protein Levels" ) + 
  facet_wrap(~Protein)



rawquantIP<-read_csv("IPafterCHXQuant.csv")

rawquantIP<-rawquantIP%>%filter(Protein != "RBM8A")
normvalues<-tibble()
for(i in unique(rawquantIP$Replicate)){
  print(i)
  for(j in unique(rawquantIP$Treatment)){
    print(j)
    for(k in unique(rawquantIP$Sample)){
      filtereddf<-rawquantIP%>%filter(Replicate==i & Treatment==j & Sample == k)
      ref_value <- filtereddf$Signal[1]
      filtereddf<-filtereddf%>%mutate(
        NormLevels = Signal/ref_value)
      normvalues<-rbind(normvalues,filtereddf)}}}
normvalues$Treatment <- factor(normvalues$Treatment, levels = c("Control", "CHX"))
library(ggpubr)
normvalues<-normvalues %>% mutate(NormLevels = 1 / NormLevels)%>%
  filter(Protein == "MAGOH"& Sample=="IP")

a<-normvalues %>%
  ggplot(aes(Protein, NormLevels, fill = Treatment)) + 
  stat_summary(fun = mean, geom = "bar", position = position_dodge(width = 0.6), width = 0.6) +
  stat_summary(fun.data = mean_cl_boot, geom = "errorbar", width = 0.3, position = position_dodge(width = 0.6)) +
  geom_point(aes(colour = Replicate, group = Treatment), 
             shape = 1, 
             position = position_dodge(width = 0.6)) +
  stat_compare_means(
    aes(group = Treatment),
    method = "wilcox.test", 
    label = "p.format",
    label.y = 18, size = 3,tip.length = 0.005
  ) +
  scale_fill_brewer(palette = "Set2") +
  scale_color_brewer(palette = "Dark2") +
  theme_linedraw() +theme(aspect.ratio = 3)+
  labs(title = "IP : FLAG-PYM1", y = "Normalized Protein Levels") + 
  theme(panel.grid.major = element_line(linetype = "blank"),
        panel.grid.minor = element_line(linetype = "blank"))
a

compare_means(NormLevels ~ Treatment, data = normvalues, method = "wilcox.test")
wilcox.test(NormLevels ~ Treatment, data = normvalues)

rawquantIP<-read_csv("IPafterCHXQuant.csv")

rawquantIP<-rawquantIP%>%filter(Protein != "RBM8A")
rawquantIP<-rawquantIP%>%filter(Sample=="TotalExtract")
normvalues<-tibble()
for(i in unique(rawquantIP$Replicate)){
  print(i)
  for(j in unique(rawquantIP$Treatment)){
    print(j)
    filtereddf<-rawquantIP%>%filter(Replicate==i & Treatment==j)
    ref_value <- filtereddf %>% filter(Protein == "TUBULIN") %>% pull(Signal)
    filtereddf<-filtereddf%>%mutate(
      NormLevels = Signal/ref_value)
    normvalues<-rbind(normvalues,filtereddf)}}
normvalues$Treatment <- factor(normvalues$Treatment, levels = c("Control", "CHX"))
library(ggpubr)



b<-normvalues %>%
  filter(Protein !="TUBULIN") %>%mutate(NormLevels=NormLevels*10)%>%
  ggplot(aes(Protein, NormLevels, fill = Treatment)) + 
  stat_summary(fun = mean, geom = "bar", position = position_dodge(width = 0.6), width = 0.6) +
  stat_summary(fun.data = mean_cl_boot, geom = "errorbar", width = 0.3, position = position_dodge(width = 0.6)) +
  geom_point(aes(colour = Replicate, group = Treatment), 
             shape = 1, 
             position = position_dodge(width = 0.6)) +
  stat_compare_means(
    aes(group = Treatment),
    method = "t.test", paired = TRUE,
    label = "p.format",
    label.y = 2.5, size = 3,tip.length = 0.005
  ) +
  scale_fill_brewer(palette = "Set2") +
  scale_color_brewer(palette = "Dark2") +
  theme_linedraw() +theme(aspect.ratio = 2,panel.grid.major = element_line(linetype = "blank"),
                          panel.grid.minor = element_line(linetype = "blank"))+scale_y_log10()+
  labs(title = "Total Extract", y = "Normalized Protein Levels")
b

compare_means(NormLevels ~ Treatment, group.by = "Protein", data = normvalues, method = "wilcox.test")
wilcox.test(NormLevels ~ Treatment, data = normvalues,subset = Protein %in% c("MAGOH"))
wilcox.test(NormLevels ~ Treatment, data = normvalues,subset = Protein %in% c("PYM1"))

library(patchwork)
(b + a) + plot_layout(guides = "collect")


