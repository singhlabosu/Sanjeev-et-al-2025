###Figure1 plots

library(tidyverse)
library(ggrepel)
library(ggpubr)
library(patchwork)
library(rstatix)
setwd("~/OneDrive - The Ohio State University/Singhlab/Bioinfo/PYMpaperFiguresUpdate/Figure1/")

###reading files in ####
E117RvsWT<-read_csv("../DESeq2 for PYM samples/all_NT_res.csv")
E117RvsWT<-na.omit(E117RvsWT)
CHXvsWT<-read_csv("../DESeq2 for PYM samples/all_WT_res.csv")
CHXvsWT<-CHXvsWT%>%mutate(log2FoldChange = -1*log2FoldChange)#michaels email from 2/2 explains that the deseq was performed untreated/treated 
CHXvsWT<-na.omit(CHXvsWT)
allcolors<-c("Up" = "orange3", "Down" = "purple3", "Unchanged" = "grey60","Up" = "red4", "Down" = "blue4")
n=1
for (id in c("CHX","E117R")) {
  print(id)
  data <- get(paste0(id, "vsWT"))
  data<-data%>%mutate(
    Comparison = id,
    significant = case_when( padj<0.05 & log2FoldChange>0 ~ "Up",
                             padj<0.05 & log2FoldChange<0 ~ "Down",
                             TRUE ~ "Unchanged"),
    alpha_val = case_when(significant == "Down" ~ 0.9,
                          significant == "Up" ~ 0.9,
                          significant == "Unchanged"  ~ 0.3),
    log10padj=log10(padj))
  count_data <- data %>%
    filter(significant %in% c("Up", "Down")) %>%
    mutate(side = ifelse(log2FoldChange < 0, "left", "right")) %>%
    count(significant, side) %>%
    mutate(label = paste0("n=", n))
  if (id == "E117R") {  
    plotcolors <- allcolors[1:3]}
  else{plotcolors <- allcolors[3:5]}
  #plot:
  assign(paste0("f",n) ,ggplot(data, aes(x = log2FoldChange, y = -log10padj)) +
    geom_point(aes(color = significant, alpha = alpha_val), size = 1.5) +
    scale_color_manual(values = plotcolors) +
    scale_alpha_identity() +
    coord_cartesian(xlim = c(-6, 6), ylim = c(0, 20)) +
    geom_text(data = count_data, aes(
      x = ifelse(side == "left", -6, 6),
      y = 1,
      label = label,
      color = significant),
      inherit.aes = FALSE,
      hjust = ifelse(count_data$side == "left", 0, 1),
      size = 2)+
    labs(x = paste0("log2FC ",id,"/WT")) + 
    theme_linedraw()+
    theme(legend.position = "none",
          panel.background = element_blank(),   # No background color in the panel
          plot.background = element_blank(),
          panel.grid = element_blank(),
          panel.grid.major = element_line(linetype = "blank"),
          panel.grid.minor = element_line(linetype = "blank"),
          axis.title = element_text(size = 7),
          axis.text = element_text(size = 5),
          aspect.ratio = 6/14)
  )
  n=n+1
}

###Reading TE, NucCytoRanks and Halflives####

TE<-read.table("translation_efficiencies.txt")%>%
  mutate(LogTE = log10(V2))%>%
  select(1,3)
NucCyto<-read_tsv("nuclear fractions and percentiles.txt")%>%
  mutate(Cpercentile = 100 - Percentile)%>%
  select(1,4)
HalfLives<-read_tsv("Agarwal_RNA_halflives.txt",skip = 1)[,c(1,3)]

##Making 1st column Geneid for all dataframes
for (i in c("CHXvsWT", "E117RvsWT", "NucCyto", "TE", "HalfLives")) {
  df <- get(i)               # Access the data frame
  names(df)[1] <- "Geneid"   # Rename the first column
  assign(i, df)              # Reassign
}

allcolors<-c("orange3","grey60","purple3", "blue4","grey60","red4")

testresultsdf<-data_frame()
n=3
for (id in c("CHX","E117R")){
  for (param in c("NucCyto","TE","HalfLives")){
    print(paste0(id,"_",param))
    Yvar = case_when( param == "NucCyto" ~ "Cytoplasmic rank percentile",
                      param == "TE" ~ "log10 Normalized Ribosome Occupancy",
                      param == "HalfLives" ~ "log half-life (Aggregate)")
    df1 <- get(paste0(id, "vsWT"))
    df2 <- get(param)
    merged_df <- inner_join(df1, df2, by = "Geneid")
    names(merged_df)[8] <- Yvar
    merged_df <- merged_df %>%
      filter(is.finite(!!sym(Yvar)))
    
    if (id == "E117R") {
      Q1<-quantile(merged_df$log2FoldChange,probs = 0.25)
      Q2<-quantile(merged_df$log2FoldChange,probs = 0.75)
      merged_df<-merged_df %>%
        mutate(Xvar = case_when(#naming the column Xvar instead of E117RvsWT to keep the naming same for CHXvsWT
          log2FoldChange > paste0(Q2) ~ "Up",#"Up in \nE117R",
          log2FoldChange < paste0(Q1) ~ "Down", #"Down in \nE117R",
          TRUE ~ "No \nChange",
        ))
      plotcolors <- allcolors[1:3]}
    else{
      merged_df<-merged_df %>%
        mutate(Xvar = case_when(#naming the column Xvar instead of CHXvsWT to keep the naming same
          log2FoldChange > 0 & padj < 0.05  ~ "Up",#"Up in \nCHX",
          log2FoldChange < 0 & padj < 0.05 ~ "Down", #"Down in \nCHX",
          TRUE ~ "No \nChange",
        ))
      plotcolors <- allcolors[4:6]}
    
    medians <- merged_df %>% # Calculate medians for each group
      group_by(Xvar) %>%
      dplyr::summarize(median_value = median(!!sym(Yvar), na.rm = TRUE),
                       n = n())%>%
      #mutate(count = paste0("\nn=", n),
      mutate(count = paste0("\n(", n,")"),       
             y_label_pos = case_when( param == "NucCyto" ~ -4,
                                      param == "TE" ~ -5,
                                      param == "HalfLives" ~ -20))
    
    pair_test_df <- merged_df %>%
      wilcox_test(as.formula(paste0("`",Yvar, "` ~ Xvar")), p.adjust.method = "bonferroni") %>%
      mutate(
        y_position = case_when( param == "NucCyto" ~ seq(102, 102 + 6.6*2, by = 6.6),
                                param == "TE" ~ seq(0.2, 0.2 + 0.33*2, by = 0.33),
                                param == "HalfLives" ~ seq(17, 17 + 2.2*2, by = 2.2)),
        xmin       = group1,
        xmax       = group2,
        annotation = sprintf("%.2e", p.adj) )  
    
    testresultsdf<-bind_rows(testresultsdf,pair_test_df)
    
    assign(paste0("f",n) , ggplot(merged_df,aes(x=Xvar,y=!!sym(Yvar)))+
      geom_boxplot(aes(fill=Xvar) ,varwidth = TRUE, alpha=0.8,
                   outlier.size = 0.5,       # Smaller outlier dots
                   outlier.stroke = 0.2,     # Thinner outline for outliers
                   size = 0.3)+
      theme_classic()+
      theme(panel.border = element_blank(),       # No box around the panel
            panel.background = element_blank(),   # No background color in the panel
            plot.background = element_blank(),
            legend.position="none",aspect.ratio = 3/1,
            axis.line = element_line(size = 0.3),
            #axis.text.x = element_text(vjust = 0.5,angle = 65),
            axis.title = element_text(size = 7),
            axis.text = element_text(size = 5),
            axis.title.x=element_blank()) +
      scale_fill_manual(values=plotcolors)+
      geom_signif(
          data = pair_test_df,
          mapping = aes(xmin = xmin, xmax = xmax, y_position = y_position, annotations = annotation),
          manual = TRUE,
          inherit.aes = FALSE,
          tip_length = 0.005,
          size = 0.15,
          textsize = 1.8) +
      geom_text(data = medians, aes(x = Xvar, y = y_label_pos, 
                                    label = sprintf("%.2f", median_value), 
                                    color = Xvar), 
                inherit.aes = FALSE, size = 1.8, vjust = -0.2) +
      geom_text(data = medians, aes(x = Xvar, y = y_label_pos, 
                                    label = sprintf("%s", count)), 
                  inherit.aes = FALSE, size = 1.8, vjust = 0.5, color = "black") +
      scale_color_manual(values = plotcolors)
      )
    n=n+1
  }}

c<-f1|f2
d<-f3+f4+f5+f6+f7+f8+plot_layout(nrow=1)

layout <- c(
  area(3,1,4,6),area(5,1,8,6))

ggsave("F1patch1_final2.pdf",
       c/d+plot_layout(design = layout)+
         plot_annotation(tag_levels = "a")  &
         theme(plot.tag = element_text(face = 'bold'),
               panel.background = element_blank(),   
               plot.background = element_blank(),
               plot.margin = unit(c(0, 0, 0, 0), "cm")),
       width = 180, 
       height = 180, 
       units = "mm")

testresultsdf <- testresultsdf %>%
  mutate(across(where(is.character), ~ str_replace_all(., "\n", " ")))
write_tsv(testresultsdf,"StatResultsFig1.tsv")
