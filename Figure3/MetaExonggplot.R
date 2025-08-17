library(ggplot2)
library(dplyr)
library(readr)
library(tidyr)
library(tidyverse)
setwd("~/OneDrive - The Ohio State University/Singhlab/Bioinfo/PYMpaperFigures/")

# Define colors for each sample
colors_base <- c("red", "red3", "blue", "orange2", "orange", "#026690", 
                 "purple", "#56862e", "#909090", "#aaaaaa", "#c3c3c3") 

samples <- list(c("WT-CHX", "WT-NT", "E117R-CHX", "E117R-NT"))
names(samples) <- "CombinedSamples"
tests <- list(c("_CSD_exon_NL100nt", 100))

mod3 <- "_exon_E100_MANE_and_lncRNA_correct_min100_ann3p_dist_bin.dat"
mod5 <- "_exon_E100_MANE_and_lncRNA_correct_min100_ann5p_dist_bin.dat"

out3 <- "_3pDistributionNew.pdf"
out5 <- "_5pDistributionNew.pdf"

xlabel <- "Distance from End of Exon (nt)"
ylabel <- "Total Read Depth (Gene Normalized) Normalized to 1"
title <- "End of Exon Read Distribution"

working <- "~/OneDrive - The Ohio State University/Singhlab/Bioinfo/PYMpaperFigures/Figure2/UpdatedMetaExon3/"

# Function to load and process data
load_data <- function(file_path, sample_name, region) {
  temp_data <- read.table(file_path, header = FALSE, sep = "\t")
  norm_data <- temp_data[,1] / sum(temp_data[,1])  # Normalize
  data.frame(Position = seq_along(norm_data), ReadDepth = norm_data, Sample = sample_name)
}

# Process and plot 3' and 5' distributions
for (test in tests) {
  region <- as.integer(test[2])
  
  for (sample_group in names(samples)) {
    all_data_3p <- data.frame()
    all_data_5p <- data.frame()
    
    for (sample in samples[[sample_group]]) {
      file_3p <- paste0(working, sample, mod3)
      file_5p <- paste0(working, sample, mod5)
      
      data_3p <- load_data(file_3p, sample, region)
      data_5p <- load_data(file_5p, sample, region)
      
      all_data_3p <- bind_rows(all_data_3p, data_3p)
      all_data_5p <- bind_rows(all_data_5p, data_5p)
    }
    
    # Assign colors
    all_data_3p$Color <- factor(all_data_3p$Sample, levels = samples[[sample_group]], labels = colors_base[1:length(samples[[sample_group]])])
    all_data_5p$Color <- factor(all_data_5p$Sample, levels = samples[[sample_group]], labels = colors_base[1:length(samples[[sample_group]])])
    
    # Plot 3' distribution
    p3 <- ggplot(all_data_3p, aes(x = -Position, y = ReadDepth, color = Sample)) +
      geom_line(size = 1) +
      scale_color_manual(values = colors_base[1:length(samples[[sample_group]])]) +
      labs(title = title, x = xlabel, y = ylabel) +
      theme_minimal() +
      theme(panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            legend.position = "top")
    
    ggsave(filename = paste0( working, sample_group, out3), plot = p3, width = 11, height = 8.5, device = "pdf")
    
    # Plot 5' distribution
    p5 <- ggplot(all_data_5p, aes(x = Position, y = ReadDepth, color = Sample)) +
      geom_line(size = 1) +
      scale_color_manual(values = colors_base[1:length(samples[[sample_group]])]) +
      labs(title = title, x = xlabel, y = ylabel) +
      theme_minimal() +
      theme(panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            legend.position = "top")
    
    ggsave(filename = paste0(working, sample_group, out5), plot = p5, width = 11, height = 8.5, device = "pdf")
  }
}

ggplot(all_data_5p, aes(x = Position, y = ReadDepth, color = Sample)) +
  geom_line(size = 1) +
  coord_cartesian(xlim = c(0,100), ylim = c(0,0.01))+
  scale_color_manual(values = colors_base[1:length(samples[[sample_group]])]) +
  labs(title = title, x = xlabel, y = ylabel) +
  theme_light() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "top")

ggplot(all_data_3p, aes(x = -Position, y = ReadDepth, color = Sample)) +
  geom_line(size = 1) +
  scale_color_manual(values = colors_base[1:length(samples[[sample_group]])]) +
  coord_cartesian(xlim = c(-100,0), ylim = c(0,0.01))+
  labs(title = title, x = xlabel, y = ylabel) +
  theme_light() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "top")
  



