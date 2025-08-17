##Comparing RNA levels of HeLa Cells and HEK293 cells

#Reading in HeLa transcript level information
#Two options: Gehring lab and Hentze lab
#CASC3 dataset
library(tidyverse)
for (i in 1:3){
  dfname<-paste0("Rep", i)
  filename<-paste0("~/OneDrive - The Ohio State University/Singhlab/Bioinfo/CASC3_RNPS1_SMG67_comparison/Kallisto_CASC3/WT",i,"/abundance.tsv")
  assign(dfname, read_tsv(filename)%>% dplyr::select("target_id", "tpm"))
}

for (i in 1:3){
  dfname<-paste0("Rep", i)
  filename<-paste0("~/Downloads/HauerHelaKallisto/ERR120144",i+6,"/abundance.tsv")
  assign(dfname, read_tsv(filename)%>% dplyr::select("target_id", "tpm"))
}

#Rename tpm columns before merging to avoid conflicts
Rep1 <- Rep1 %>% dplyr::rename(tpm1 = tpm)
Rep2 <- Rep2 %>% dplyr::rename(tpm2 = tpm)
Rep3 <- Rep3 %>% dplyr::rename(tpm3 = tpm)

#Join all three dataframes by 'target_id'
Hela <- Rep1 %>%
  inner_join(Rep2, by = "target_id") %>%
  inner_join(Rep3, by = "target_id")

#Compute the average TPM
Hela <- Hela %>%
  mutate(avg_tpm_Hela = (tpm1 + tpm2 + tpm3) / 3)

#Reading in HEK transcript level information
for (i in 1:3){
  dfname<-paste0("Rep", i)
  path<-paste0("~/OneDrive - The Ohio State University/Singhlab/Bioinfo/PYMpaperFiguresUpdate/Figure5/Kallisto/MS22_0",i, "*", "/abundance.tsv")
  filename<-Sys.glob(path)
  assign(dfname, read_tsv(filename)%>% dplyr::select("target_id", "tpm"))
}


#Rename tpm columns before merging to avoid conflicts
Rep1 <- Rep1 %>% dplyr::rename(tpm1 = tpm)
Rep2 <- Rep2 %>% dplyr::rename(tpm2 = tpm)
Rep3 <- Rep3 %>% dplyr::rename(tpm3 = tpm)

#Join all three dataframes by 'target_id'
HEK <- Rep1 %>%
  inner_join(Rep2, by = "target_id") %>%
  inner_join(Rep3, by = "target_id")

#Compute the average TPM
HEK <- HEK %>%
  mutate(avg_tpm_HEK = (tpm1 + tpm2 + tpm3) / 3)

length(intersect(Hela$target_id,HEK$target_id))
##Combine Hela and HEK
CombinedTPMs<-inner_join(HEK[,c(1,5)],Hela[,c(1,5)])

CombinedTPMs <- CombinedTPMs %>%
  mutate(target_id=str_split(target_id,"\\.",simplify = T)[,1])
CombinedTPMs <- CombinedTPMs %>%
  mutate(fold_change = pmax(avg_tpm_HEK / avg_tpm_Hela, avg_tpm_Hela / avg_tpm_HEK),
         within_1.5x = fold_change <= 1.5)

CombinedTPMs0.1 <- CombinedTPMs %>%
  filter(avg_tpm_HEK >0.1, avg_tpm_Hela >0.1)
CombinedTPMs0.01 <- CombinedTPMs %>%
  filter(avg_tpm_HEK >0.01, avg_tpm_Hela >0.01)
CombinedTPMs10 <- CombinedTPMs %>%
  filter(avg_tpm_HEK >10, avg_tpm_Hela >10)


#get MANE select isoforms
MANEiso<-read_tsv("~/OneDrive - The Ohio State University/Singhlab/Bioinfo/PYMpaperFiguresUpdate/MANEexons.tsv")

CombinedTPMs0.1 <- CombinedTPMs0.1 %>% filter(target_id %in% MANEiso$ensembl_transcript_id)
CombinedTPMs0.01 <- CombinedTPMs0.01 %>% filter(target_id %in% MANEiso$ensembl_transcript_id)
CombinedTPMs10 <- CombinedTPMs10 %>% filter(target_id %in% MANEiso$ensembl_transcript_id)

library(ggpmisc)
library(ggpubr)
library(ggpointdensity)

CombinedTPMs0.1%>%ggplot(aes(avg_tpm_HEK,avg_tpm_Hela))+
  geom_pointdensity(shape=1)+
  scale_colour_gradient(low = "grey60", high = "black", na.value = NA)+
  scale_x_log10(labels = scales::comma)+scale_y_log10(labels = scales::comma)+theme_linedraw()+
  stat_cor(aes(label = after_stat(rr.label)),digits = 3,label.x = 0.5, label.y = 4.5)+
  geom_smooth(method = "lm", se = FALSE, color = "red", linetype = "dashed", linewidth = 0.5)+
  # geom_line(data = line_data, aes(x = avg_tpm_HEK, y = avg_tpm_Hela),
  #           inherit.aes = FALSE, color = "red", linewidth = 0.4, linetype = "dashed") +
  theme(aspect.ratio = 1,panel.grid.major = element_line(linetype = "blank"),
        panel.grid.minor = element_line(linetype = "blank"),legend.position = "none")+
        #,axis.text = element_text(size = 5),axis.title = element_text(size = 7))
  labs(title = "Average tpm of MANE select isoforms")


CombinedTPMs0.1 <- CombinedTPMs0.1 %>%
  mutate(within_2x = fold_change <= 2)
CombinedTPMs10 <- CombinedTPMs10 %>%
  mutate(within_2x = fold_change <= 2)

CombinedTPMs0.1%>%ggplot(aes(avg_tpm_HEK,avg_tpm_Hela))+
  geom_pointdensity(aes(color = within_2x),shape=1)+
  scale_color_manual(values = c("TRUE" = "green4", "FALSE" = "grey20"))+
  #scale_colour_gradient(low = "grey60", high = "black", na.value = NA)+
  scale_x_log10(labels = scales::comma)+scale_y_log10(labels = scales::comma)+theme_linedraw()+
  stat_cor(aes(label = after_stat(rr.label)),digits = 3,label.x = 0.5, label.y = 4.5)+
  geom_smooth(method = "lm", se = FALSE, color = "red", linetype = "dashed", linewidth = 0.5)+
  theme(aspect.ratio = 1,panel.grid.major = element_line(linetype = "blank"),
        panel.grid.minor = element_line(linetype = "blank"),legend.position = "none")+
  #,axis.text = element_text(size = 5),axis.title = element_text(size = 7))
  labs(title = "Average tpm of MANE select isoforms", caption = "n=7920")

##Filtering SE genes list
SEgenes<-MANEiso%>%filter(rank==1)
CombinedTPMs0.1%>%filter(target_id%in% SEgenes$ensembl_transcript_id)%>%dim()
CombinedTPMs0.1%>%filter(target_id%in% SEgenes$ensembl_transcript_id)%>%
  ggplot(aes(avg_tpm_HEK,avg_tpm_Hela))+
  geom_pointdensity(aes(color = within_2x),shape=1)+
  scale_color_manual(values = c("TRUE" = "green4", "FALSE" = "grey20"))+
  #scale_colour_gradient(low = "grey60", high = "black", na.value = NA)+
  scale_x_log10(labels = scales::comma)+scale_y_log10(labels = scales::comma)+theme_linedraw()+
  stat_cor(aes(label = after_stat(rr.label)),digits = 3,label.x = 0.5, label.y = 4.5)+
  geom_smooth(method = "lm", se = FALSE, color = "red", linetype = "dashed", linewidth = 0.5)+
  theme(aspect.ratio = 1,panel.grid.major = element_line(linetype = "blank"),
        panel.grid.minor = element_line(linetype = "blank"),legend.position = "none")+
  #,axis.text = element_text(size = 5),axis.title = element_text(size = 7))
  labs(title = "Average tpm of MANE select single exon isoforms", caption = "n=183")

names(CombinedTPMs0.1)[1]<-"ensembl_transcript_id"
CombinedTPMs0.1<-inner_join(CombinedTPMs0.1,MANEisoforms)

SEgeneswithin2X<-CombinedTPMs0.1%>%filter(within_2x == T &
                                            rank == 1)%>%dplyr::select(ensembl_gene_id)
write_tsv(SEgeneswithin2X,"~/OneDrive - The Ohio State University/Singhlab/Bioinfo/PYMpaperFiguresUpdate/Figure4/SEgeneswithin2X.tsv")

d<-SEgenesScatter%>%filter(Clip != 0)%>%filter(Geneid %in% SEgeneswithin2X$ensembl_gene_id)%>%
  ggplot(aes(Clip,RIPiT, label = Name))+
  geom_point(shape = 1, alpha =0.8, color="black",size=0.8)+geom_label_repel(size = 2)+
  geom_point(data = SEgenesScatter[SEgenesScatter$Name != "",], color = "red")+
  theme_linedraw()+
  labs(x= "CLIP-seq",y="RIPiT-seq")+
  theme(panel.grid.major = element_line(linetype = "blank"),
        panel.grid.minor = element_line(linetype = "blank"),
        aspect.ratio = 1)+
  stat_cor(aes(label = after_stat(rr.label)),label.x = 0, label.y = 2.7)+
  geom_smooth(method = lm, color = "grey", linewidth = 0.4, linetype="dashed", se = F)+
  scale_x_log10()+scale_y_log10()
d+ labs(title = "SE genes scatter: Clip vs RIPiT",
        subtitle = "Only genes within 2X")




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
b<-SEgenesScatter%>%
  ggplot(aes(RNA,RIPiT, label = Name))+
  geom_point(shape = 1, alpha =0.8, color="black", size=0.8)+geom_label_repel(size = 2)+
  geom_point(data = SEgenesScatter[SEgenesScatter$Name != "",], color = "red")+
  theme_linedraw()+
  labs(x= "RNA-seq",y="RIPiT-seq")+
  theme(panel.grid.major = element_line(linetype = "blank"),
        panel.grid.minor = element_line(linetype = "blank"),
        aspect.ratio = 1)+
  stat_cor(aes(label = after_stat(rr.label)),label.x = 0, label.y = 2.7)+
  geom_smooth(method = lm, color = "grey", linewidth = 0.4, linetype="dashed", se = F)+
  scale_x_log10()+scale_y_log10()
b+labs(title = "SE genes scatter: RNA vs RIPiT")
