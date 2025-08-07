##Reading Kovalak etal supplementary table to identify SE genes without introns
setwd("~/OneDrive - The Ohio State University/Singhlab/Bioinfo/PYMpaperFiguresUpdate/Figure4/Kovalak_etal/")
library(tidyverse)

library(biomaRt)
listMarts(host='https://apr2020.archive.ensembl.org')
ensembl100=useMart(host='https://apr2020.archive.ensembl.org', biomart='ENSEMBL_MART_ENSEMBL', dataset='hsapiens_gene_ensembl')
attributes <- listAttributes(ensembl100)
filters <- listFilters(ensembl100)

allgenes <- getBM(attributes = c("ensembl_gene_id","ensembl_transcript_id","chromosome_name",
                                 "external_gene_name","transcript_biotype"),
                  mart = ensembl100,
                  useCache = FALSE)

AnnInt<-read_tsv("13059_2021_2309_MOESM2_ESM.txt",skip = 5)
NovInt<-read_tsv("13059_2021_2309_MOESM3_ESM.txt",skip = 10)

filteredgenes<-allgenes%>%filter(!external_gene_name %in% AnnInt$Gene_Name)
filteredgenes<-filteredgenes%>%filter(!external_gene_name %in% NovInt$`Gene_Name(s)`)
filteredgenes<-filteredgenes%>%filter(transcript_biotype=="protein_coding")

table(filteredgenes$chromosome_name)
table(filteredgenes$transcript_biotype)

##Chacking if slcn counts have lncRNA
slcn<-read_csv("../../SLCN_counts/slcn-trimmed-last-allcans-final_NT_res.csv")
slcn$geneid<-str_split(slcn$...1, "\\-",simplify = T)[,1]
slcn$region<-str_split(slcn$...1, "\\-",simplify = T)[,2]

slcngenes<-unique(slcn$geneid)

allgenes <- getBM(attributes = c("ensembl_gene_id","gene_biotype"),
                  values = slcngenes,
                  filters = "ensembl_gene_id",
                  mart = ensembl100,
                  useCache = FALSE)
table(allgenes$gene_biotype)#######ONLY protein codiing genes

MANEselectgenes <- getBM(attributes = c("ensembl_gene_id","ensembl_transcript_id","transcript_mane_select"),
                  mart = ensembl100,
                  useCache = FALSE)
MANEselectgenes <- subset(MANEselectgenes, MANEselectgenes$transcript_mane_select != "")
Exons<-getBM(attributes = c("ensembl_gene_id","ensembl_transcript_id","ensembl_exon_id","rank","transcript_length"),
             filters = "ensembl_transcript_id",
             values = MANEselectgenes$ensembl_transcript_id,
             mart = ensembl100,
             useCache = FALSE)
Exons<-Exons%>%arrange(ensembl_transcript_id,desc(rank))%>%group_by(ensembl_transcript_id)%>%dplyr::slice(1)
SEgenes<-Exons%>%filter(rank == 1)
SEgenesNointrons<-read_tsv("single-exon_no_annotated_no_novel.bed", col_names = F)
SEgenesNointrons$ensembl_gene_id<-str_split(SEgenesNointrons$X4,"\\_",simplify = T)[,1]
intersect(SEgenes$ensembl_gene_id,SEgenesNointrons$ensembl_gene_id)

SEgenes_MANE_noInt<-inner_join(SEgenes,SEgenesNointrons[,"ensembl_gene_id"])
write_tsv(SEgenes_MANE_noInt,"SEgenes_MANE_noInt.tsv")

slcn<-slcn%>%filter(baseMean>10)
slcn<-slcn%>%filter(geneid%in%SEgenes_MANE_noInt$ensembl_gene_id)

NoInt3UTR<-read_tsv("Kovalak_etal/3UTR_no_annotatedIntrons_no_novel.bed", col_names = F)
NoInt3UTR$ensembl_gene_id<-str_split(NoInt3UTR$X4,"\\_",simplify = T)[,1]




