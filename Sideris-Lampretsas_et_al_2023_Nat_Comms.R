# Mouse Microglia - Neuroinflammation and Pain RNA-Seq_KCL_Analysis
# author: George Sideris-Lampretsas, Dr Gopuraja Dharmalingam, Dr Archana Bajpai and Dr David Collier
# date: 22/12/2021
  

# Packages

library(rtracklayer)
library(ggplot2)
library(Rsamtools)
library(tidyverse)
library(DESeq2)
library(AnnotationDbi)
library(org.Mm.eg.db)
library(DOSE)
library(pathview)
library(clusterProfiler)
library(ReactomePA)
library(RColorBrewer)
library(ggrepel)
library(pheatmap)
library(ggvenn)


# unique colors
colors25 = c(
  "dodgerblue2", "#E31A1C", # red
  "green4",
  "#6A3D9A", # purple
  "#FF7F00", # orange
  "black", "gold1",
  "skyblue2", "#FB9A99", # lt pink
  "palegreen2",
  "#CAB2D6", # lt purple
  "#FDBF6F", # lt orange
  "gray70", "khaki2",
  "maroon", "orchid1", "deeppink1", "blue1", "steelblue4",
  "darkturquoise", "green1", "yellow4", "yellow3",
  "darkorange4", "brown"
)

colors16 = c("dodgerblue2", "#E31A1C", "green4", "#6A3D9A", "#FF7F00", "black", "gold1", "skyblue2", "palegreen2", "#FDBF6F", "gray70", "maroon", "orchid1", "darkturquoise", "darkorange4", "brown")



featcountsummary <- read.table("C:~/Microglia_RNASeq.txt.summary", sep = "\t", header = T)

colnames(featcountsummary) = gsub("~.Microglia_RNA_seq.vendor_data.raw_data.R\\d+.|.bam", "", colnames(featcountsummary))
rownames(featcountsummary) <- featcountsummary[,1]
featcountsummary <- t(featcountsummary[,-1])

sampleinfo <- data.frame(
  Sample = paste0("R", c(1:28)),
  Condition = c(rep("TASTPM_KBxN", 4), rep("TASTPM_Control", 3), rep("WT_KBxN", 4), rep("WT_Control", 3),
                rep("TASTPM_KBxN", 4), rep("TASTPM_Control", 4), rep("WT_KBxN", 4), rep("WT_Control", 3)),
  Day = c(rep("Day5", 14), rep("Day30", 14))
)

sampleinfo$ConditionDay <- paste(sampleinfo$Condition, sampleinfo$Day, sep = "_")


featcountsummary <- merge(featcountsummary, sampleinfo[,-4], by.x = 0, by.y = 1)
featcountsummary <- featcountsummary[order(featcountsummary$Row.names),]
DT::datatable(featcountsummary[,c(1,16:17,2,3,9,10,15)])


# Annotation - Retrieve annotation for Ens v104 Mouse from GTF

gtf <- import.gff("~/Mus_musculus.GRCm39.105.gtf/Mus_musculus.GRCm39.105.gtf")
GeneSym <- data.frame(EnsID = gtf$gene_id,  Symbol = gtf$gene_name)
GeneSym <- GeneSym[!duplicated(GeneSym$EnsID),]


### Quality Assessment

genecounts = read.table("~/Microglia_RNASeq.txt", header = T, row.names=1)

colnames(genecounts) = gsub("~.Microglia_RNA_seq.vendor_data.raw_data.R\\d+.|.bam", "", colnames(genecounts))



# WT KBxN Day 5 vs WT Control Day 5

vsd1.sub_WK5 <- vsd2[ , vsd2$ConditionDay %in% c("WT_KBxN_Day5","WT_Control_Day5") ]
plotPCA(vsd1.sub_WK5, "ConditionDay")+ geom_text(aes(label = vsd1.sub_WK5$Sample))

WTKBxNDay5_vs_WTControlDay5 = results(dds1,contrast=c("ConditionDay", "WT_KBxN_Day5","WT_Control_Day5"))

##### alpha = 0.05
summary(WTKBxNDay5_vs_WTControlDay5, alpha=0.05)

DESeq2::plotMA(WTKBxNDay5_vs_WTControlDay5, main="WT KBxN Day 5 vs WT Control Day 5")


WTKBxNDay5_vs_WTControlDay5 = as.data.frame(WTKBxNDay5_vs_WTControlDay5)

WTKBxNDay5_vs_WTControlDay5 <- merge(WTKBxNDay5_vs_WTControlDay5, GeneSym, by.x = 0, by.y = 1, all.x = T)

WTKBxNDay5_vs_WTControlDay5 = WTKBxNDay5_vs_WTControlDay5[order(WTKBxNDay5_vs_WTControlDay5$padj,decreasing=F),]

WTKBxNDay5_vs_WTControlDay5$Diff_exprs <- "NO"

WTKBxNDay5_vs_WTControlDay5$Diff_exprs[WTKBxNDay5_vs_WTControlDay5$log2FoldChange > 1 & WTKBxNDay5_vs_WTControlDay5$padj < 0.05] <- "UP"

WTKBxNDay5_vs_WTControlDay5$Diff_exprs[WTKBxNDay5_vs_WTControlDay5$log2FoldChange < -1 & WTKBxNDay5_vs_WTControlDay5$padj < 0.05] <- "DOWN"

WTKBxNDay5_vs_WTControlDay5 <- WTKBxNDay5_vs_WTControlDay5 %>% 
  filter(!is.na(padj))

WTKBxNDay5_vs_WTControlDay5_sign <- WTKBxNDay5_vs_WTControlDay5[(WTKBxNDay5_vs_WTControlDay5$padj<0.05),]

DT::datatable(WTKBxNDay5_vs_WTControlDay5[,c(1,2,3,7,8,9)])


### Volcano plots

pwk5_final <- ggplot(WTKBxNDay5_vs_WTControlDay5) +
  geom_text_repel(aes(x = log2FoldChange, y = -log10(padj), label = ifelse(Symbol %in% c("Traf1",  "Il1r2", "Lgals1", "Lgals3", "Anxa1", "Mmp8", "Mmp9", "Bcl2", "Itgal", "Cd93", "Ccl7"), WTKBxNDay5_vs_WTControlDay5$Symbol, "" )),  size = 6.7, max.overlaps = Inf) +
  geom_point(aes(x = log2FoldChange, y = -log10(padj), colour = Diff_exprs)) +
  ggtitle("WT K/BxN vs WT Control Day 5") +
  xlab("log2 fold change") +
  ylab("-log10 adjusted p-value") +
  theme_light() + 
  geom_vline(xintercept=-1, size=0.5, linetype="dashed", colour="black") +
  geom_vline(xintercept= 1, size=0.5, linetype="dashed", colour="black") +
  geom_hline(yintercept= 1.301, size=0.5, linetype="dashed", colour="black") +
  theme(legend.position = "none",
        plot.title = element_text(size = rel(1.5), hjust = 0.5),
        axis.title = element_text(size = rel(1.25)))

print(pwk5_final + scale_colour_manual(values = c("darkviolet", "grey21","forestgreen")))



### Pathway Analysis

colnames(WTKBxNDay5_vs_WTControlDay5)[1] <- "ENSEMBL"

annotations_orgDb_WK5 <- AnnotationDbi::select(org.Mm.eg.db,
                                               keys = WTKBxNDay5_vs_WTControlDay5$Symbol,
                                               columns = c("ENSEMBL", "ENTREZID", "GENENAME"),
                                               keytype = "SYMBOL")

annotations_orgDb_WK5 <- annotations_orgDb_WK5 %>% filter(!is.na(annotations_orgDb_WK5$ENSEMBL))


### Gene Ontology

res_ids_WK5 <- inner_join(WTKBxNDay5_vs_WTControlDay5, annotations_orgDb_WK5, by=c("ENSEMBL"="ENSEMBL"))
res_ids_WK5 <- res_ids_WK5 %>% select(-"Symbol")

allOE_WK5 <- filter(res_ids_WK5, res_ids_WK5$log2FoldChange != 0)
allOE_WK5_genes <- as.character((allOE_WK5$ENSEMBL))

#####using p = 0.05
sigOE_WK5 <- filter(res_ids_WK5, res_ids_WK5$padj< 0.05)
sigOE_WK5_genes <- as.character(sigOE_WK5$ENSEMBL)

###RUN GO enrichment analysis
WK5_GO <- enrichGO(gene = sigOE_WK5_genes,
                   universe = allOE_WK5_genes,
                   keyType = "ENSEMBL",
                   OrgDb = org.Mm.eg.db,
                   ont = "BP",
                   pvalueCutoff = 1,
                   pAdjustMethod = "BH",
                   qvalueCutoff = 1,
                   readable = TRUE)

cluster_summary_WK5 <- data.frame(WK5_GO)
dotplot(WK5_GO, showCategory = 10)


### Reactome

#####using FDR = 0.05
res_entrez_WK5 <- dplyr::filter(res_ids_WK5, ENTREZID != "NA") %>% filter(res_ids_WK5$padj < 0.05)

res_entrez_WK5 <- res_entrez_WK5[which(duplicated(res_entrez_WK5$ENTREZID) == F),]

fold_changes_WK5 <- res_entrez_WK5$log2FoldChange
names(fold_changes_WK5) <- res_entrez_WK5$ENTREZID
fold_changes_WK5 <- sort(fold_changes_WK5, decreasing = TRUE)

gseaREACT_WK5 <- gsePathway(geneList = fold_changes_WK5, # ordered named vector of fold changes (Entrez IDs are the associated names)
                            organism = "mouse", # supported organisms listed below
                            minGSSize = 10, # minimum gene set size (# genes in set) - change to test more sets or recover sets with fewer # genes
                            pvalueCutoff = 1, # padj cutoff value
                            verbose = FALSE)

gseaREACT_results_WK5 <- gseaREACT_WK5@result
#GSEA
dotplot(gseaREACT_WK5, showCategory = 10)



# TASTPM KBxN Day 5 vs TASTPM Control Day 5

vsd1.sub_TK5 <- vsd2[ , vsd2$ConditionDay %in% c("TASTPM_KBxN_Day5","TASTPM_Control_Day5") ]

plotPCA(vsd1.sub_TK5, "ConditionDay")+ geom_text(aes(label = vsd1.sub_TK5$Sample))

TASTPMKBxNDay5_vs_TASTPMControlDay5 = results(dds1,contrast=c("ConditionDay", "TASTPM_KBxN_Day5","TASTPM_Control_Day5"))



##### alpha = 0.05
summary(TASTPMKBxNDay5_vs_TASTPMControlDay5, alpha=0.05)

DESeq2::plotMA(TASTPMKBxNDay5_vs_TASTPMControlDay5, main="TASTPM KBxN Day 5 vs TASTPM Control Day 5")

TASTPMKBxNDay5_vs_TASTPMControlDay5 = as.data.frame(TASTPMKBxNDay5_vs_TASTPMControlDay5)

TASTPMKBxNDay5_vs_TASTPMControlDay5 <- merge(TASTPMKBxNDay5_vs_TASTPMControlDay5, GeneSym, by.x = 0, by.y = 1, all.x = T)

TASTPMKBxNDay5_vs_TASTPMControlDay5 = TASTPMKBxNDay5_vs_TASTPMControlDay5[order(TASTPMKBxNDay5_vs_TASTPMControlDay5$padj,decreasing=F),]

TASTPMKBxNDay5_vs_TASTPMControlDay5$Diff_exprs <- "NO"
TASTPMKBxNDay5_vs_TASTPMControlDay5$Diff_exprs[TASTPMKBxNDay5_vs_TASTPMControlDay5$log2FoldChange > 1 & TASTPMKBxNDay5_vs_TASTPMControlDay5$padj < 0.05] <- "UP"
TASTPMKBxNDay5_vs_TASTPMControlDay5$Diff_exprs[TASTPMKBxNDay5_vs_TASTPMControlDay5$log2FoldChange < -1 & TASTPMKBxNDay5_vs_TASTPMControlDay5$padj < 0.05] <- "DOWN"


TASTPMKBxNDay5_vs_TASTPMControlDay5 <- TASTPMKBxNDay5_vs_TASTPMControlDay5 %>% 
  filter(!is.na(padj))


TASTPMKBxNDay5_vs_TASTPMControlDay5_sign <- TASTPMKBxNDay5_vs_TASTPMControlDay5[(TASTPMKBxNDay5_vs_TASTPMControlDay5$padj<0.05),]


DT::datatable(TASTPMKBxNDay5_vs_TASTPMControlDay5[,c(1,2,3,7,8,9)])


### Volcano plots

ptk5_final <- ggplot(TASTPMKBxNDay5_vs_TASTPMControlDay5)  +
  geom_text_repel(aes(x = log2FoldChange, y = -log10(padj), label = ifelse(Symbol %in% c("Ifitm2","Ifitm3", "Ifitm6", "Lgals3", "Pkm", "Pfkp", "Ccl7", "Gpi1", "Atp8b4", "Clec4n", "Socs3", "Stat4"), TASTPMKBxNDay5_vs_TASTPMControlDay5$Symbol, "" )),  size = 6.7, max.overlaps = Inf) +
  geom_point(aes(x = log2FoldChange, y = -log10(padj), colour = Diff_exprs)) +
  ggtitle("TASTPM K/BxN vs TASTPM Control Day 5") +
  xlab("log2 fold change") +
  ylab("-log10 adjusted p-value") +
  theme_light() + 
  geom_vline(xintercept=-1, size=0.5, linetype="dashed", colour="black") +
  geom_vline(xintercept= 1, size=0.5, linetype="dashed", colour="black") +
  geom_hline(yintercept= 1.301, size=0.5, linetype="dashed", colour="black") +
  theme(legend.position = "none",
        plot.title = element_text(size = rel(1.5), hjust = 0.5),
        axis.title = element_text(size = rel(1.25)))

print(ptk5_final + scale_colour_manual(values = c("dodgerblue2", "grey21","firebrick3")))




### Pathway Analysis
colnames(TASTPMKBxNDay5_vs_TASTPMControlDay5)[1] <- "ENSEMBL"

annotations_orgDb_TK5 <- AnnotationDbi::select(org.Mm.eg.db,
                                               keys = TASTPMKBxNDay5_vs_TASTPMControlDay5$Symbol,
                                               columns = c("ENSEMBL", "ENTREZID", "GENENAME"),
                                               keytype = "SYMBOL")

annotations_orgDb_TK5 <- annotations_orgDb_TK5 %>% filter(!is.na(annotations_orgDb_TK5$ENSEMBL))


### Gene Ontology
res_ids_TK5 <- inner_join(TASTPMKBxNDay5_vs_TASTPMControlDay5, annotations_orgDb_TK5, by=c("ENSEMBL"="ENSEMBL"))
res_ids_TK5 <- res_ids_TK5 %>% select(-"Symbol")

allOE_TK5 <- filter(res_ids_TK5, res_ids_TK5$log2FoldChange != 0)
allOE_TK5_genes <- as.character((allOE_TK5$ENSEMBL))

#####using p = 0.05
sigOE_TK5 <- filter(res_ids_TK5, res_ids_TK5$padj< 0.05)
sigOE_TK5_genes <- as.character(sigOE_TK5$ENSEMBL)

###RUN GO enrichment analysis
TK5_GO <- enrichGO(gene = sigOE_TK5_genes,
                   universe = allOE_TK5_genes,
                   keyType = "ENSEMBL",
                   OrgDb = org.Mm.eg.db,
                   ont = "BP",
                   pvalueCutoff = 1,
                   pAdjustMethod = "BH",
                   qvalueCutoff = 1,
                   readable = TRUE)

cluster_summary_TK5 <- data.frame(TK5_GO)
dotplot(TK5_GO, showCategory = 10)


### Reactome

#####using FDR = 0.05
res_entrez_TK5 <- dplyr::filter(res_ids_TK5, ENTREZID != "NA") %>% filter(res_ids_TK5$padj < 0.05)
res_entrez_TK5 <- res_entrez_TK5[which(duplicated(res_entrez_TK5$ENTREZID) == F),]
fold_changes_TK5 <- res_entrez_TK5$log2FoldChange
names(fold_changes_TK5) <- res_entrez_TK5$ENTREZID
fold_changes_TK5 <- sort(fold_changes_TK5, decreasing = TRUE)

gseaREACT_TK5 <- gsePathway(geneList = fold_changes_TK5, # ordered named vector of fold changes (Entrez IDs are the associated names)
                            organism = "mouse", # supported organisms listed below
                            minGSSize = 10, # minimum gene set size (# genes in set) - change to test more sets or recover sets with fewer # genes
                            pvalueCutoff = 1, # padj cutoff value
                            verbose = FALSE)

gseaREACT_results_TK5 <- gseaREACT_TK5@result
#GSEA
dotplot(gseaREACT_TK5, showCategory = 10)


## KBxN only genes
WK5_significant <- WTKBxNDay5_vs_WTControlDay5 %>% filter(WTKBxNDay5_vs_WTControlDay5$padj<0.05 & abs(WTKBxNDay5_vs_WTControlDay5$log2FoldChange) > 1) 

TK5_significant <- TASTPMKBxNDay5_vs_TASTPMControlDay5 %>% filter(TASTPMKBxNDay5_vs_TASTPMControlDay5$padj<0.05 &  abs(TASTPMKBxNDay5_vs_TASTPMControlDay5$log2FoldChange) > 1)

KBxN_D5_only_genes <- intersect(WK5_significant$Symbol, TK5_significant$Symbol)


WT_KBxN_D5_only_genes <- setdiff(WK5_significant$Symbol, TK5_significant$Symbol)

WT_KBxN_D5_only <- WK5_significant %>% filter(WK5_significant$Symbol %in% WT_KBxN_D5_only_genes)


TASTPM_KBxN_D5_only_genes <- setdiff(TK5_significant$Symbol, WK5_significant$Symbol)

TASTPM_KBxN_D5_only <- TK5_significant %>% filter(TK5_significant$Symbol %in% TASTPM_KBxN_D5_only_genes)


### Comparison between logFC WT vs TASTPM K/BxN Day 5

WTKBxNDay5_vs_WTControlDay5$Genotype <- "NA"
WTKBxNDay5_vs_WTControlDay5$Genotype[WTKBxNDay5_vs_WTControlDay5$Symbol %in% KBxN_D5_only_genes] <- "COMMON"
WTKBxNDay5_vs_WTControlDay5$Genotype[WTKBxNDay5_vs_WTControlDay5$Symbol %in% WT_KBxN_D5_only_genes] <- "WT"
WTKBxNDay5_vs_WTControlDay5$Genotype[WTKBxNDay5_vs_WTControlDay5$Symbol %in% TASTPM_KBxN_D5_only_genes] <- "TASTPM"

TASTPMKBxNDay5_vs_TASTPMControlDay5$Genotype <- "NA"
TASTPMKBxNDay5_vs_TASTPMControlDay5$Genotype[TASTPMKBxNDay5_vs_TASTPMControlDay5$Symbol %in% KBxN_D5_only_genes] <- "COMMON"
TASTPMKBxNDay5_vs_TASTPMControlDay5$Genotype[TASTPMKBxNDay5_vs_TASTPMControlDay5$Symbol %in% WT_KBxN_D5_only_genes] <- "WT"
TASTPMKBxNDay5_vs_TASTPMControlDay5$Genotype[TASTPMKBxNDay5_vs_TASTPMControlDay5$Symbol %in% TASTPM_KBxN_D5_only_genes] <- "TASTPM"


significant_FC_Day5 <- inner_join(WTKBxNDay5_vs_WTControlDay5, TASTPMKBxNDay5_vs_TASTPMControlDay5, by=c("ENSEMBL"="ENSEMBL")) %>% rename(Genotype = Genotype.x)

significant_FC_Day5_WTK5 <- significant_FC_Day5[((abs(significant_FC_Day5$log2FoldChange.x) > 1 & significant_FC_Day5$padj.x < 0.05 ) | (abs(significant_FC_Day5$log2FoldChange.y) > 1 & significant_FC_Day5$padj.y < 0.05)),]

significant_FC_Day5_WTK5 <- filter(significant_FC_Day5_WTK5, significant_FC_Day5_WTK5$Genotype != "NA")


p3 <- ggplot(significant_FC_Day5_WTK5, aes(x = significant_FC_Day5_WTK5$log2FoldChange.x, y =  significant_FC_Day5_WTK5$log2FoldChange.y, colour = Genotype)) +
  geom_point() +
  geom_hline(yintercept=1, size=0.5,linetype="dashed", colour="black") +
  geom_hline(yintercept=-1, size=0.5,linetype="dashed", colour="black") +
  geom_vline(xintercept=1, size=0.5,linetype="dashed", colour="black") +
  geom_vline(xintercept=-1, size=0.5,linetype="dashed", colour="black") +
  geom_abline(slope = 1, intercept=0, size=0.5, linetype="dashed", colour="grey") +
  geom_label_repel(aes(label = ifelse(significant_FC_Day5_WTK5$Symbol.x %in% c("Lgals3"), significant_FC_Day5_WTK5$Symbol.x, "" )), show_guide=F, max.overlaps = Inf) + #show_guide eliminates a's in the legend
  xlab("logFC WT K/BxN Day 5") +
  ylab("logFC TASTPM K/BxN Day 5")+
  theme_classic()

print(p3 + scale_colour_manual(values = c("royalblue","firebrick3",  "forestgreen")))




