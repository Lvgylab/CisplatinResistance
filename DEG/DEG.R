#CDDP25ug/ml vs con 20h
rm(list = ls())
exprSet <- data.table::fread("data/exprSet_counts_cddp25vscon_20H.csv",data.table = F)
exprSet<-exprSet[1:7]
colnames(exprSet)[2:7] <- c("CON_1","CON_2","CON_3","CDDP_1","CDDP_2","CDDP_3")
metadata <- data.table::fread("data/metadata_cddp25vscon.txt",data.table = F,header = F)
rownames(metadata) <- metadata[,1]
metadata <- metadata[colnames(exprSet)[-1],]
metadata$group <- ifelse(grepl("CDDP",metadata$V2),"treat","control")
metadata = metadata[,c(1,3)]
colnames(metadata) <- c("sample","group")
save(metadata,file = "output/metadata_cddp25vscon_20H.Rdata")
library(DESeq2)
dds <-DESeqDataSetFromMatrix(countData=exprSet, 
                             colData=metadata, 
                             design=~group,
                             tidy=TRUE) 
dds <- dds[rowSums(counts(dds))>1,]
vsd <- vst(dds, blind = FALSE)
exprSet_vst <- as.data.frame(assay(vsd))
test <- exprSet_vst[1:10,1:6]
save(exprSet_vst,file = "output/exprSet_vst_cddp25vscon_20h.Rdata")
dds <- DESeq(dds)
contrast=c("group", "treat", "control")
dd1 <- results(dds, contrast=contrast, alpha = 0.05)
plotMA(dd1, ylim=c(-5,5))
dd2 <- lfcShrink(dds,contrast=contrast, res=dd1,type="ashr")
plotMA(dd2, ylim=c(-5,5))
library(dplyr)
library(tibble)
library(tidyr)
res <- dd2 %>% 
  as.data.frame() %>% 
  rownames_to_column("gene_id") 
library(AnnotationDbi)
library(org.Hs.eg.db)
res$symbol <- mapIds(org.Hs.eg.db,
                     keys=res$gene_id,
                     column="SYMBOL",
                     keytype="ENSEMBL",
                     multiVals="first")
res$entrez <- mapIds(org.Hs.eg.db,
                     keys=res$gene_id,
                     column="ENTREZID",
                     keytype="ENSEMBL",
                     multiVals="first")
colnames(res) <- c("gene_id","baseMean","logFC","lfcSE","P.Value","adj.P.Val","gene","entrez")
save(res,file = "output/DEseq2_CDDP25vscon_Diff_20h.Rdata")
write.csv(res,file = "output/DEseq2_CDDP25vscon_20H_Diff.csv")

colnames(res)[7]="symbol"
library(dplyr)
ensemble_symbol = res %>% 
  dplyr::select(gene_id,symbol) %>% 
  filter(symbol !="") 
exprSet <- cbind(gene_id=rownames(exprSet_vst),exprSet_vst)
exprSet <- merge(ensemble_symbol,exprSet,by="gene_id")
exprSet <- exprSet %>% 
  dplyr::select(-gene_id) %>% 
  mutate(newcolumn = rowMeans(.[,-1])) %>% 
  arrange(desc(newcolumn)) %>% 
  distinct(symbol,.keep_all = T) %>% 
  dplyr::select(-newcolumn)
rownames(exprSet) <- exprSet[,1]
exprSet <- exprSet[,-1]
exprSet <- t(exprSet)
exprSet <- as.data.frame(exprSet)
test <- exprSet[,1:10]
exprSet <- cbind(group= metadata$group,exprSet)
test <- exprSet[,1:10]
save(exprSet,file = "output/exprSet_tidy_CDDP25vscon_20h.Rdata")

my_comparisons <- list(
  c("treat", "control")
)
genelist=c("MFN1","MFN2","OPA1", "DNM1L","MFF","FIS1")#mitochondrial dynamics related genes
genelist=c("ROCK1","RHOA","ACTB","PFN1","LIMK1","CFL1","INF2","MYH9","MYH10","MYLK")#cytoskeleton associated genes
genelist=c("MTOR","PRKAA1","ULK1","BCL2L1","BECN1","ATG5","ATG12","SQSTM1", "MAP1LC3B","LAMP2")#autophagy associated genes
genelist=c("SEC62","RETREG1","RTN3", "TEX264","CCPG1","ATL3","CANX","SQSTM1", "MAP1LC3B")#ER-phagy associated genes
genelist=c("XBP1","ATF4","DDIT3", "ATF6","HSPA5","ERN1","EIF2AK3")#ER stress associated genes
data <- exprSet[,c("group",genelist)]
library(tidyr)
library(ggplot2)
library(ggpubr)
data <- data %>% 
  pivot_longer(cols=-1,
               names_to= "gene",
               values_to = "expression")
data$gene=factor(data$gene,levels=genelist)
ggplot(data = data,aes(x=group,y=expression,fill=group))+
  geom_boxplot()+
  geom_jitter()+
  theme_bw()+
  facet_grid(.~gene)+
  scale_fill_manual(values=c("#FFD034", "#0057e7"))+
  stat_compare_means(comparisons = my_comparisons, label = "p.signif", method = "t.test")
graph2tif(file="output/CDDP25/boxplot_mitochondrial_dynamics.tiff")
graph2tif(file="output/CDDP25/boxplot_regulation_cell_skeleton.tiff")
graph2tif(file="output/CDDP25/boxplot_autophagy.tiff")
graph2tif(file="output/CDDP25/boxplot_ERphagy.tiff")
graph2tif(file="output/CDDP25/boxplot_ERstress.tiff")

####heatmap
TPM = read.csv(file = "data/TPM_ALL_20h.csv",check.names = F)
TPM = filter(TPM,symbol!="")
heatgene_up <- res %>% 
  filter(adj.P.Val< 0.05) %>% 
  filter(logFC >2) %>% 
  arrange(desc(logFC)) %>% 
  arrange(adj.P.Val)
heatgene_up_gene = heatgene_up[1:250,7]
heatgene_down <- res %>% 
  filter(adj.P.Val< 0.05) %>% 
  filter(logFC <(-2))%>% 
  arrange(desc(abs(logFC))) %>% 
  arrange(adj.P.Val)
heatgene_down_gene = heatgene_down[1:250,7]
heatgene=c(heatgene_up_gene,heatgene_down_gene)
TPM=distinct(TPM,symbol,.keep_all=TRUE)
rownames(TPM)=TPM$symbol
heatTPM=TPM[heatgene,]
heatTPM=filter(heatTPM,symbol!="")
heatTPM=heatTPM[,-c(1,11)]
library(pheatmap)
metadata <- data.table::fread("data/metadata_TPM_all_20h.txt",data.table = F,header = F)
annotation_col <- data.frame(group=metadata$V2)
rownames(annotation_col) <- metadata$V1
ann_colors = list(group = c(control = "#72C4AF", CDDP25 = "#51B7DC",CDDP50 = "#6655A2") )
pheatmap(heatTPM, 
         cluster_rows = TRUE,
         cluster_cols = F,
         annotation_col =annotation_col, 
         annotation_legend=TRUE, 
         show_rownames = F,
         show_colnames = F,
         scale = "row", 
         color =colorRampPalette(c("#0c2c84","#225ea8","#1d91c0","#41b6c4","#7fcdbb","#c7e9b4","#ffffcc"))(100),
         annotation_colors = ann_colors,
         cellwidth = 25, cellheight = 0.4,
         fontsize = 10)
library(export)
graph2tif(file="output/heatmap_logfc2.tiff")

####KEGG barplot
library(dplyr)
diffgene <- res %>% 
  filter(symbol!="") %>% 
  filter(adj.P.Val< 0.05) %>% 
  filter(abs(logFC) >1) %>% 
  arrange(desc(abs(logFC))) %>% 
  arrange(adj.P.Val) 
write.csv(diffgene,file="output/CDDP25_20H_diffgene_ALL_logfc1.csv")

library(clusterProfiler)
EGG <- enrichKEGG(gene = diffgene$entrez,
                  organism = 'hsa',
                  pvalueCutoff = 0.05,
                  use_internal_data =T)
KEGG_read <-  setReadable(EGG, 'org.Hs.eg.db', 'ENTREZID')
KEGG_df <- as.data.frame(KEGG_read)

library(stringr)
KEGG<- KEGG_df[1:10,]
ggplot(KEGG, aes(x=Description,y=Count,fill=p.adjust))+
  geom_bar(stat = "identity")+
  scale_x_discrete(limits=rev(KEGG[,2]),labels=function(x) str_wrap(x,width = 30))+
  coord_flip()+
  labs(x="",y="",title = "KEGG")+
  scale_fill_gradient(low="#0057e7",high ="#FFD500")+
  theme(panel.background = element_rect(fill = "transparent",color = "gray"),
        axis.text.y = element_text(color = "black",size=12),
        axis.text.x = element_text(color = "black",size=12))
library(export)
graph2tif(file="output/CDDP25/KEGG_barplot.tiff")

####network of various pathway and process
library(dplyr)
diffgene <- res %>% 
  filter(symbol!="") %>% 
  filter(adj.P.Val< 0.05) %>% 
  filter(abs(logFC) >2) %>% 
  arrange(desc(abs(logFC))) %>% 
  arrange(adj.P.Val) 
save(diffgene,file = "output/diffgene_CDDP25vscon_LOGFC2.Rdata")
write.csv(diffgene,file = "output/diffgene_CDDP25vscon_LOGFC2.csv")
#Select differently expressed genes that (logFC) >2 & adj.P.Val<0.001, following by enrichment analysis using Metascape
#(https://metascape.org/gp/index.html#/main/step1), the result is visualized by cytoscape.

####GSEA
allDiff=res
nrow(allDiff)
allDiff=na.omit(allDiff)

library(dplyr)
allDiff=as.data.frame(allDiff)
allDiff=allDiff %>% distinct(symbol, .keep_all = TRUE)
rownames(allDiff)=allDiff$symbol
head(allDiff)
allDiff=allDiff[,2:6]

library(dplyr)
library(tibble)
data <- allDiff %>%
  rownames_to_column("symbol") %>% 
  dplyr::select(symbol,logFC,adj.P.Val) %>%
  dplyr::mutate(
    offset = logFC,
    offset_abs = abs(logFC),
  ) %>% 
  dplyr::mutate(rank = rank(-offset), priority = offset) %>%
  dplyr::select(priority, rank, symbol) %>% 
  column_to_rownames("symbol")

#Rdata
library(Pi)
index <- c("GWAS2EF", "GWAS_LD", "IlluminaHumanHT",
           "IlluminaOmniExpress", "ig.DO",
           "ig.EF", "ig.GOBP", "ig.GOCC", "ig.GOMF", "ig.HPCM", "ig.HPMA",
           "ig.HPMI", "ig.HPPA",
           "ig.MP", "org.Hs.eg", "org.Hs.egDGIdb", "org.Hs.egDO", "org.Hs.egGOBP",
           "org.Hs.egGOCC", "org.Hs.egGOMF", "org.Hs.egHPCM", "org.Hs.egHPMA",
           "org.Hs.egHPMI",
           "org.Hs.egHPPA", "org.Hs.egMP", "org.Hs.egMsigdbC1",
           "org.Hs.egMsigdbC2BIOCARTA",
           "org.Hs.egMsigdbC2CGP", "org.Hs.egMsigdbC2CPall",
           "org.Hs.egMsigdbC2CP",
           "org.Hs.egMsigdbC2KEGG", "org.Hs.egMsigdbC2REACTOME",
           "org.Hs.egMsigdbC3MIR",
           "org.Hs.egMsigdbC3TFT", "org.Hs.egMsigdbC4CGN", "org.Hs.egMsigdbC4CM",
           "org.Hs.egMsigdbC5BP", "org.Hs.egMsigdbC5CC", "org.Hs.egMsigdbC5MF",
           "org.Hs.egMsigdbC6", "org.Hs.egMsigdbC7", "org.Hs.egMsigdbH",
           "org.Hs.egPS",
           "org.Hs.egSF", "org.Hs.egPfam", "org.Hs.string", "org.Hs.PCommons_DN",
           "org.Hs.PCommons_UN")
dir.create("ontology_Rdata/")
for (i in index) {
  print(i)
  geneset = xRDataLoader(RData=i)
  assign(i,geneset)
  save(list = i,file = paste0("ontology_Rdata/",i,".Rdata"))
}


# GSEA using MsigdbH
####################
### GSEA 分析
library(Pi)
eGSEA <-
  Pi::xPierGSEA(
    data,
    ontology = "MsigdbC2KEGG",
    size.range = c(10, 500),
    nperm = 20000,
    fast = F,
    RData.location = paste0(getwd(),"/ontology_Rdata"))
eGSEA$df_summary$nes <- -(eGSEA$df_summary$nes)
test=eGSEA$df_summary
ls_gp <- xGSEAdotplot(eGSEA, top=10, signature=T,colormap = "Navy-MediumBlue-HotPink-DeepPink-MediumVioletRed",
                      subtitle ="both",leading = T,leading.size = 2)
ls_gp
library(export)
graph2tif(file="output/CDDP25/GSEA_Pi_KEGG_REG_CYTOSKELETON.tif")
graph2ppt(file="output/CDDP25/GSEA_Pi_KEGG_REG_CYTOSKELETON.pptx")

#CDDP50ug/ml vs con 20h
rm(list = ls())
exprSet <- data.table::fread("data/exprSet_counts_cddp50vscon_20H.csv",data.table = F)
exprSet<-exprSet[1:7]
colnames(exprSet)[2:7] <- c("CON_1","CON_2","CON_3","CDDP_1","CDDP_2","CDDP_3")
metadata <- data.table::fread("data/metadata_cddp50vscon.txt",data.table = F,header = F)
rownames(metadata) <- metadata[,1]
metadata <- metadata[colnames(exprSet)[-1],]
metadata$group <- ifelse(grepl("CDDP",metadata$V2),"treat","control")
metadata = metadata[,c(1,3)]
colnames(metadata) <- c("sample","group")
save(metadata,file = "output/metadata_cddp50vscon_20H.Rdata")
library(DESeq2)
dds <-DESeqDataSetFromMatrix(countData=exprSet, 
                             colData=metadata, 
                             design=~group,
                             tidy=TRUE) 
dds <- dds[rowSums(counts(dds))>1,]
vsd <- vst(dds, blind = FALSE)
exprSet_vst <- as.data.frame(assay(vsd))
test <- exprSet_vst[1:10,1:6]
save(exprSet_vst,file = "output/exprSet_vst_cddp50vscon.Rdata")
dds <- DESeq(dds)
contrast=c("group", "treat", "control")
dd1 <- results(dds, contrast=contrast, alpha = 0.05)
plotMA(dd1, ylim=c(-5,5))
dd2 <- lfcShrink(dds,contrast=contrast, res=dd1,type="ashr")
plotMA(dd2, ylim=c(-5,5))
library(dplyr)
library(tibble)
library(tidyr)
res <- dd2 %>% 
  as.data.frame() %>% 
  rownames_to_column("gene_id") 
library(AnnotationDbi)
library(org.Hs.eg.db)
res$symbol <- mapIds(org.Hs.eg.db,
                     keys=res$gene_id,
                     column="SYMBOL",
                     keytype="ENSEMBL",
                     multiVals="first")
res$entrez <- mapIds(org.Hs.eg.db,
                     keys=res$gene_id,
                     column="ENTREZID",
                     keytype="ENSEMBL",
                     multiVals="first")
colnames(res) <- c("gene_id","baseMean","logFC","lfcSE","P.Value","adj.P.Val","gene","entrez")
save(res,file = "output/DEseq2_CDDP50vscon_20H_Diff.Rdata")
write.csv(res,file = "output/DEseq2_CDDP50vscon_20H_Diff.csv")

colnames(res)[7]="symbol"
library(dplyr)
ensemble_symbol = res %>% 
  dplyr::select(gene_id,symbol) %>% 
  filter(symbol !="") 
exprSet <- cbind(gene_id=rownames(exprSet_vst),exprSet_vst)
exprSet <- merge(ensemble_symbol,exprSet,by="gene_id")
exprSet <- exprSet %>% 
  dplyr::select(-gene_id) %>% 
  mutate(newcolumn = rowMeans(.[,-1])) %>% 
  arrange(desc(newcolumn)) %>% 
  distinct(symbol,.keep_all = T) %>% 
  dplyr::select(-newcolumn)
rownames(exprSet) <- exprSet[,1]
exprSet <- exprSet[,-1]
exprSet <- t(exprSet)
exprSet <- as.data.frame(exprSet)
test <- exprSet[,1:10]
exprSet <- cbind(group= metadata$group,exprSet)
test <- exprSet[,1:10]
save(exprSet,file = "output/exprSet_tidy_CDDP50vscon_20h.Rdata")

my_comparisons <- list(
  c("treat", "control")
)
genelist=c("MFN1","MFN2","OPA1", "DNM1L","MFF","FIS1")#mitochondrial dynamics related genes
genelist=c("ROCK1","RHOA","ACTB","PFN1","LIMK1","CFL1","INF2","MYH9","MYH10","MYLK")#cytoskeleton associated genes
genelist=c("MTOR","PRKAA1","ULK1","BCL2L1","BECN1","ATG5","ATG12","SQSTM1", "MAP1LC3B","LAMP2")#autophagy associated genes
genelist=c("SEC62","RETREG1","RTN3", "TEX264","CCPG1","ATL3","CANX","SQSTM1", "MAP1LC3B")#ER-phagy associated genes
genelist=c("XBP1","ATF4","DDIT3", "ATF6","HSPA5","ERN1","EIF2AK3")#ER stress associated genes
data <- exprSet[,c("group",genelist)]
library(tidyr)
library(ggplot2)
library(ggpubr)
data <- data %>% 
  pivot_longer(cols=-1,
               names_to= "gene",
               values_to = "expression")
data$gene=factor(data$gene,levels=genelist)
ggplot(data = data,aes(x=group,y=expression,fill=group))+
  geom_boxplot()+
  geom_jitter()+
  theme_bw()+
  facet_grid(.~gene)+
  scale_fill_manual(values=c("#FFD034", "#0057e7"))+
  stat_compare_means(comparisons = my_comparisons, label = "p.signif", method = "t.test")
graph2tif(file="output/CDDP25/boxplot_mitochondrial_dynamics.tiff")
graph2tif(file="output/CDDP25/boxplot_regulation_cell_skeleton.tiff")
graph2tif(file="output/CDDP25/boxplot_autophagy.tiff")
graph2tif(file="output/CDDP25/boxplot_ERphagy.tiff")
graph2tif(file="output/CDDP25/boxplot_ERstress.tiff")

library(dplyr)
diffgene <- res %>% 
  filter(symbol!="") %>% 
  filter(adj.P.Val< 0.05) %>% 
  filter(abs(logFC) >1) %>% 
  arrange(desc(abs(logFC))) %>% 
  arrange(adj.P.Val) 
write.csv(diffgene,file = "output/diffgene_CDDP50vscon_20h_LOGFC1.csv")

#CDDP25ug/ml vs con 6h
rm(list = ls())
exprSet <- data.table::fread("data/exprSet_counts_cddp25vscon_6h.csv",data.table = F)
exprSet<-exprSet[1:7]
colnames(exprSet)[2:7] <- c("CON_1","CON_2","CON_3","CDDP_1","CDDP_2","CDDP_3")
metadata <- data.table::fread("data/metadata_cddp25vscon _6h.txt",data.table = F,header = F)
rownames(metadata) <- metadata[,1]
metadata <- metadata[colnames(exprSet)[-1],]
metadata$group <- ifelse(grepl("CDDP",metadata$V2),"treat","control")
metadata = metadata[,c(1,3)]
colnames(metadata) <- c("sample","group")
save(metadata,file = "output/metadata_cddp25vscon_6H.Rdata")
library(DESeq2)
dds <-DESeqDataSetFromMatrix(countData=exprSet, 
                             colData=metadata, 
                             design=~group,
                             tidy=TRUE) 
dds <- dds[rowSums(counts(dds))>1,]
vsd <- vst(dds, blind = FALSE)
exprSet_vst <- as.data.frame(assay(vsd))
test <- exprSet_vst[1:10,1:6]
save(exprSet_vst,file = "output/exprSet_vst_cddp25vscon_6h.Rdata")
dds <- DESeq(dds)
contrast=c("group", "treat", "control")
dd1 <- results(dds, contrast=contrast, alpha = 0.05)
plotMA(dd1, ylim=c(-5,5))
dd2 <- lfcShrink(dds,contrast=contrast, res=dd1,type="ashr")
plotMA(dd2, ylim=c(-5,5))
library(dplyr)
library(tibble)
library(tidyr)
res <- dd2 %>% 
  as.data.frame() %>% 
  rownames_to_column("gene_id") 
library(AnnotationDbi)
library(org.Hs.eg.db)
res$symbol <- mapIds(org.Hs.eg.db,
                     keys=res$gene_id,
                     column="SYMBOL",
                     keytype="ENSEMBL",
                     multiVals="first")
res$entrez <- mapIds(org.Hs.eg.db,
                     keys=res$gene_id,
                     column="ENTREZID",
                     keytype="ENSEMBL",
                     multiVals="first")
colnames(res) <- c("gene_id","baseMean","logFC","lfcSE","P.Value","adj.P.Val","gene","entrez")
save(res,file = "output/DEseq2_CDDP25vscon_6H_Diff.Rdata")

colnames(res)[7]="symbol"
library(dplyr)
ensemble_symbol = res %>% 
  dplyr::select(gene_id,symbol) %>% 
  filter(symbol !="") 
exprSet <- cbind(gene_id=rownames(exprSet_vst),exprSet_vst)
exprSet <- merge(ensemble_symbol,exprSet,by="gene_id")
exprSet <- exprSet %>% 
  dplyr::select(-gene_id) %>% 
  mutate(newcolumn = rowMeans(.[,-1])) %>% 
  arrange(desc(newcolumn)) %>% 
  distinct(symbol,.keep_all = T) %>% 
  dplyr::select(-newcolumn)
rownames(exprSet) <- exprSet[,1]
exprSet <- exprSet[,-1]
exprSet <- t(exprSet)
exprSet <- as.data.frame(exprSet)
test <- exprSet[,1:10]
exprSet <- cbind(group= metadata$group,exprSet)
test <- exprSet[,1:10]
save(exprSet,file = "output/exprSet_tidy_CDDP25vscon_6h.Rdata")

my_comparisons <- list(
  c("treat", "control")
)
genelist=c("MFN1","MFN2","OPA1", "DNM1L","MFF","FIS1")#mitochondrial dynamics related genes
genelist=c("ROCK1","RHOA","ACTB","PFN1","LIMK1","CFL1","INF2","MYH9","MYH10","MYLK")#cytoskeleton associated genes
genelist=c("MTOR","PRKAA1","ULK1","BCL2L1","BECN1","ATG5","ATG12","SQSTM1", "MAP1LC3B","LAMP2")#autophagy associated genes
genelist=c("SEC62","RETREG1","RTN3", "TEX264","CCPG1","ATL3","CANX","SQSTM1", "MAP1LC3B")#ER-phagy associated genes
genelist=c("XBP1","ATF4","DDIT3", "ATF6","HSPA5","ERN1","EIF2AK3")#ER stress associated genes
data <- exprSet[,c("group",genelist)]
library(tidyr)
library(ggplot2)
library(ggpubr)
data <- data %>% 
  pivot_longer(cols=-1,
               names_to= "gene",
               values_to = "expression")
data$gene=factor(data$gene,levels=genelist)
ggplot(data = data,aes(x=group,y=expression,fill=group))+
  geom_boxplot()+
  geom_jitter()+
  theme_bw()+
  facet_grid(.~gene)+
  scale_fill_manual(values=c("#FFD034", "#0057e7"))+
  stat_compare_means(comparisons = my_comparisons, label = "p.signif", method = "t.test")
graph2tif(file="output/CDDP25/boxplot_mitochondrial_dynamics.tiff")
graph2tif(file="output/CDDP25/boxplot_regulation_cell_skeleton.tiff")
graph2tif(file="output/CDDP25/boxplot_autophagy.tiff")
graph2tif(file="output/CDDP25/boxplot_ERphagy.tiff")
graph2tif(file="output/CDDP25/boxplot_ERstress.tiff")

####normalize and adjust batch effect
library("sva")
exprSet <- read.table("data/rawcount1.txt",comment.char="!",stringsAsFactors=F,header=T)
rownames(exprSet) <- exprSet[,1]
exprSet <- exprSet[,-1]
pheno<- read.table("data/batch.txt", header=T, sep="\t")
lograwcount1<- log2(exprSet+1)
batch <-  pheno$batch
adjustrawcount1<- ComBat(lograwcount1,batch = batch)
write.csv(adjustrawcount1,file="data/CDDP 20H 6H counts adjust2.csv",quote=F, sep="\t", row.names=T, col.names=T)

library("sva")
exprSet <- read.table("data/rawcount2.txt",comment.char="!",stringsAsFactors=F,header=T)
rownames(exprSet) <- exprSet[,1]
exprSet <- exprSet[,-1]
adjustrawcount2<- log2(exprSet+1)
write.csv(adjustrawcount2,file="data/otherCells_CDDP counts adjust.csv",quote=F, sep="\t", row.names=T, col.names=T)

###GSVA for KEGG
rm(list = ls())
exprSet=data.table::fread("data/CDDP 20H 6H counts adjust.csv",data.table = F)
res <- exprSet
library(dplyr)
library(tibble)
library(tidyr)
library(AnnotationDbi)
library(org.Hs.eg.db)
res$symbol <- mapIds(org.Hs.eg.db,
                     keys=res$V1,
                     column="SYMBOL",
                     keytype="ENSEMBL",
                     multiVals="first")
res$entrez <- mapIds(org.Hs.eg.db,
                     keys=res$V1,
                     column="ENTREZID",
                     keytype="ENSEMBL",
                     multiVals="first")

res <- res %>% 
  dplyr::select(-V1) %>% 
  filter(symbol !="") %>% 
  mutate(newcolumn = rowMeans(.[,-c(22:24)])) %>% 
  arrange(desc(newcolumn)) %>% 
  distinct(symbol,.keep_all = T) %>% 
  dplyr::select(-newcolumn)

rownames(res) <- res[,"symbol"]
res <- res[,1:21]

library(clusterProfiler)
library(GSVA)
kegggmt <- read.gmt("ontology_Rdata/c2.cp.kegg.v7.1.symbols.gmt")
colnames(kegggmt)
kegg_list = split(kegggmt$gene, kegggmt$term)
kegg1 <- gsva(expr=as.matrix(res), kegg_list, kcdf="Gaussian",method = "gsva",parallel.sz=1)

library(GSEABase)
keggSet <- getGmt("ontology_Rdata/c2.cp.kegg.v7.1.symbols.gmt")
kegg2 <- gsva(expr=as.matrix(res), keggSet, kcdf="Gaussian",method = "gsva",parallel.sz=1)

library(pheatmap)
annotation_col = data.frame(
  Time = factor(c(rep("6H", 12),rep("20H", 9)),levels = c("6H","20H")), 
  Concentration = factor(c(rep("CDDP0", 3),rep("CDDP25", 3),rep("CDDP50", 3),rep("CDDP100", 3),
                           rep("CDDP0", 3),rep("CDDP25", 3),rep("CDDP50", 3)),levels =c("CDDP0","CDDP25","CDDP50","CDDP100"))
)
rownames(annotation_col) = colnames(kegg2)

nm <- rownames(kegg2)
nm <- gsub("KEGG_","",nm)
rownames(kegg2) <- nm
pheatmap(kegg2[1:100,],
         cluster_rows = TRUE,
         cluster_cols = FALSE,
         annotation_col =annotation_col,
         annotation_legend=TRUE, 
         show_rownames = T,
         scale = "row",
         color =colorRampPalette(c("purple", "white","orange"))(20),
         cellwidth = 10, cellheight = 5,
         fontsize = 5)
library(export)
graph2tif(file="GSVA_KEGG_CDDP25VSCON_100.tif", width=10,height=20)
save(kegg2,file="GSVA_KEGG.csv")

pheatmap(kegg2[101:186,],
         cluster_rows = TRUE,
         cluster_cols = FALSE,
         annotation_col =annotation_col,
         annotation_legend=TRUE, 
         show_rownames = T,
         scale = "row",
         color =colorRampPalette(c("purple", "white","orange"))(20),
         cellwidth = 10, cellheight = 5,
         fontsize = 5)
library(export)
graph2tif(file="GSVA_KEGG_CDDP25VSCON_101-186.tif", width=10,height=20)

Rapid_r <- read.table("data/rapid response model.txt",comment.char="!",stringsAsFactors=F,header=F)
rapid <-Rapid_r$V1
rownames(kegg2) <- nm
pheatmap(kegg2[rapid,],
         cluster_rows = TRUE,
         cluster_cols = FALSE,
         annotation_col =annotation_col,
         annotation_legend=TRUE, 
         show_rownames = T,
         scale = "row",
         color =colorRampPalette(c("purple", "white","orange"))(20),
         cellwidth = 10, cellheight = 5,
         fontsize = 5)
library(export)
graph2tif(file="GSVA_KEGG_RBE_rapid_model.tif", width=10,height=8)

Weak_r <- read.table("data/weak response model.txt",comment.char="!",stringsAsFactors=F,header=F)
weak <-Weak_r$V1
rownames(kegg2) <- nm
pheatmap(kegg2[weak,],
         cluster_rows = TRUE,
         cluster_cols = FALSE,
         annotation_col =annotation_col,
         annotation_legend=TRUE, 
         show_rownames = T,
         scale = "row",
         color =colorRampPalette(c("purple", "white","orange"))(20),
         cellwidth = 10, cellheight = 5,
         fontsize = 5)
library(export)
graph2tif(file="GSVA_KEGG_RBE_weak_model.tif", width=10,height=10)

Continous_r <- read.table("data/continuous response model.txt",comment.char="!",stringsAsFactors=F,header=F)
countinous <-Continous_r$V1
rownames(kegg2) <- nm
pheatmap(kegg2[countinous,],
         cluster_rows = TRUE,
         cluster_cols = FALSE,
         annotation_col =annotation_col,
         annotation_legend=TRUE, 
         show_rownames = T,
         scale = "row",
         color =colorRampPalette(c("purple", "white","orange"))(20),
         cellwidth = 10, cellheight = 5,
         fontsize = 5)
library(export)
graph2tif(file="GSVA_KEGG_RBE_countinous_model.tif", width=12,height=15)


####GSVA for KEGG (other cells)
rm(list = ls())
exprSet=data.table::fread("data/otherCells_CDDP counts adjust.csv",data.table = F)
res <- exprSet
library(dplyr)
library(tibble)
library(tidyr)
library(AnnotationDbi)
library(org.Hs.eg.db)
res$symbol <- res$V1
res$entrez <- mapIds(org.Hs.eg.db,
                     keys=res$V1,
                     column="ENTREZID",
                     keytype="SYMBOL",
                     multiVals="first")

res <- res %>% 
  dplyr::select(-V1) %>% 
  filter(symbol !="") %>% 
  mutate(newcolumn = rowMeans(.[,-c(46:47)])) %>% 
  arrange(desc(newcolumn)) %>% 
  distinct(symbol,.keep_all = T) %>% 
  dplyr::select(-newcolumn)

rownames(res) <- res[,"symbol"]
res <- res[,1:45]

library(clusterProfiler)
library(GSVA)
kegggmt <- read.gmt("ontology_Rdata/c2.cp.kegg.v7.1.symbols.gmt")
colnames(kegggmt)
kegg_list = split(kegggmt$gene, kegggmt$term)
kegg1 <- gsva(expr=as.matrix(res), kegg_list, kcdf="Gaussian",method = "gsva",parallel.sz=1)

library(GSEABase)
keggSet <- getGmt("ontology_Rdata/c2.cp.kegg.v7.1.symbols.gmt")
kegg2 <- gsva(expr=as.matrix(res), keggSet, kcdf="Gaussian",method = "gsva",parallel.sz=1)

colnames(res)
library(pheatmap)
annotation_col = data.frame(
  Time = factor(c(rep("0H", 3),rep("6H", 3),rep("20H", 3),rep("0H", 3),rep("6H", 3),rep("20H", 3),
                  rep("0H", 3),rep("6H", 3),rep("20H", 9),rep("0H", 3),rep("6H", 3),rep("20H", 6)),levels = c("0H","6H","20H")), 
  Concentration = factor(c(rep("CDDP0", 3),rep("CDDP25", 3),rep("CDDP10", 3),rep("CDDP0", 3),
                           rep("CDDP25", 3),rep("CDDP10", 3),rep("CDDP0", 3),rep("CDDP25", 3),rep("CDDP10", 3),
                           rep("CDDP25", 3),rep("CDDP50", 3),rep("CDDP0", 3),rep("CDDP25", 6),rep("CDDP50", 3)),levels =c("CDDP0","CDDP10","CDDP25","CDDP50"))
)
rownames(annotation_col) = colnames(kegg2)

nm <- rownames(kegg2)
nm <- gsub("KEGG_","",nm)
rownames(kegg2) <- nm
pheatmap(kegg2[1:100,],
         cluster_rows = TRUE,
         cluster_cols = FALSE,
         annotation_col =annotation_col,
         annotation_legend=TRUE, 
         show_rownames = T,
         scale = "row",
         color =colorRampPalette(c("purple", "white","orange"))(20),
         cellwidth = 10, cellheight = 5,
         fontsize = 5)
library(export)
graph2tif(file="GSVA_KEGG_othercells_1-100.tif", width=15,height=20)
save(kegg2,file="GSVA_KEGG.csv")

pheatmap(kegg2[101:186,],
         cluster_rows = TRUE,
         cluster_cols = FALSE,
         annotation_col =annotation_col,
         annotation_legend=TRUE, 
         show_rownames = T,
         scale = "row",
         color =colorRampPalette(c("purple", "white","orange"))(20),
         cellwidth = 10, cellheight = 5,
         fontsize = 5)
library(export)
graph2tif(file="GSVA_KEGG_othercells_101-186.tif", width=15,height=20)

Rapid_r <- read.table("data/rapid response model.txt",comment.char="!",stringsAsFactors=F,header=F)
rapid <-Rapid_r$V1
rownames(kegg2) <- nm
pheatmap(kegg2[rapid,],
         cluster_rows = TRUE,
         cluster_cols = FALSE,
         annotation_col =annotation_col,
         annotation_legend=TRUE, 
         show_rownames = T,
         scale = "row",
         color =colorRampPalette(c("purple", "white","orange"))(20),
         cellwidth = 10, cellheight = 5,
         fontsize = 5)
library(export)
graph2tif(file="GSVA_KEGG_othercells_rapid_model.tif", width=12,height=8)

Weak_r <- read.table("data/weak response model.txt",comment.char="!",stringsAsFactors=F,header=F)
weak <-Weak_r$V1
rownames(kegg2) <- nm
pheatmap(kegg2[weak,],
         cluster_rows = TRUE,
         cluster_cols = FALSE,
         annotation_col =annotation_col,
         annotation_legend=TRUE, 
         show_rownames = T,
         scale = "row",
         color =colorRampPalette(c("purple", "white","orange"))(20),
         cellwidth = 10, cellheight = 5,
         fontsize = 5)
library(export)
graph2tif(file="GSVA_KEGG_Othercells_weak_model.tif", width=12,height=10)

Continous_r <- read.table("data/continuous response model.txt",comment.char="!",stringsAsFactors=F,header=F)
countinous <-Continous_r$V1
rownames(kegg2) <- nm
pheatmap(kegg2[countinous,],
         cluster_rows = TRUE,
         cluster_cols = FALSE,
         annotation_col =annotation_col,
         annotation_legend=TRUE, 
         show_rownames = T,
         scale = "row",
         color =colorRampPalette(c("purple", "white","orange"))(20),
         cellwidth = 10, cellheight = 5,
         fontsize = 5)
library(export)
graph2tif(file="GSVA_KEGG_OtherCells_countinous_model.tif", width=12,height=15)

####GSVA for specific GO_CC
kegg_selected <- data.table::fread("data/GO_CC_list_CDDP_final.csv",data.table = F)

library(tidyr)
data <- kegg_selected %>% 
  pivot_longer(cols=-1,
               names_to= "fake",
               values_to = "gene") %>% 
  filter(gene !="") %>%
  as.data.frame()
colnames(data)[1] <- "term"
data <- data[,-2]

exprSet=data.table::fread("data/CDDP 20H 6H counts adjust.csv",data.table = F)
res <- exprSet
library(dplyr)
library(tibble)
library(tidyr)
library(AnnotationDbi)
library(org.Hs.eg.db)

res$symbol <- mapIds(org.Hs.eg.db,
                     keys=res$V1,
                     column="SYMBOL",
                     keytype="ENSEMBL",
                     multiVals="first")

res$entrez <- mapIds(org.Hs.eg.db,
                     keys=res$V1,
                     column="ENTREZID",
                     keytype="ENSEMBL",
                     multiVals="first")

res <- res %>% 
  dplyr::select(-V1) %>% 
  filter(symbol !="") %>% 
  mutate(newcolumn = rowMeans(.[,-c(22:24)])) %>% 
  arrange(desc(newcolumn)) %>% 
  distinct(symbol,.keep_all = T) %>% 
  dplyr::select(-newcolumn)

rownames(res) <- res[,"symbol"]
res <- res[,1:21]

library(clusterProfiler)
library(GSVA)
kegg_list = split(data$gene, data$term)
kegg1 <- gsva(expr=as.matrix(res), kegg_list, kcdf="Gaussian",method = "gsva",parallel.sz=1)
reord <- unique(data$term)
kegg_reo <- kegg1[reord,]
library(pheatmap)
annotation_col = data.frame(
  Time = factor(c(rep("6H", 12),rep("20H", 9)),levels = c("6H","20H")), 
  Concentration = factor(c(rep("CDDP0", 3),rep("CDDP25", 3),rep("CDDP50", 3),rep("CDDP100", 3),
                           rep("CDDP0", 3),rep("CDDP25", 3),rep("CDDP50", 3)),levels =c("CDDP0","CDDP25","CDDP50","CDDP100"))
)
rownames(annotation_col) = colnames(kegg_reo)

library(RColorBrewer)
pal1<-brewer.pal(12,'Paired')
pal2<-brewer.pal(8,'Set2')
mycolor<-c(pal1[8:11])

ann_colors = list(Time = c("6H" = "#79bd9a", "20H" = "#614ad3"),
                  Concentration=c(CDDP0="#FF7F00",CDDP25= "#CAB2D6",CDDP50= "#6A3D9A",
                                  CDDP100="#FFFF99"))


nm <- rownames(kegg_reo)
nm <- gsub("GOCC_","",nm)
rownames(kegg_reo) <- nm
pheatmap(kegg_reo,
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         annotation_col =annotation_col,
         annotation_colors = ann_colors,
         annotation_legend=TRUE, 
         show_rownames = T,
         scale = "row",
         color =colorRampPalette(c("#3300CC","#3399FF","white","#FF3333","#CC0000"))(20),
         cellwidth = 12, cellheight = 12,
         fontsize = 10)
library(export)
graph2tif(file="output/GSVA_GO_CC_CDDP_specific.tif")
graph2ppt(file="output/GSVA_GO_CC_CDDP_specific.pptx")
save(kegg_reo,file="output/GSVA_GO_CC_CDDP_specific.csv")

####STEM genelist
rm(list = ls())
CDDP25diff <- data.table::fread("data/Diffgene_CDDP25vscon_20h_LOGFC1.csv",data.table = F)
CDDP50diff <- data.table::fread("data/Diffgene_CDDP50vscon_20h_LOGFC1.csv",data.table = F)
diff<-rbind(CDDP25diff,CDDP50diff)
diffgene<-as.data.frame(diff$symbol)
colnames(diffgene)[1]<-"symbol"
diffgene<-distinct(diffgene,symbol,.keep_all=TRUE)
write.csv(diffgene,file = "data/diffgene_CDDP25_50_20HforSTEM.csv")
library("sva")
exprSet <- data.table::fread("data/count_table_CDDP_20H.csv",data.table = F)
exprSet <- distinct(exprSet,GeneSymbol,.keep_all=TRUE)
exprSet <- filter(exprSet,GeneSymbol !="")
rownames(exprSet) <- exprSet[,11]
exprSet <- exprSet[,2:10]
lograwcount1<- log2(exprSet+1)  
STEMexpr <- lograwcount1[diffgene$symbol,]
STEMexpr<- filter(STEMexpr,!is.na(L.CON.1))
write.csv(STEMexpr,file = "data/STEM_gene_expression_log2.csv")
  
  