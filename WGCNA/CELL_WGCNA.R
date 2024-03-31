rm(list = ls())
library(WGCNA)
raw_counts <- data.table::fread('data/count_table_CDDP_ALL2.csv', data.table = F)
rownames(raw_counts) <- raw_counts[,1]
raw_counts <- raw_counts[,-1]
raw_counts<-raw_counts[,1:21]
datExpr0 <- as.data.frame(t(raw_counts))
### 1.goodSamplesGenes 
gsg = goodSamplesGenes(datExpr0, verbose = 3)
gsg$allOK
if (!gsg$allOK){
  if (sum(!gsg$goodGenes)>0) 
      if (sum(!gsg$goodSamples)>0) 
      printFlush(paste("Removing samples:", paste(rownames(datExpr0)[!gsg$goodSamples], collapse = ", ")));
    datExpr0 = datExpr0[gsg$goodSamples, gsg$goodGenes]
}
raw_counts <- as.data.frame(t(datExpr0))
library(RColorBrewer)
pheno <- data.table::fread('data/sample_metadata_CDDP_ALL2.csv',data.table = F)
num_conditions <- nlevels(as.factor(pheno$condition))
pal <- colorRampPalette(brewer.pal(num_conditions, "Set1"))(num_conditions)
cond_colors <- pal[as.integer(as.factor(pheno$condition))]

library(gplots)
heatmap.2(cor(raw_counts), 
          RowSideColors=cond_colors,
          trace='none', 
          main='Sample correlations (raw)')
### 2.filtering genes
low_count_mask <- rowSums(raw_counts) < 100*ncol(raw_counts)
sum(low_count_mask)
raw_counts <- raw_counts[!low_count_mask,]
log_counts <- log2(raw_counts  + 1)
log_counts <- log_counts[apply(log_counts, 1, var) > 0,]

rownames(pheno) <- pheno$sample_id
index <- intersect(pheno$sample_id,colnames(log_counts))
pheno <- pheno[index,]
datExpr = as.data.frame(t(log_counts[,index]))

sampleTree2 = hclust(dist(datExpr), method = "average")
traitColors = numbers2colors(as.numeric(as.factor(pheno$condition)), signed = FALSE)
plotDendroAndColors(sampleTree2, 
                    traitColors,
                    groupLabels = names(pheno)[2], 
                    main = "Sample dendrogram and trait heatmap")

### 3.batch effect
color <- as.numeric(as.factor(pheno$condition))
datExprt <-t(datExpr)
boxplot(datExprt,outline=FALSE, notch=T, las=2,col=color)

library(limma) 
datExprt=normalizeBetweenArrays(datExprt)
par(mar = rep(0, 4))
boxplot(datExprt,outline=FALSE, notch=T, las=2,col=color)

library(RColorBrewer)
num_conditions <- nlevels(as.factor(pheno$condition))
pal <- colorRampPalette(brewer.pal(num_conditions, "Set1"))(num_conditions)
cond_colors <- pal[as.integer(as.factor(pheno$condition))]

library(gplots)
par(mar = rep(4, 4))
heatmap.2(cor(datExprt), 
          RowSideColors=cond_colors,
          trace='none', 
          main='Sample correlations (batch effect)',
          margins = c(7, 7))

datExpr <- t(datExprt)
#save(file = "remove_batch_CDDP_ALL_strict3.2.Rdata")
### 
load(file = "remove_batch_CDDP_ALL_strict3.2.Rdata")
dim(datExpr)
collectGarbage()
### 设置多线程
enableWGCNAThreads()
powers = c(c(1:10), seq(from = 12, to=40, by=2))
length(powers)
sft = pickSoftThreshold(datExpr, powerVector = powers,networkType="signed", verbose = 5)
test <- sft$fitIndices
sizeGrWindow(9, 5)
par(mfrow = c(1,2))
cex1 = 0.9
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red")

abline(h=0.8,col="red")
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
abline(h=100,col="red")
library(export)
graph2tif(file="output/scalepower_strict3.2.tiff")

cor <- WGCNA::cor
net = blockwiseModules(datExpr, power = 26,
                       TOMType = "signed", minModuleSize = 30,
                       reassignThreshold = 0, mergeCutHeight = 0.25,
                       numericLabels = TRUE, pamRespectsDendro = FALSE,
                       saveTOMs = TRUE,
                       saveTOMFileBase = "myTOM",
                       verbose = 3)
cor<-stats::cor
table(net$colors)

sizeGrWindow(12, 9)
mergedColors = labels2colors(net$colors)
table(mergedColors)
# Plot the dendrogram and the module colors underneath
plot(net$dendrograms[[1]])
plotDendroAndColors(net$dendrograms[[1]], 
                    mergedColors[net$blockGenes[[1]]],
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
library(export)
graph2ppt(file="output/geneModule_vst.pptx")
graph2tif(file="output/geneModule_vst.tiff")

# Recalculate module eigengenes
moduleColors = labels2colors(net$colors)
MEs = moduleEigengenes(datExpr, moduleColors)$eigengenes
# Add the weight to existing module eigengenes
MET = orderMEs(MEs)
# Plot the relationships among the eigengenes and the trait
sizeGrWindow(5,7.5)
par(cex = 0.9)
plotEigengeneNetworks(MET, "", marDendro = c(0,4,1,2), marHeatmap = c(3,4,1,2), cex.lab = 0.8, xLabelsAngle= 90)
graph2ppt(file="output/ModuleModule_correlation_vst.pptx")
graph2tif(file="output/ModuleModule_correlation_vst.tiff")

moduleLabels = net$colors
moduleColors = labels2colors(net$colors)
geneTree = net$dendrograms[[1]]
#save(moduleLabels, moduleColors, geneTree, file = "module_detection_CDDP_ALL_strict3.2.RData")

###Module-trait relationships
rm(list = ls())
load(file = "remove_batch_CDDP_ALL_strict3.2.Rdata")
load(file = "module_detection_CDDP_ALL_strict3.2.RData")
module_eigengenes <- MEs
test <- moduleEigengenes(datExpr, moduleColors)
MEs0 = moduleEigengenes(datExpr, moduleColors)$eigengenes
MEs = orderMEs(MEs0)
write.csv(MEs,file="output/ME_across_samples.csv")
datTraits <- pheno
datTraits$condition <- factor(datTraits$condition,levels = c("CON_6H","CDDP25_6H","CDDP50_6H","CDDP100_6H",
                                                             "CON_20H","CDDP25_20H","CDDP50_20H"))


design=model.matrix(~0+ datTraits$condition)
design <- as.data.frame(design)
colnames(design)=levels(datTraits$condition)
design
design$Concentration <- as.numeric(datTraits$CDDP_concentration)
design$Time <- as.numeric(datTraits$CDDP_time)

moduleTraitCor = cor(MEs, design, use = "p")
nSamples <- nrow(datExpr)
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples)

sizeGrWindow(10,6)
textMatrix =  paste(signif(moduleTraitCor, 2), "\n(",
                    signif(moduleTraitPvalue, 1), ")", sep = "");
dim(textMatrix) = dim(moduleTraitCor)
par(mar = c(6, 8.5, 3, 3));
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = colnames(design),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.5,
               zlim = c(-1,1),
               main = paste("Module-trait relationships"))
graph2ppt(file="output/ModuleTrait_labledHeatmap2_strict3.2.pptx")
graph2tif(file="output/ModuleTrait_labledHeatmap2_strict3.2.tiff")

### identify genes invloved in different modules
colorlevels=unique(moduleColors)
for (i in c(1:length(colorlevels))) 
{
  print(i)
  module=colorlevels[[i]]
  moduleGenes = moduleColors==module
  gene=data.frame(ID=colnames(datExpr)[moduleGenes])
  gene <- as.character(gene$ID)
  library(clusterProfiler)
  gene = bitr(gene, fromType="ENSEMBL", toType=c("ENTREZID","SYMBOL"), OrgDb="org.Hs.eg.db")
  moduleName=paste0("ME",module,"List_strict3.2.csv")
  write.csv(gene,file = moduleName)
}

module = "turquoise"###also other module genes
moduleGenes = moduleColors==module
datExpr <- as.data.frame(datExpr)
gene=data.frame(ID=colnames(datExpr)[moduleGenes])
gene <- as.character(gene$ID)
library(clusterProfiler)
gene = bitr(gene, fromType="ENSEMBL", toType=c("ENTREZID","SYMBOL"), OrgDb="org.Hs.eg.db")
head(gene)
dd <- datExpr[,moduleGenes]
library(dplyr)
library(tibble)
library(tidyr)
test <- dd[1:10,1:10]
dd <- dd %>% 
  rownames_to_column("sample") %>% 
  mutate(group = datTraits$condition) %>% 
  dplyr::select(group,everything()) %>% 
  pivot_longer(cols = 3:ncol(.),
               names_to = "gene",
               values_to = "expression") 
library(ggplot2)
library(RColorBrewer)
ggplot(dd,aes(x = group,y = expression)) +
  geom_boxplot(aes(fill=group))+
  scale_fill_brewer(palette="YlGnBu")
library(export)
graph2ppt(file="output/MEturquoise_geneexpression.pptx")
graph2tif(file="output/MEturquoise_geneexpression.tiff")


####Heatmap using Module eigengene value of different samples
hetadata=MEs[,1:9]
colnames(MEs)
annotation_col = data.frame(
  Module=factor(c("green","magenta","turquoise","brown","blue","yellow","red","black","pink")))
rownames(annotation_col) <- colnames(MEs)[1:9]
ann_colors = list(
  Module = c(green="green",magenta="magenta",turquoise="turquoise",brown="brown",blue="blue",
             yellow="yellow",red="red",black="black",pink="pink")
)

head(ann_colors)

color.key<-c("#3300CC","#3399FF","white","#FF3333","#CC0000")
library(pheatmap)
pheatmap(hetadata, 
         cluster_rows = F,
         cluster_cols = F,
         annotation_col =annotation_col, 
         annotation_legend=TRUE, 
         show_rownames = T,
         show_colnames = F,
         scale = "row", 
         color =colorRampPalette(color.key)(50),
         annotation_colors = ann_colors,
         cellwidth = 20, cellheight =12,
         fontsize = 10)

library(export)
graph2ppt(file="output/Modulesample_hetamap_strict3.pptx")
graph2tif(file="output/Modulesample_hetamap_strict3.tiff")


####egigene expression across different groups
library(ggplot2)
library(ggpubr)
library(export)
colorlevels=unique(moduleColors)
for (i in c(1:length(colorlevels))) 
{
  module=colorlevels[[i]]
  dd <- data.frame(group= datTraits$condition,
                   module=MEs[,paste0("ME",module)])
  p1=ggboxplot(dd, x = "group", y = "module",
               fill = module, add="jitter", palette = "jama",lwd=0.8)+
    stat_compare_means(aes(label = paste0("p = ", ..p.format..)),label.y = 0.4, size = 6)+
    theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1,size = 12))+
    theme(axis.text.y = element_text(size = 12))+
    labs(x=NULL,y=NULL)
  filename1=paste0("output2/Egiengene_",module,"_strict3.pptx")
  filename2=paste0("output2/Egiengene_",module,"_strict3.tiff")
  graph2tif(p1,file=filename2,width=6,height=6)
  graph2ppt(p1,file=filename1,width=6,height=6)
}



