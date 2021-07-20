library(WGCNA)
library(reshape2)
library(stringr)
options(stringsAsFactors = FALSE)
enableWGCNAThreads()
#WGCNA
rm(list = ls())
options(stringsAsFactors = F)
a = a <- read.csv("E:/R/workspace/WGCNA-JJ.csv", header=T)
datTraits=read.csv("E:/R/workspace/datTraits-JJ.csv", header=T)
row.names(a) <- a$ID
a <- a[,-1]


##Screen the top 75% of the genes with the median absolute deviation > 0.01
m.vars=apply(a,1,var)
datExpr0=a[which(m.vars>quantile(m.vars, probs = seq(0, 1, 0.25))[4]),]

m.mad <- apply(a,1,mad)
dataExprVar <- a[which(m.mad > 
                         max(quantile(m.mad, probs=seq(0, 1, 0.25))[2],0.01)),]
WGCNA_matrix = t(datExpr0)
dataExpr0 <- WGCNA_matrix
dataExpr <- dataExpr0

## Detect missing values
gsg = goodSamplesGenes(dataExpr, verbose = 3)
if (!gsg$allOK){
  if (sum(!gsg$goodGenes)>0) 
    printFlush(paste("Removing genes:", 
                     paste(names(dataExpr)[!gsg$goodGenes], collapse = ",")));
  if (sum(!gsg$goodSamples)>0) 
    printFlush(paste("Removing samples:", 
                     paste(rownames(dataExpr)[!gsg$goodSamples], collapse = ",")));
  dataExpr = dataExpr[gsg$goodSamples, gsg$goodGenes]
}

nGenes = ncol(dataExpr)
nSamples = nrow(dataExpr)
dim(dataExpr)
head(dataExpr)[,1:8]

#calculate β value
powers = c(c(1:10), seq(from = 12, to=30, by=2))
sft = pickSoftThreshold(dataExpr, powerVector = powers, verbose = 5)
par(mfrow = c(1,2))
cex1 = 0.9
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"))
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red")
abline(h=0.8,col="red")

plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
sft
#so β=16

#build a co-expression matrix
corType = "pearson"
exprMat <- "E:/R/results/ICVClean.txt" 
type="unsigned"
net = blockwiseModules(dataExpr, power = 12, maxBlockSize = 12000,
                       TOMType = type, minModuleSize = 100,
                       reassignThreshold = 0, mergeCutHeight = 0.25,
                       numericLabels = TRUE, pamRespectsDendro = FALSE,
                       saveTOMs=TRUE, corType = corType, 
                       maxPOutliers=0, loadTOMs=TRUE,
                       saveTOMFileBase = paste0(exprMat, ".tom"),
                       verbose = 3)
table(net$colors)
moduleLabels = net$colors
moduleColors = labels2colors(moduleLabels)
table(moduleColors)
plotDendroAndColors(net$dendrograms[[1]], moduleColors[net$blockGenes[[1]]],
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
MEs = net$MEs
MEs_col = MEs
colnames(MEs_col) = paste0("ME", labels2colors(
  as.numeric(str_replace_all(colnames(MEs),"ME",""))))
MEs_col = orderMEs(MEs_col)
plotEigengeneNetworks(MEs_col, "Eigengene adjacency heatmap", 
                      marDendro = c(3,4,2,4),
                      marHeatmap = c(3,4,2,2), plotDendrograms = T, 
                      xLabelsAngle = 90)

#co-traits association analysis
design=model.matrix(~0+ traitData$tissue)
colnames(design)=levels(as.factor(traitData$tissue))
moduleColors <- labels2colors(net$colors)
MEs0 = moduleEigengenes(dataExpr, moduleColors)$eigengenes
MEs = orderMEs(MEs0)
moduleTraitCor = cor(MEs, design , use = "p");
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples)
sizeGrWindow(10,6)
textMatrix = paste(signif(moduleTraitCor, 2), "\n(",
                   signif(moduleTraitPvalue, 1), ")", sep = "")
dim(textMatrix) = dim(moduleTraitCor)
par(mar = c(6, 8.5, 3, 3))
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = colnames(design),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = greenWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.5,
               zlim = c(-1,1),
               main = paste("Module-trait relationships"))


##Visual network
load(net$TOMFiles[1], verbose=T)
TOM <- as.matrix(TOM)
dissTOM = 1-TOM
plotTOM = dissTOM^7
diag(plotTOM) = NA
TOMplot(plotTOM, net$dendrograms, moduleColors, 
        main = "Network heatmap plot, all genes")

#Randomly select 400 genes to perform heatmap
nSelect = 400
set.seed(10)
select = sample(nGenes, size = nSelect)
selectTOM = dissTOM[select, select]
selectTree = hclust(as.dist(selectTOM), method = "average")
selectColors = moduleColors[select]
sizeGrWindow(9,9)
plotDiss = selectTOM^7;
diag(plotDiss) = NA;
TOMplot(plotDiss, selectTree, selectColors, main = "Network heatmap plot, selected genes")

#extract genes from candidate module
module='pink'
probes=colnames(datExpr)
inModule=(moduleColors==module)
modProbes=probes[inModule]
head(modProbes)
write.csv(modProbes, file = 'E:/R/results/module-gene/case-control/ICV/pink.csv')

#GO/KEGG analysis
setwd('E:/R/results/module-gene/case-control/ICV/')
library(clusterProfiler)
library(org.Bt.eg.db)

rt <- read.csv("pink.csv",header = T)

eg <- bitr(rt$id,fromType = "SYMBOL",toType = "ENTREZID",OrgDb = org.Bt.eg.db)  
dim(eg)  

ego <- enrichGO(eg$ENTREZID,org.Bt.eg.db,ont = "ALL",pvalueCutoff = 0.05,qvalueCutoff = 0.05,readable = TRUE)  #GO
write.csv(ego,"E:/R/results/WGCNA/pink-GO.csv",row.names = F) 

kk <- enrichKEGG(eg$ENTREZID,pvalueCutoff = 0.05,qvalueCutoff = 0.05,organism = "bta") #kegg
write.csv(kk,"E:/R/results/WGCNA/pink-KEGG.csv",row.names = F) 
