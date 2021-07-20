# Display the current working directory
rm(list = ls())
options(stringsAsFactors = F)
getwd();
# If necessary, change the path below to the directory where the data files are stored. 
# "." means current directory. On Windows use a forward slash / instead of the usual \.
workingDir = "E:/R/workspace/WGCNA/consensus module/ICV";
setwd(workingDir); 

# Load the WGCNA and stringr package
library(WGCNA);
library(stringr);
# The following setting is important, do not omit.
options(stringsAsFactors = FALSE);
#Read in the female liver data 
case= read.csv("ICV-case.csv");
# Take a quick look at what is in the data set:
dim(case);
names(case);

#convert data
datExpr0 = as.data.frame(t(case[, c(2:7)]));
names(datExpr0) = case$ID;
rownames(datExpr0) = names(case)[c(2:7)];
dim(datExpr0)

#Check that all genes and samples have sufficiently low numbers of missing values.
gsg = goodSamplesGenes(datExpr0, verbose = 3);
gsg$allOK

if (!gsg$allOK)
{
  # Optionally, print the gene and sample names that were removed:
  if (sum(!gsg$goodGenes)>0) 
    printFlush(paste("Removing genes:", paste(names(datExpr0)[!gsg$goodGenes], collapse = ", ")));
  if (sum(!gsg$goodSamples)>0) 
    printFlush(paste("Removing samples:", paste(rownames(datExpr0)[!gsg$goodSamples], collapse = ", ")));
  # Remove the offending genes and samples from the data:
  datExpr0 = datExpr0[gsg$goodSamples, gsg$goodGenes]
}

sampleTree = hclust(dist(datExpr0), method = "average");
# Plot the sample tree: Open a graphic output window of size 12 by 9 inches
# The user should change the dimensions if the window is too large or too small.
sizeGrWindow(12,9)

#Sample clustering to detect outliers
par(cex = 0.6);
par(mar = c(0,4,2,0))
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5, 
     cex.axis = 1.5, cex.main = 2)

# Plot a line to show the cut
abline(h = 15, col = "red");
# Determine cluster under the line
clust = cutreeStatic(sampleTree, cutHeight = 15, minSize = 10)
table(clust)
# clust 1 contains the samples we want to keep.
keepSamples = (clust==1)
datExpr = datExpr0[keepSamples, ]
datExpr = datExpr0
nGenes = ncol(datExpr)
nSamples = nrow(datExpr)

#Read the traits file
traitData = read.csv("ClinicalTraits-case.csv");
dim(traitData)
names(traitData)

# remove columns that hold information we do not need.
allTraits = traitData
dim(allTraits)
names(allTraits)

# Form a data frame analogous to expression data that will hold the clinical traits.

caseSamples = rownames(datExpr);
traitRows = match(caseSamples, allTraits$ID);
dataTraits = allTraits[traitRows, -1];
rownames(dataTraits) = allTraits[traitRows, 1];

collectGarbage();

# Re-cluster samples
sampleTree2 = hclust(dist(datExpr), method = "average")
# Convert traits to a color representation: white means low, red means high, grey means missing entry
traitColors = numbers2colors(datExpr, signed = FALSE);
# Plot the sample dendrogram and the colors underneath.
plotDendroAndColors(sampleTree2, traitColors,
                    groupLabels = names(dataTraits), 
                    main = "Sample dendrogram and trait heatmap")

# Choose a set of soft-thresholding powers
powers = c(c(1:10), seq(from = 12, to=20, by=2))
# Call the network topology analysis function
sft = pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)
sft
# Plot the results:
sizeGrWindow(9, 5)
par(mfrow = c(1,2));
cex1 = 0.9;
# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));


text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");
# this line corresponds to using an R^2 cut-off of h
abline(h=0.80,col="red")
# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")

#Build a co-expression matrix
net = blockwiseModules(datExpr, power = 9,
                       TOMType = "unsigned", minModuleSize = 100,
                       reassignThreshold = 0, mergeCutHeight = 0.25,
                       numericLabels = TRUE, pamRespectsDendro = FALSE,
                       saveTOMs = TRUE,
                       saveTOMFileBase = "CASETOM", 
                       verbose = 3)

# open a graphics window
sizeGrWindow(12, 9)
# Convert labels to colors for plotting
mergedColors = labels2colors(net$colors)
# Plot the dendrogram and the module colors underneath
plotDendroAndColors(net$dendrograms[[1]], mergedColors[net$blockGenes[[1]]],
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
moduleLabels = net$colors
moduleColors = labels2colors(net$colors)
MEs = net$MEs;
geneTree = net$dendrograms[[1]];
save(MEs, moduleLabels, moduleColors, geneTree, 
     file = "CASE-02-networkConstruction-auto.RData")

#Same-character association analysis
design=model.matrix(~0+ traitData$tissue)
colnames(design)=levels(as.factor(traitData$tissue))
moduleColors <- labels2colors(net$colors)
MEs0 = moduleEigengenes(datExpr, moduleColors)$eigengenes
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

#Inter-module correlation
MEs = net$MEs
MEs_col = MEs
colnames(MEs_col) = paste0("ME", labels2colors(
  as.numeric(str_replace_all(colnames(MEs),"ME",""))))
MEs_col = orderMEs(MEs_col)
plotEigengeneNetworks(MEs_col, "Eigengene adjacency heatmap", 
                      marDendro = c(3,4,2,4),
                      marHeatmap = c(3,4,2,2), plotDendrograms = T, 
                      xLabelsAngle = 90)

#save the case-analysis result
lnames = load(file = "CASE-02-networkConstruction-auto.RData");
lnames

#extract genes
module='black'
probes=colnames(datExpr)
inModule=(moduleColors==module)
modProbes=probes[inModule]
head(modProbes)
write.csv(modProbes, file = 'E:/R/results/module-gene/case-control/JJ/black.csv')

#GO/KEGG analysis
setwd('E:/R/results/module-gene/case-control/JJ/')
library(clusterProfiler)
library(org.Bt.eg.db)

rt <- read.csv("black.csv",header = T)

eg <- bitr(rt$id,fromType = "SYMBOL",toType = "ENTREZID",OrgDb = org.Bt.eg.db)  
dim(eg)  

ego <- enrichGO(eg$ENTREZID,org.Bt.eg.db,ont = "ALL",pvalueCutoff = 0.05,qvalueCutoff = 0.05,readable = TRUE)  #GO
write.csv(ego,"E:/R/results/WGCNA/JJ-GO.csv",row.names = F) 

kk <- enrichKEGG(eg$ENTREZID,pvalueCutoff = 0.05,qvalueCutoff = 0.05,organism = "bta")
write.csv(kk,"E:/R/results/WGCNA/JJ-KEGG.csv",row.names = F) 

