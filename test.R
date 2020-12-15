library( "DESeq2" )
library('org.Hs.eg.db')
library('pheatmap')
library('RColorBrewer')
library('BiocManager')
library(ggplot2)

contpluse2aData <- read.table(file="Galaxy51-[featureCounts_on_data_46__Counts].tabular", sep = "\t", header = T)
colnames(contpluse2aData)[1:2] <- c('Gene_name', 'contpluse2a')

contpluse2bData <- read.table(file="Galaxy53-[featureCounts_on_data_47__Counts].tabular", sep = "\t", header = T)
colnames(contpluse2bData)[1:2] <- c('Gene_name', 'contpluse2b')

cont.e2aData <- read.table(file="Galaxy55-[featureCounts_on_data_48__Counts].tabular", sep = "\t", header = T)
colnames(cont.e2aData)[1:2] <- c('Gene_name', 'cont.e2a')

cont.e2bData <- read.table(file="Galaxy57-[featureCounts_on_data_49__Counts].tabular", sep = "\t", header = T)
colnames(cont.e2bData)[1:2] <- c('Gene_name', 'cont.e2b')

countData <- Reduce(function(x, y) merge(x, y, by = 1, all=TRUE), list(contpluse2aData, contpluse2bData, cont.e2aData, cont.e2bData))

mapping <-select(org.Hs.eg.db, as.character(countData$Gene_name), keytype = "ENTREZID", column="SYMBOL")
countData <- merge(countData, mapping, by.x = 'Gene_name', by.y = "ENTREZID", all = 'TRUE')
missing <- is.na(countData$SYMBOL)
countData <- countData[!missing,]

o <- order(rowSums(countData[,c(2:5)]), decreasing=TRUE)
countData <- countData[o,]

row.names(countData) <- countData$SYMBOL
countData$SYMBOL <- NULL
countData$Gene_name <- NULL

metaData <- read.table(file="metadata.tsv", sep = "\t", header = T)

#Create DESeq2Dataset object
dds <- DESeqDataSetFromMatrix(countData = countData, colData = metaData, design = ~ 1)
#Normalization
dds <- estimateSizeFactors(dds)
#normalized_counts <- counts(dds, normalized=TRUE)

logcounts <- log2( counts(dds, normalized=TRUE) + 1 )
pc <- prcomp( t( logcounts ) )

design(dds) <- ~ replicate + treatment

dds <- DESeq( dds )
res <- results( dds )

head(results(dds, contrast=c("treatment","Veh","E2")))

resBigFC <- results(dds, lfcThreshold=1, altHypothesis="greaterAbs")
plotMA(resBigFC, ylim=c(-5,5))
abline(h=c(-1,1),lwd=5)

#Top genes, sorted by pvalue
resSort <- res[order(res$pvalue),]
head(resSort)

#Significant genes
resSigUp = res[ which(res$padj < 0.05 & res$log2FoldChange > 0.58), ]
resSigDown = res[ which(res$padj < 0.05 & res$log2FoldChange < -0.58), ]
resSig = rbind(resSigUp, resSigDown)

sink("resSig.txt")
cat(sapply(rownames(resSig), toString), sep="\n")
sink()


rld <- rlog(dds, blind=FALSE)

#heatmap
sampleDists <- dist(t(assay(rld)))
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(rld$sample, rld$treatment, rld$replicate, sep="-" )
colnames(sampleDistMatrix) <- paste(rld$sample, rld$treatment, rld$replicate, sep="-" )
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix, clustering_distance_rows=sampleDists, clustering_distance_cols=sampleDists, col=colors)

#PCA
plotPCA(rld, intgroup = c("treatment", "replicate"))

#Gene heatmap
geneVars <- rowVars(assay(rld))
geneVarsOrdered <- order(geneVars, decreasing = TRUE)
topVarGenes <- head(geneVarsOrdered, 50)
mat <- assay(rld)[ topVarGenes, ]
mat <- mat - rowMeans(mat)
df <- as.data.frame(colData(rld)[,c("treatment","replicate")])
clear_col_names <- paste( rld$sample, rld$treatment, rld$replicate, sep=".")
topGenesHeatmap <- pheatmap(mat, annotation_col=df, labels_col = clear_col_names)
