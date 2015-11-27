source("https://bioconductor.org/biocLite.R")
biocLite("EDASeq")
library(EDASeq)

# fetch GC content and length of longest transcript for each gene

library("biomaRt")
mart <- useMart(biomart="ENSEMBL_MART_ENSEMBL", host="grch37.ensembl.org", path="/biomart/martservice" ,dataset="hsapiens_gene_ensembl") # GRCh37, v75
gclength <- getBM(attributes=c("ensembl_gene_id", "percentage_gc_content", "transcript_length"), mart=mart)
gclength <- gclength[order(gclength$ensembl_gene_id, -gclength$transcript_length),]
gclength <- gclength[!duplicated(gclength$ensembl_gene_id),]
gclength <- data.frame(row.names = gclength$ensembl_gene_id, gc=gclength$percentage_gc_content, length=gclength$transcript_length)

# set up sample table

# iAMP
samples <- read.delim("/mnt/projects/iamp/data/qlucore/annotations.txt")
count.files <- read.delim("/mnt/projects/iamp/results/anduril/execute/_deseq_TumorVsNormal_countFiles_array1/array/_index")
sampleTable <- merge(count.files, samples[,c("Name", "Subtype")], by.x="Key", by.y="Name", all.x = T)
names(sampleTable) <- c("Key", "File", "group")

# Fiona
samples <- read.delim("/mnt/projects/fiona/data/samples.csv")
count.files <- read.delim("/mnt/projects/fiona/results/anduril/execute/_deseq_oeERvsEmpty_counts_array1/array/_index")
sampleTable <- merge(count.files, samples[,c("Alias", "Batch")], by.x="Key", by.y="Alias", all.x = T)
sampleTable$group <- as.factor(gsub("oe(ER|Em|RHD).*", "\\1", sampleTable$Key))

# IKAROS
sampleTable <- read.delim("/mnt/projects/ikaros/data/samples.csv")
rownames(sampleTable) <- sampleTable$Alias
counts <- read.delim("/mnt/projects/ikaros/results/anduril/execute/htseqExprMatrix/countArray/all.csv", row.names = 1)
sampleTable <- sampleTable[sampleTable$Alias %in% colnames(counts),]
counts <- counts[,rownames(sampleTable)]
sampleTable$group <- as.factor(sampleTable$Subtype)
colors <- list(IK6="red", IKD="blue", IKC="black")

# get count matrix via DESeq

library(DESeq2)
deseq <- DESeqDataSetFromHTSeqCount(sampleTable = sampleTable, directory = "/", design = ~group)
counts <- counts(deseq)

# EDAseq

expressed <- apply(counts, 1, function(x) sum(x) > 10)
common <- intersect(rownames(counts[expressed,]), rownames(gclength))
data <- newSeqExpressionSet(
  counts=as.matrix(counts[common,]), 
  featureData = gclength[common,], 
  phenoData = data.frame(conditions=sampleTable$group, color=unlist(colors[sampleTable$Subtype]), row.names = colnames(counts))
)
data

boxplot(data, main="Expression per sample", las=2, cex.axis=0.7)
MDPlot(data, c(1,3), main="MA plot")
meanVarPlot(data, log=TRUE, main="Mean-variance relationship")

biasPlot(data, "gc", log=TRUE, ylim=c(0,8), col=c("red", "blue", "darkgray"), lwd=1.5, main="GC bias per sample (before norm.)")
legend("topleft", levels(data@phenoData@data$conditions), fill=c("red", "blue", "darkgray"))
#biasPlot(data, "length", log=TRUE, ylim=c(1,9), main="Length bias per sample")

pdf("/mnt/projects/ikaros/results/edaseq.pdf", height=10, width=15)
par(mfrow=c(1,2))
biasPlot(data[,c(1,2)], "gc", log=TRUE, ylim=c(1,7), main="GC bias per sample")
lfc <- log(counts(data)[,1]+0.1) - log(counts(data)[,2]+0.1)
biasBoxplot(lfc, fData(data)$gc, main="log2FC by gc", las=2)
dev.off()

# normalization

dataWithin <- withinLaneNormalization(data, "gc", which="full")
dataNorm <- betweenLaneNormalization(dataWithin, which="full")

biasPlot(dataNorm, "gc", log=TRUE, ylim=c(1,7), main="GC bias per sample (after norm.)")

dataOffset <- withinLaneNormalization(data, "gc", which="full", offset=TRUE)
dataOffset <- betweenLaneNormalization(dataOffset, which="full", offset=TRUE)

normFactors <- exp(-1 * dataOffset@assayData$offset)
normFactors <- normFactors / exp(rowMeans(log(normFactors)))

# DESeq using raw counts
#----

deseq <- DESeqDataSetFromMatrix(
  countData = dataNorm@assayData$counts, 
  colData = data.frame(row.names=sampleTable$Key, group=sampleTable$group), 
  design = ~group
)
deseq <- estimateSizeFactors(deseq)
deseq <- estimateDispersions(deseq)
deseq <- nbinomWaldTest(deseq)
#res <- as.data.frame(results(deseq, contrast=c("group", "IK6", "IKC"), cooksCutoff=FALSE)) # IKAROS
res <- as.data.frame(results(deseq, contrast=c("group", "iAMP", "ER"), cooksCutoff=FALSE)) # iAMP
res <- res[order(res$padj),]
res <- merge(res, gclength, by="row.names")
rownames(res) <- res$Row.names ; res <- res[,-1]
res$gcbin <- cut(res$gc, seq(0, 100, by=10))
res$lengthbin <- cut(res$length, c(100,200,500,1000,5000,10000,20000))
nrow(res[!is.na(res$padj) & res$padj <= 1e-2 & res$pvalue <= 1e-5,])

# DESeq using normalized counts
#----

deseq.norm <- DESeqDataSetFromMatrix(
  countData = dataNorm@assayData$normalizedCounts, 
  colData = data.frame(row.names=sampleTable$Key, group=sampleTable$group), 
  design = ~group
)
sizeFactors(deseq.norm) <- rep(1, nrow(sampleTable))
deseq.norm <- estimateDispersions(deseq.norm)
deseq.norm <- nbinomWaldTest(deseq.norm)
#res.norm <- as.data.frame(results(deseq.norm, contrast=c("group", "IK6", "IKC"), cooksCutoff=FALSE)) # IKAROS
res.norm <- as.data.frame(results(deseq.norm, contrast=c("group", "iAMP", "ER"), cooksCutoff=FALSE)) # iAMP
res.norm <- res.norm[order(res.norm$padj),]
res.norm <- merge(res.norm, gclength, by="row.names")
rownames(res.norm) <- res.norm$Row.names ; res.norm <- res.norm[,-1]
res.norm$gcbin <- cut(res.norm$gc, seq(0, 100, by=10))
res.norm$lengthbin <- cut(res.norm$length, c(100,200,500,1000,5000,10000,20000))
nrow(res.norm[!is.na(res.norm$padj) & res.norm$padj <= 1e-2 & res.norm$pvalue <= 1e-5,])

# log2fc by GC bin

pdf("/mnt/projects/ikaros/results/edaseq.fcByGCbin.pdf", width=10, height=15)
par(mfrow=c(2,1))
boxplot(log2FoldChange~gcbin, res, ylim=c(-5, 5), main="IK6 vs. IKC DESeq2 FC w/o EDAseq GC norm")
abline(h=0, col="lightgray", lty=2)
boxplot(log2FoldChange~gcbin, res.norm, ylim=c(-5, 5), main="IK6 vs. IKC DESeq2 FC with EDAseq GC norm")
abline(h=0, col="lightgray", lty=2)
dev.off()

# log2fc by gene length

pdf("/mnt/projects/iamp/results/edaseq.fcByGeneLength.pdf", width=10, height=15)
par(mfrow=c(2,1))
boxplot(abs(log2FoldChange)~lengthbin, res, ylim=c(0, 5), ylab="abs(log2(fc))", main="iAMP vs. ER DESeq2 FC w/o EDAseq GC norm")
abline(h=0.5, col="gray", lty=2)
boxplot(abs(log2FoldChange)~lengthbin, res.norm, ylim=c(0, 5), ylab="abs(log2(fc))", main="iAMP vs. ER DESeq2 FC with EDAseq GC norm")
abline(h=0.5, col="gray", lty=2)
dev.off()


# p-value by GC bin

pdf("/mnt/projects/iamp/results/edaseq.pvalueByGCbin.pdf", width=10, height=15)
par(mfrow=c(2,1))
boxplot(log(pvalue, 10)~gcbin, res, ylim=c(0, -10), main="iAMP vs. ER DESeq2 p-values w/o EDAseq GC norm")
abline(h=-5, col="gray", lty=2)
boxplot(log(pvalue, 10)~gcbin, res.norm, ylim=c(0, -10), main="iAMP vs. ER DESeq2 p-values with EDAseq GC norm")
abline(h=-5, col="gray", lty=2)
dev.off()

# p-value by gene length

pdf("/mnt/projects/iamp/results/edaseq.pvalueByGeneLength.pdf", width=10, height=15)
par(mfrow=c(2,1))
boxplot(log(pvalue, 10)~lengthbin, res, ylim=c(0, -4), ylab="log10(p-value)", main="iAMP vs. ER DESeq2 p-values w/o EDAseq GC norm")
abline(h=-0.8, col="gray", lty=2)
boxplot(log(pvalue, 10)~lengthbin, res.norm, ylim=c(0, -4), ylab="log10(p-value)", main="iAMP vs. ER DESeq2 p-values with EDAseq GC norm")
abline(h=-0.8, col="gray", lty=2)
dev.off()

# scatter plots

pdf("/mnt/projects/iamp/results/edaseq.pdf", width=30, height=20)
par(mfrow=c(5,7))
for (i in 1:ncol(dataNorm@assayData$counts)) {
  plot(dataNorm@assayData$counts[,i], dataNorm@assayData$normalizedCounts[,i], cex=0.2, log="xy")
}
dev.off()

sig <- !is.na(res[common, "padj"]) & res[common, "padj"] < 0.1
sig.norm <- !is.na(res.norm[common, "padj"]) & res.norm[common, "padj"] < 0.1

par(mfrow=c(1,1))
plot(res[common, "log2FoldChange"][sig | sig.norm], res.norm[common, "log2FoldChange"][sig | sig.norm], cex=0.3, col=rgb(0,0,0,0.6))
abline(0, 1, col="red")
abline(h=0, col="lightgray", lty=2)
abline(v=0, col="lightgray", lty=2)

plot(-log(res[common, "padj"], 10), -log(res.norm[common, "padj"], 10), cex=0.3, col=rgb(0,0,0,0.8))
abline(0, 1, col="red")
abline(h=10, col="lightgray", lty=2)
abline(v=10, col="lightgray", lty=2)

