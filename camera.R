library(edgeR)
library(biomaRt)

mart = useMart(biomart="ENSEMBL_MART_ENSEMBL", host="grch37.ensembl.org", path="/biomart/martservice", dataset="hsapiens_gene_ensembl") # GRCh37, v75
hgnc <- getBM(attributes=c("hgnc_symbol", "ensembl_gene_id"), mart=mart)

# get sample annotation and count matrix

sampleAnn <- read.delim("/mnt/projects/ikaros/data/samples.csv", row.names="Alias", stringsAsFactors = FALSE)
x <- read.delim("/mnt/projects/ikaros/results/anduril/execute/htseqExprMatrix/countArray/all.csv", row.names=1)

# parse gene sets from GMT file, convert to index lists

msigdb <- list()
indexsets <- list()
conn <- file("/mnt/projects/generic/data/msigdb5.0/msigdb.v5.0.symbols.gmt", open="r")
msigdb.raw <- readLines(conn)
close(conn)
for (i in 1:length(msigdb.raw)) {
  cols <- unlist(strsplit(msigdb.raw[i], "\t"))
  gs.hgnc <- cols[c(-1,-2)]
  gs.ensembl <- na.omit(hgnc$ensembl_gene_id[match(gs.hgnc, hgnc$hgnc_symbol)])
  msigdb[[cols[1]]] <- list()
  msigdb[[cols[1]]][["URL"]] <- cols[2]
  msigdb[[cols[1]]][["genes"]] <- paste(cols[3:length(cols)], collapse=",")
  indexsets[[cols[1]]] <- which(rownames(x) %in% gs.ensembl)
}
msigdb <- do.call(rbind.data.frame, msigdb)

# define sample groups

samples.ikn <- rownames(sampleAnn[sampleAnn$IKdn == "yes" & sampleAnn$Xeno != "yes" & rownames(sampleAnn) != "IK6_2",])
samples.ikd <- rownames(sampleAnn[sampleAnn$Group == "IKD" & sampleAnn$Xeno != "yes" & !rownames(sampleAnn) %in% samples.ikn,])
samples.ikc <- rownames(sampleAnn[sampleAnn$Group == "IKC" & sampleAnn$Xeno != "yes",])

# IKN vs. IKC

group <- as.factor(sampleAnn[c(samples.ikn, samples.ikc), "IKdn"])
y <- DGEList(counts=x[,c(samples.ikn, samples.ikc)], group=group)
y <- calcNormFactors(y)
y <- voom(y)
design <- model.matrix(~group)

result.IKNvsIKC <- camera(y, indexsets, design)
hist(result.IKNvsIKC$PValue)

# (try mroast as well...)

y <- DGEList(counts=x[,c(samples.ikn, samples.ikc)], group=group)
y <- calcNormFactors(y)
y <- estimateDisp(y, design, robust=TRUE)
result.IKNvsIKC.mroast <- mroast(y, indexsets, design, nrot=10000)
hist(result.IKNvsIKC.mroast$PValue)


# IKD vs. IKC

group <- as.factor(sampleAnn[c(samples.ikd, samples.ikc), "Group"])
y <- DGEList(counts=x[,c(samples.ikd, samples.ikc)], group=group)
y <- calcNormFactors(y)
y <- voom(y)
design <- model.matrix(~group)

result.IKDvsIKC <- camera(y, indexsets, design)
hist(result.IKDvsIKC$PValue)

# IKM vs. IKC

group <- as.factor(gsub("(IK6|IKD)", "IKM", sampleAnn[c(samples.ikn, samples.ikd, samples.ikc), "Group"]))
y <- DGEList(counts=x[,c(samples.ikn, samples.ikd, samples.ikc)], group=group)
y <- calcNormFactors(y)
y <- voom(y)
design <- model.matrix(~group)

result.IKMvsIKC <- camera(y, indexsets, design)
hist(result.IKMvsIKC$PValue)

# control: enrichment in randomized samples

#group <- as.factor(sample(sampleAnn[c(samples.ikn, samples.ikc), "IKdn"]))
#y <- DGEList(counts=x[,c(samples.ikn, samples.ikc)], group=group)
#y <- calcNormFactors(y)
#y <- voom(y)
#design <- model.matrix(~group)

#result.random <- camera(y, indexsets, design)
#hist(result.random$PValue)

# merge and annotate results, write out

result <- merge(result.IKNvsIKC, result.IKDvsIKC, by=c("row.names", "NGenes"), suffixes=c(".IKNvsIKC", ".IKDvsIKC"), all=T)
result <- merge(result, result.IKMvsIKC, by.x=c("Row.names", "NGenes"), by.y=c("row.names", "NGenes"), all=T)
names(result)[names(result)=="Row.names"] <- "Geneset"
names(result)[11:14] <- paste0(names(result)[11:14], ".IKMvsIKC")
result <- merge(result, msigdb, by.x="Geneset", by.y="row.names")
result <- result[!grepl("(V\\$|MODULE_|_UNKNOWN)", result$Geneset),]

write.table(result, "/mnt/projects/ikaros/results/enriched-gene-sets.camera.tsv", row.names=F, col.names=T, sep="\t", quote=F, na="")

