# get gene gc content and length

library(biomaRt)
mart = useMart(biomart="ENSEMBL_MART_ENSEMBL", host="grch37.ensembl.org", path="/biomart/martservice", dataset="hsapiens_gene_ensembl") # GRCh37, v75
hgnc <- getBM(attributes=c("hgnc_symbol", "ensembl_gene_id"), mart=mart)
gclength <- getBM(attributes=c("ensembl_gene_id", "percentage_gc_content", "transcript_length"), mart=mart)
gclength <- gclength[order(gclength$ensembl_gene_id, -gclength$transcript_length),]
gclength <- gclength[!duplicated(gclength$ensembl_gene_id),]
gclength <- data.frame(id=gclength$ensembl_gene_id, gc=gclength$percentage_gc_content, length=gclength$transcript_length)

# read and filter differentially expressed genes

IKNp.vs.IKCp <- read.delim("/mnt/projects/ikaros/results/anduril/execute/deseqAnnotated_IKNp_vs_IKCp/table.csv", stringsAsFactors = FALSE)
IKNp.vs.IKCp <- merge(IKNp.vs.IKCp, gclength, by.x="ids", by.y="id")
IKNp.vs.IKCp <- IKNp.vs.IKCp[!is.na(IKNp.vs.IKCp$q) & IKNp.vs.IKCp$p <= 1e-2,]

# scatter plot

length <- log10(IKNp.vs.IKCp$length)
signif <- -log10(IKNp.vs.IKCp$p)
logfc <- IKNp.vs.IKCp$fc

plot(length, signif, ylim=c(min(signif),15), xlab="log10 transcript length in bp", ylab="-log10(p)", main="IKAROS IKNp vs. IKCp")
fit <- lm(signif~length)
abline(fit, col="red")
test <- cor.test(signif,length,method="spearman")

plot(length, logfc, xlab="log10 transcript length in bp", ylab="log2fc", main="IKAROS IKNp vs. IKCp", cex=0.3, col=rgb(0,0,0,0.5))
fit <- lm(logfc~length)
abline(fit, col="red")
test <- cor.test(signif,logfc,method="spearman")
text(2, 0, sprintf("R=%.2f/p=%.3g", test$estimate, test$p.value), adj=0)
