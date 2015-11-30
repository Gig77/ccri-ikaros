#source("http://bioconductor.org/biocLite.R")
#biocLite("hgu133plus2.db")
#biocLite("hgu133plus2probe")
#biocLite("hgu133plus2cdf")
library(annotate)
library(hgu133plus2.db)
library(hgu133plus2probe)
library(hgu133plus2cdf)

library("biomaRt")
mart <- useMart(biomart="ENSEMBL_MART_ENSEMBL", host="grch37.ensembl.org", path="/biomart/martservice" ,dataset="hsapiens_gene_ensembl") # GRCh37, v75
hgnc <- getBM(attributes=c("hgnc_symbol", "ensembl_gene_id"), mart=mart)
hgnc <- hgnc[hgnc$hgnc_symbol != "" & grepl("ENSG", hgnc$ensembl_gene_id),]

denBoer.tableS2 <- read.delim("/mnt/projects/ikaros/data/public/denBoer_2013_suppTableS2_Affymetrix U133 probesets classifier that was used for the identification of BCR-ABL1-like cases.txt", stringsAsFactors = F)
denBoer.tableS2.hgnc <- sort(unique(getSYMBOL(denBoer.tableS2$U133_.probesets, "hgu133plus2.db")))
denBoer.tableS2.ensembl <- merge(data.frame(hgnc=denBoer.tableS2.hgnc), hgnc, by.x="hgnc", by.y="hgnc_symbol")
write.table(denBoer.tableS2.ensembl[, c("ensembl_gene_id", "hgnc")], "/mnt/projects/ikaros/results/qlucore/denBoer_2013_suppTableS2_Affymetrix U133 probesets classifier that was used for the identification of BCR-ABL1-like cases.qlucore.txt", row.names = F, col.names=T, quote=F, sep="\t")

mullighan.tableS7a <- read.delim("/mnt/projects/ikaros/data/public/harvey_mullighan_2010_suppTable7a_ph-like-signature_probesets-used-in-clustering.txt", stringsAsFactors = F)
mullighan.tableS7a.hgnc <- sort(unique(getSYMBOL(mullighan.tableS7a$HC_P9906, "hgu133plus2.db")))
mullighan.tableS7a.ensembl <- merge(data.frame(hgnc=mullighan.tableS7a.hgnc), hgnc, by.x="hgnc", by.y="hgnc_symbol")
write.table(mullighan.tableS7a.ensembl[, c("ensembl_gene_id", "hgnc")], "/mnt/projects/ikaros/results/qlucore/harvey_mullighan_2010_suppTable7a_ph-like-signature_probesets-used-in-clustering.qlucore.txt", row.names = F, col.names=T, quote=F, sep="\t")

mullighan.tableS15 <- read.delim("/mnt/projects/ikaros/data/public/harvey_mullighan_2010_suppTable15_Top 100 Rank Order Genes Defining ROSE Cluster R8.txt", stringsAsFactors = F)
mullighan.tableS15.ensembl <- merge(data.frame(hgnc=unique(mullighan.tableS15$Symbol[!is.na(mullighan.tableS15$Symbol) & mullighan.tableS15$Symbol != ""])), hgnc, by.x="hgnc", by.y="hgnc_symbol")
write.table(mullighan.tableS15.ensembl[, c("ensembl_gene_id", "hgnc")], "/mnt/projects/ikaros/results/qlucore/harvey_mullighan_2010_suppTable15_Top 100 Rank Order Genes Defining ROSE Cluster R8.qlucore.txt", row.names = F, col.names=T, quote=F, sep="\t")

mullighan.tableS18 <- read.delim("/mnt/projects/ikaros/data/public/harvey_mullighan_2010_suppTable18_Genes Common to Rank Order and BCR-ABL1-like Signature.txt", stringsAsFactors = F)
mullighan.tableS18.up <- mullighan.tableS18$BCR_ABL_up[!is.na(mullighan.tableS18$BCR_ABL_up) & mullighan.tableS18$BCR_ABL_up != ""]
mullighan.tableS18.dn <- mullighan.tableS18$BCR_ABL_down[!is.na(mullighan.tableS18$BCR_ABL_down) & mullighan.tableS18$BCR_ABL_down != ""]
mullighan.tableS18.ensembl <- merge(data.frame(hgnc=unique(c(mullighan.tableS18.up, mullighan.tableS18.dn))), hgnc, by.x="hgnc", by.y="hgnc_symbol")
write.table(mullighan.tableS18.ensembl[, c("ensembl_gene_id", "hgnc")], "/mnt/projects/ikaros/results/qlucore/harvey_mullighan_2010_suppTable18_Genes Common to Rank Order and BCR-ABL1-like Signature.qlucore.txt", row.names = F, col.names=T, quote=F, sep="\t")

roberts.tableS10 <- read.delim("/mnt/projects/ikaros/data/public/roberts_2014_suppTable10.Ph-like.vs.other.txt", stringsAsFactors = F)
roberts.tableS10 <- roberts.tableS10[order(roberts.tableS10$Ph.like.vs.Other.qval),]
roberts.tableS10.up <- roberts.tableS10$Gene[roberts.tableS10$Ph.like.vs.Other.logFC >= 1][1:100]
roberts.tableS10.dn <- roberts.tableS10$Gene[roberts.tableS10$Ph.like.vs.Other.logFC <= -1][1:100]
roberts.tableS10.ensembl <- merge(data.frame(hgnc=unique(c(roberts.tableS10.up, roberts.tableS10.dn))), hgnc, by.x="hgnc", by.y="hgnc_symbol")
write.table(roberts.tableS10.ensembl[, c("ensembl_gene_id", "hgnc")], "/mnt/projects/ikaros/results/qlucore/roberts_2014_suppTable10.Ph-like.vs.other.qlucore.txt", row.names = F, col.names=T, quote=F, sep="\t")

# Gene list from Dagmar Schinnerl
# http://appft.uspto.gov/netacgi/nph-Parser?Sect1=PTO2&Sect2=HITOFF&u=%2Fnetahtml%2FPTO%2Fsearch-adv.html&r=1&p=1&f=G&l=50&d=PG01&S1=ph-like.TTL.&OS=ttl/ph-like&RS=TTL/ph-like
mullighan.patent <- c("IGJ", "SPATS2L", "MUC4", "CRLF2", "CA6", "NRXN3", "BMPR1B", "GPR110", "CHN2", "SEMA6A", "PON2", "SLC2A5", "S100Z", "TP53INP1", "IFITM1")
mullighan.patent.ensembl <- merge(data.frame(hgnc=mullighan.patent), hgnc, by.x="hgnc", by.y="hgnc_symbol")
write.table(mullighan.patent.ensembl[, c("ensembl_gene_id", "hgnc")], "/mnt/projects/ikaros/results/qlucore/mullighan_patent_phlike-15.qlucore.txt", row.names = F, col.names=T, quote=F, sep="\t")
