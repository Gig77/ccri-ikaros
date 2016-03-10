msigdb <- read.delim("/mnt/projects/generic/data/msigdb5.0/msigdb.v5.0.symbols.gmt", stringsAsFactors = F, header = F)
custom <- read.delim("/mnt/projects/ikaros/data/ikaros_curated_genesets_gsea.gmt", stringsAsFactors = F, header = F)
NvsC.dis <- read.delim("/mnt/projects/ikaros/results/anduril/execute/deseqAnnotated_IKNp_vs_IKCp/table.csv", stringsAsFactors = F)
NvsC.val <- read.delim("/mnt/projects/ikaros/results/anduril_validation/execute/deseqAnnotated_IKNv_vs_IKCv/table.csv", stringsAsFactors = F)

plotCorrelation <- function(m, title) {
  x <- m$fc.dis
  y <- m$fc.val
  fit <- lm(y~x)
  test <- cor.test(x,y,method="spearman")
  
  plot(x, y, xlab = "log2FC Discovery cohort", ylab = "log2FC Validation cohort", ylim=c(min(0, y), max(y)), xlim=c(min(0, x), max(x)), cex=0.4, pch=19, main=title)
  text(x, y-0.03, m$Gene, cex = 0.4)
  abline(fit, col="red")
  text((min(x)+max(x))/2, max(y), sprintf("R=%.2f/p=%.3g", test$estimate, test$p.value), adj=0.5)
}

#-------------------
# top-N DEGs
#-------------------

topN <- 50

NvsC.dis.topN <- NvsC.dis[!is.na(NvsC.dis$p) & NvsC.dis$p <= 1e-8,]
NvsC.dis.topN <- NvsC.dis.topN[order(NvsC.dis.topN$fc, decreasing=T),]
NvsC.dis.topN <- NvsC.dis.topN[c(1:topN,(nrow(NvsC.dis.topN)-topN+1):nrow(NvsC.dis.topN)),]


merged <- merge(NvsC.dis.topN[,c("ids", "Gene", "fc")], NvsC.val[,c("ids", "Gene", "fc")], by=c("ids", "Gene"), suffixes=c(".dis", ".val"))

pdf("/mnt/projects/ikaros/results/correlation-discovery-vs-validation.pdf")
plotCorrelation(merged, "Top-50 up- and downregulated genes\n from discovery cohort (IKNp vs. IKCp)")
dev.off()

#-------------------
# PID FAK Pathway
#-------------------

pid.fak.genes <- grep(".", as.character(msigdb[msigdb$V1=="PID_FAK_PATHWAY",c(-1, -2)]), value=TRUE)
pid.fak.expr <- NvsC.dis[NvsC.dis$Gene %in% pid.fak.genes & !is.na(NvsC.dis$q) & (NvsC.dis$meanA >= 100 | NvsC.dis$meanB >= 100), c("ids", "Gene", "fc")]
pid.fak.expr <- merge(pid.fak.expr, NvsC.val[,c("ids", "Gene", "fc")], by=c("ids", "Gene"), suffixes=c(".dis", ".val"))
plotCorrelation(pid.fak.expr, "PID_FAK_PATHWAY")

#-------------------
# PID FAK Pathway
#-------------------

pid.vegfr.genes <- grep(".", as.character(msigdb[msigdb$V1=="PID_VEGFR1_2_PATHWAY",c(-1, -2)]), value=TRUE)
pid.vegfr.expr <- NvsC.dis[NvsC.dis$Gene %in% pid.vegfr.genes & !is.na(NvsC.dis$q) & (NvsC.dis$meanA >= 100 | NvsC.dis$meanB >= 100), c("ids", "Gene", "fc")]
pid.vegfr.expr <- merge(pid.vegfr.expr, NvsC.val[,c("ids", "Gene", "fc")], by=c("ids", "Gene"), suffixes=c(".dis", ".val"))
plotCorrelation(pid.vegfr.expr, "PID_VEGFR1_2_PATHWAY")

#-------------------
# FERREIROS_2013_IK6_OVEREXPR_48H_DN
#-------------------

ferreiros2013.IK6.48h.dn.genes <- grep(".", as.character(custom[custom$V1=="FERREIROS_2013_IK6_OVEREXPR_48H_DN",c(-1, -2, -3)]), value=TRUE)
ferreiros2013.IK6.48h.dn.expr <- NvsC.dis[NvsC.dis$Gene %in% ferreiros2013.IK6.48h.dn.genes & !is.na(NvsC.dis$q) & (NvsC.dis$meanA >= 100 | NvsC.dis$meanB >= 100), c("ids", "Gene", "fc")]
ferreiros2013.IK6.48h.dn.expr <- merge(ferreiros2013.IK6.48h.dn.expr, NvsC.val[,c("ids", "Gene", "fc")], by=c("ids", "Gene"), suffixes=c(".dis", ".val"))
plotCorrelation(ferreiros2013.IK6.48h.dn.expr, "FERREIROS_2013_IK6_OVEREXPR_48H_DN")

#-------------------
# FERREIROS_2013_IK6_OVEREXPR_48H_UP
#-------------------

ferreiros2013.IK6.48h.up.genes <- grep(".", as.character(custom[custom$V1=="FERREIROS_2013_IK6_OVEREXPR_48H_UP",c(-1, -2, -3)]), value=TRUE)
ferreiros2013.IK6.48h.up.expr <- NvsC.dis[NvsC.dis$Gene %in% ferreiros2013.IK6.48h.up.genes & !is.na(NvsC.dis$q) & (NvsC.dis$meanA >= 100 | NvsC.dis$meanB >= 100), c("ids", "Gene", "fc")]
ferreiros2013.IK6.48h.up.expr <- merge(ferreiros2013.IK6.48h.up.expr, NvsC.val[,c("ids", "Gene", "fc")], by=c("ids", "Gene"), suffixes=c(".dis", ".val"))
plotCorrelation(ferreiros2013.IK6.48h.up.expr, "FERREIROS_2013_IK6_OVEREXPR_48H_UP")

#-------------------
# SCHWICKERT_2014_IKAROS_REPRESSED
#-------------------

schwickert2014.ikaros.repressed.genes <- grep(".", as.character(custom[custom$V1=="SCHWICKERT_2014_IKAROS_REPRESSED",c(-1, -2, -3)]), value=TRUE)
schwickert2014.ikaros.repressed.expr <- NvsC.dis[NvsC.dis$Gene %in% schwickert2014.ikaros.repressed.genes & !is.na(NvsC.dis$q) & (NvsC.dis$meanA >= 100 | NvsC.dis$meanB >= 100), c("ids", "Gene", "fc")]
schwickert2014.ikaros.repressed.expr <- merge(schwickert2014.ikaros.repressed.expr, NvsC.val[,c("ids", "Gene", "fc")], by=c("ids", "Gene"), suffixes=c(".dis", ".val"))
plotCorrelation(schwickert2014.ikaros.repressed.expr, "SCHWICKERT_2014_IKAROS_REPRESSED")

#-------------------
# IACOBUCCI_2012_IKAROS_DELETED_UP
#-------------------

iacobucci2012.ikaros.deleted.up.genes <- grep(".", as.character(custom[custom$V1=="IACOBUCCI_2012_IKAROS_DELETED_UP",c(-1, -2, -3)]), value=TRUE)
iacobucci2012.ikaros.deleted.up.expr <- NvsC.dis[NvsC.dis$Gene %in% iacobucci2012.ikaros.deleted.up.genes & !is.na(NvsC.dis$q) & (NvsC.dis$meanA >= 100 | NvsC.dis$meanB >= 100), c("ids", "Gene", "fc")]
iacobucci2012.ikaros.deleted.up.expr <- merge(iacobucci2012.ikaros.deleted.up.expr, NvsC.val[,c("ids", "Gene", "fc")], by=c("ids", "Gene"), suffixes=c(".dis", ".val"))
plotCorrelation(iacobucci2012.ikaros.deleted.up.expr, "IACOBUCCI_2012_IKAROS_DELETED_UP")

#-------------------
# IACOBUCCI_2012_IKAROS_DELETED_DN
#-------------------

iacobucci2012.ikaros.deleted.dn.genes <- grep(".", as.character(custom[custom$V1=="IACOBUCCI_2012_IKAROS_DELETED_DN",c(-1, -2, -3)]), value=TRUE)
iacobucci2012.ikaros.deleted.dn.expr <- NvsC.dis[NvsC.dis$Gene %in% iacobucci2012.ikaros.deleted.dn.genes & !is.na(NvsC.dis$q) & (NvsC.dis$meanA >= 100 | NvsC.dis$meanB >= 100), c("ids", "Gene", "fc")]
iacobucci2012.ikaros.deleted.dn.expr <- merge(iacobucci2012.ikaros.deleted.dn.expr, NvsC.val[,c("ids", "Gene", "fc")], by=c("ids", "Gene"), suffixes=c(".dis", ".val"))
plotCorrelation(iacobucci2012.ikaros.deleted.dn.expr, "IACOBUCCI_2012_IKAROS_DELETED_DN")
