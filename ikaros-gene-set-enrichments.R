library(clusterProfiler)

rm(list=ls())
result <- list()
min.expr <- 20

# read gene sets (.gmt file)

gmt <- read.delim("/mnt/projects/ikaros/data/ikaros_curated_genesets_gsea.gmt", stringsAsFactors = F, header = FALSE)
term2gene <- do.call("rbind", apply(gmt, 1, function(x) {
  data.frame(term=x[1], gene=unique(x[4:length(x)][x[4:length(x)]!=""]), row.names = NULL)
}))
term2gene <- rbind(term2gene, data.frame(term="ROBERTS_2014_PHLIKE_UP_OR_DN", gene=term2gene$gene[term2gene$term %in% c("ROBERTS_2014_PHLIKE_UP", "ROBERTS_2014_PHLIKE_DN")]))
term2gene <- rbind(term2gene, data.frame(term="IACOBUCCI_2012_IKAROS_DELETED_UP_OR_DN", gene=term2gene$gene[term2gene$term %in% c("IACOBUCCI_2012_IKAROS_DELETED_UP", "IACOBUCCI_2012_IKAROS_DELETED_DN")]))

# IKNp.vs.IKCp

m <- read.delim("/mnt/projects/ikaros/results/anduril/execute/deseqAnnotated_IKNp_vs_IKCp/table.csv", stringsAsFactors = F, row.names = 1)
expressed <- unique(m$Gene[m$meanA >= min.expr | m$meanB >= min.expr])
deg <- read.delim("/mnt/projects/ikaros/results/anduril/execute/degCalled_IKNp_vs_IKCp/table.csv", stringsAsFactors = F, row.names = 1)
IKNp.vs.IKCp <- unique(deg$Gene[!is.na(deg$Gene)])

result[["IKNp.vs.IKCp"]] <- enricher(  
  gene          = IKNp.vs.IKCp, 
  pvalueCutoff  = 1, 
  pAdjustMethod = "BH", 
  universe      = expressed, 
  minGSSize     = 5, 
  qvalueCutoff  = 1, 
  TERM2GENE     = rbind(term2gene, data.frame(term="expressed", gene=expressed)), 
  TERM2NAME     = gmt[,c(1,3)]
)@result

rownames(result[["IKNp.vs.IKCp"]]) <- NULL ; colnames(result[["IKNp.vs.IKCp"]])[1] <- "GeneSet"
result[["IKNp.vs.IKCp"]] <- cbind(data.frame(comparison="IKNp.vs.IKCp"), result[["IKNp.vs.IKCp"]])

# IKMp.vs.IKCp

m <- read.delim("/mnt/projects/ikaros/results/anduril/execute/deseqAnnotated_IKMp_vs_IKCp/table.csv", stringsAsFactors = F, row.names = 1)
expressed <- unique(m$Gene[m$meanA >= min.expr | m$meanB >= min.expr])
deg <- read.delim("/mnt/projects/ikaros/results/anduril/execute/degCalled_IKMp_vs_IKCp/table.csv", stringsAsFactors = F, row.names = 1)
IKMp.vs.IKCp <- unique(deg$Gene[!is.na(deg$Gene)])

result[["IKMp.vs.IKCp"]] <- enricher(
  gene          = IKMp.vs.IKCp, 
  pvalueCutoff  = 1, 
  pAdjustMethod = "BH", 
  universe      = expressed, 
  minGSSize     = 5, 
  qvalueCutoff  = 1, 
  TERM2GENE     = rbind(term2gene, data.frame(term="expressed", gene=expressed)), 
  TERM2NAME     = gmt[,c(1,3)]
)@result

rownames(result[["IKMp.vs.IKCp"]]) <- NULL ; colnames(result[["IKMp.vs.IKCp"]])[1] <- "GeneSet"
result[["IKMp.vs.IKCp"]] <- cbind(data.frame(comparison="IKMp.vs.IKCp"), result[["IKMp.vs.IKCp"]])

# IKD.vs.IKCp

m <- read.delim("/mnt/projects/ikaros/results/anduril/execute/deseqAnnotated_IKD_vs_IKCp/table.csv", stringsAsFactors = F, row.names = 1)
expressed <- unique(m$Gene[m$meanA >= min.expr | m$meanB >= min.expr])
deg <- read.delim("/mnt/projects/ikaros/results/anduril/execute/degCalled_IKD_vs_IKCp/table.csv", stringsAsFactors = F, row.names = 1)
IKD.vs.IKCp <- unique(deg$Gene[!is.na(deg$Gene)])

result[["IKD.vs.IKCp"]] <- enricher(
  gene          = IKD.vs.IKCp, 
  pvalueCutoff  = 1, 
  pAdjustMethod = "BH", 
  universe      = expressed, 
  minGSSize     = 5, 
  qvalueCutoff  = 1, 
  TERM2GENE     = rbind(term2gene, data.frame(term="expressed", gene=expressed)), 
  TERM2NAME     = gmt[,c(1,3)]
)@result
rownames(result[["IKD.vs.IKCp"]]) <- NULL ; colnames(result[["IKD.vs.IKCp"]])[1] <- "GeneSet"
result[["IKD.vs.IKCp"]] <- cbind(data.frame(comparison="IKD.vs.IKCp"), result[["IKD.vs.IKCp"]])

# IKDclean.vs.IKCp

m <- read.delim("/mnt/projects/ikaros/results/anduril/execute/deseqAnnotated_IKDclean_vs_IKCp/table.csv", stringsAsFactors = F, row.names = 1)
expressed <- unique(m$Gene[m$meanA >= min.expr | m$meanB >= min.expr])
deg <- read.delim("/mnt/projects/ikaros/results/anduril/execute/degCalled_IKDclean_vs_IKCp/table.csv", stringsAsFactors = F, row.names = 1)
IKDclean.vs.IKCp <- unique(deg$Gene[!is.na(deg$Gene)])

result[["IKDclean.vs.IKCp"]] <- enricher(
  gene          = IKDclean.vs.IKCp, 
  pvalueCutoff  = 1, 
  pAdjustMethod = "BH", 
  universe      = expressed, 
  minGSSize     = 5, 
  qvalueCutoff  = 1, 
  TERM2GENE     = rbind(term2gene, data.frame(term="expressed", gene=expressed)), 
  TERM2NAME     = gmt[,c(1,3)]
)@result
rownames(result[["IKDclean.vs.IKCp"]]) <- NULL ; colnames(result[["IKDclean.vs.IKCp"]])[1] <- "GeneSet"
result[["IKDclean.vs.IKCp"]] <- cbind(data.frame(comparison="IKDclean.vs.IKCp"), result[["IKDclean.vs.IKCp"]])

# write results

write.table(do.call("rbind", result), "/mnt/projects/ikaros/results/ikaros-gene-sets-enrichment.tsv", sep="\t", quote=F, col.names=T, row.names=F)