library(biomaRt)
human = useMart(biomart="ENSEMBL_MART_ENSEMBL", host="grch37.ensembl.org", path="/biomart/martservice", dataset="hsapiens_gene_ensembl") # GRCh37, v75
mouse = useMart(biomart="ENSEMBL_MART_ENSEMBL", host="feb2014.archive.ensembl.org", dataset="mmusculus_gene_ensembl") # GRCm38, v75

fc.cutoff <- 20

schwickert.activated <- read.delim("/mnt/projects/ikaros/data/public/schwickert_2014_suppl_table_1_in_vitro_ikaros_activated.txt", stringsAsFactors = F)
schwickert.activated$fc <- -schwickert.activated$fc
schwickert.repressed <- read.delim("/mnt/projects/ikaros/data/public/schwickert_2014_suppl_table_1_in_vitro_ikaros_repressed.txt", stringsAsFactors = F)
schwickert <- rbind(schwickert.activated, schwickert.repressed)
schwickert <- schwickert[order(abs(schwickert$fc), decreasing=T),]
schwickert <- schwickert[!duplicated(schwickert$gene),]
schwickert$fc[schwickert$fc > fc.cutoff] <- fc.cutoff
schwickert$fc[schwickert$fc < -fc.cutoff] <- -fc.cutoff

humOrt <- getLDS(attributes = c("mgi_symbol"),
                 filters = "mgi_symbol", values=c(schwickert.activated$gene, schwickert.repressed$gene), mart = mouse,
                 attributesL = c("hgnc_symbol", "entrezgene"), martL = human)
colnames(humOrt) <- c("mgi_symbol", "hgnc", "entrezgene")
humOrt <- humOrt[humOrt$hgnc != "",]

schwickert <- merge(schwickert, humOrt[,c("mgi_symbol", "hgnc")], by.x="gene", by.y="mgi_symbol")
schwickert <- schwickert[order(abs(schwickert$fc), decreasing=T),]
schwickert <- schwickert[!duplicated(schwickert$hgnc),]

pdf("/mnt/projects/ikaros/results/schwickert2014.pdf", height=12, width=6)
par(mfrow=c(2,1))

# IKDvsIKCp
IKDvsIKCp <- read.delim("/mnt/projects/ikaros/results/anduril/execute/deseqAnnotated_IKD_vs_IKCp/table.csv", stringsAsFactors = F)
m <- merge(schwickert, IKDvsIKCp[,c("Gene", "fc", "q")], by.x="hgnc", by.y="Gene", suffixes=c(".schwickert", ".conny"))
fit <- lm(fc.conny~fc.schwickert, m)
p <- anova(fit)$'Pr(>F)'[1]
R <- summary(fit)$r.squared
plot(m$fc.schwickert, m$fc.conny, xlim=c(-fc.cutoff-1, fc.cutoff+1), main=sprintf("Conny IKDvsIKCp vs. Schwickert2014 IKinactvsCtrl\np=%.2g, R=%.2g", p, R), xlab="FC Schwickert", ylab="FC Conny", pch=20, cex=0.3, col=ifelse(is.na(m$q.conny) | m$q.conny > 0.1, "black", "red"))
abline(fit, col="red")
labels <- m[order(rowMeans(scale(m[,c("fc.schwickert", "fc.conny")])), decreasing=F), c("fc.schwickert", "fc.conny", "q.conny", "hgnc")]
with(labels[1:20,], text(fc.schwickert+0.05, fc.conny+0.05, hgnc, cex=0.4, adj=0, col=ifelse(is.na(q.conny) | q.conny > 0.1, "black", "red")))
with(labels[(nrow(labels)-20):nrow(labels),], text(fc.schwickert+0.05, fc.conny+0.05, hgnc, cex=0.3, adj=0, col=ifelse(is.na(q.conny) | q.conny > 0.1, "black", "red")))


# IKNpvsIKCp
IKNpvsIKCp <- read.delim("/mnt/projects/ikaros/results/anduril/execute/deseqAnnotated_IKNp_vs_IKCp/table.csv", stringsAsFactors = F)
m <- merge(schwickert, IKNpvsIKCp[,c("Gene", "fc", "q")], by.x="hgnc", by.y="Gene", suffixes=c(".schwickert", ".conny"))
fit <- lm(fc.conny~fc.schwickert, m)
p <- anova(fit)$'Pr(>F)'[1]
R <- summary(fit)$r.squared
plot(m$fc.schwickert, m$fc.conny, xlim=c(-fc.cutoff-1, fc.cutoff+1), main=sprintf("Conny IKNpvsIKCp vs. Schwickert2014 IKinactvsCtrl\np=%.2g, R=%.2g", p, R), xlab="FC Schwickert", ylab="FC Conny", pch=20, cex=0.3, col=ifelse(is.na(m$q.conny) | m$q.conny > 0.1, "black", "red"))
abline(fit, col="red")
labels <- m[order(rowMeans(scale(m[,c("fc.schwickert", "fc.conny")])), decreasing=F), c("fc.schwickert", "fc.conny", "q.conny", "hgnc")]
with(labels[1:20,], text(fc.schwickert+0.05, fc.conny+0.05, hgnc, cex=0.4, adj=0, col=ifelse(is.na(q.conny) | q.conny > 0.1, "black", "red")))
with(labels[(nrow(labels)-20):nrow(labels),], text(fc.schwickert+0.05, fc.conny+0.05, hgnc, cex=0.3, adj=0, col=ifelse(is.na(q.conny) | q.conny > 0.1, "black", "red")))

dev.off()