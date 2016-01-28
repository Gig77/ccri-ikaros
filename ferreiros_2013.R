library(biomaRt)
human = useMart(biomart="ENSEMBL_MART_ENSEMBL", host="grch37.ensembl.org", path="/biomart/martservice", dataset="hsapiens_gene_ensembl") # GRCh37, v75
mouse = useMart(biomart="ENSEMBL_MART_ENSEMBL", host="feb2014.archive.ensembl.org", dataset="mmusculus_gene_ensembl") # GRCm38, v75

pdf("/mnt/projects/ikaros/results/ferreiros2013.pdf", height=12, width=6)
par(mfrow=c(2,1))

# -------------------------
# 2 hours
# -------------------------

ferreiros.2h <- read.delim("/mnt/projects/ikaros/data/public/Ferreiros-Vidal_2013_SuppData1_2hrs.txt", stringsAsFactors = F)
humOrt <- getLDS(attributes = c("mgi_symbol"), filters = "mgi_symbol", values=ferreiros.2h$MGI, mart = mouse, attributesL = c("hgnc_symbol", "entrezgene"), martL = human)
colnames(humOrt) <- c("mgi_symbol", "hgnc", "entrezgene")
humOrt <- humOrt[humOrt$hgnc != "",]
ferreiros.2h <- merge(ferreiros.2h, humOrt[,c("mgi_symbol", "hgnc")], by.x="MGI", by.y="mgi_symbol")
ferreiros.2h <- ferreiros.2h[order(abs(ferreiros.2h$fc), decreasing = T),]
ferreiros.2h <- ferreiros.2h[!duplicated(ferreiros.2h$hgnc),]

# IKDvsIKCp 2h

outliers <- c("IKZF1", "SDC4", "SEMA6B")
IKDvsIKCp <- read.delim("/mnt/projects/ikaros/results/anduril/execute/deseqAnnotated_IKD_vs_IKCp/table.csv", stringsAsFactors = F)
m <- merge(ferreiros.2h, IKDvsIKCp[,c("Gene", "fc", "q")], by.x="hgnc", by.y="Gene", suffixes=c(".ferreiros.2h", ".conny"))
m <- m[!m$hgnc %in% outliers,] # remove outliers
fit <- lm(fc.conny~fc.ferreiros.2h, m)
p <- anova(fit)$'Pr(>F)'[1]
R <- summary(fit)$r.squared
plot(m$fc.ferreiros.2h, m$fc.conny, xlim=c(min(m$fc.ferreiros.2h)-0.2, max(m$fc.ferreiros.2h)+0.2), main=sprintf("Conny IKDvsIKCp vs. Ferreiros2013 oeIK6vsCtrl2h\np=%.2g, R=%.2g", p, R), xlab="FC Ferreiros-Vidal", ylab="FC Conny", pch=20, cex=0.3, col=ifelse(is.na(m$q.conny) | m$q.conny > 0.1, "black", "red"))
abline(h=0, lty=3)
abline(v=0, lty=3)
abline(fit, col="red")
labels <- m[order(rowMeans(scale(m[,c("fc.ferreiros.2h", "fc.conny")])), decreasing=F), c("fc.ferreiros.2h", "fc.conny", "q.conny", "hgnc", "peaked")]
with(labels[1:20,], text(fc.ferreiros.2h+0.02, fc.conny+0.02, hgnc, cex=0.3, adj=0, col=ifelse(is.na(q.conny) | q.conny > 0.1, "black", "red"), font=1+peaked))
with(labels[(nrow(labels)-20):nrow(labels),], text(fc.ferreiros.2h+0.02, fc.conny+0.02, hgnc, cex=0.3, adj=0, col=ifelse(is.na(q.conny) | q.conny > 0.1, "black", "red"), font=1+peaked))
text(0.1, -1.4, paste0("Removed outliers: ", paste0(outliers, collapse=", ")), adj=0, cex=0.5)

# IKNpvsIKCp 2h

outliers <- c("IKZF1")
IKNpvsIKCp <- read.delim("/mnt/projects/ikaros/results/anduril/execute/deseqAnnotated_IKNp_vs_IKCp/table.csv", stringsAsFactors = F)
m <- merge(ferreiros.2h, IKNpvsIKCp[,c("Gene", "fc", "q")], by.x="hgnc", by.y="Gene", suffixes=c(".ferreiros.2h", ".conny"))
m <- m[!m$hgnc %in% outliers,] # remove outliers
fit <- lm(fc.conny~fc.ferreiros.2h, m)
p <- anova(fit)$'Pr(>F)'[1]
R <- summary(fit)$r.squared
plot(m$fc.ferreiros.2h, m$fc.conny, xlim=c(min(m$fc.ferreiros.2h)-0.2, max(m$fc.ferreiros.2h)+0.2), main=sprintf("Conny IKNpvsIKCp vs. Ferreiros2013 oeIK6vsCtrl2h\np=%.2g, R=%.2g", p, R), xlab="FC Ferreiros-Vidal", ylab="FC Conny", pch=20, cex=0.3, col=ifelse(is.na(m$q.conny) | m$q.conny > 0.1, "black", "red"))
abline(h=0, lty=3)
abline(v=0, lty=3)
abline(fit, col="red")
labels <- m[order(rowMeans(scale(m[,c("fc.ferreiros.2h", "fc.conny")])), decreasing=F), c("fc.ferreiros.2h", "fc.conny", "q.conny", "hgnc", "peaked")]
with(labels[1:20,], text(fc.ferreiros.2h+0.02, fc.conny+0.02, hgnc, cex=0.3, adj=0, col=ifelse(is.na(q.conny) | q.conny > 0.1, "black", "red"), font=1+peaked))
with(labels[(nrow(labels)-20):nrow(labels),], text(fc.ferreiros.2h+0.02, fc.conny+0.02, hgnc, cex=0.3, adj=0, col=ifelse(is.na(q.conny) | q.conny > 0.1, "black", "red"), font=1+peaked))
text(0.1, -2, paste0("Removed outliers: ", paste0(outliers, collapse=", ")), adj=0, cex=0.5)

# -------------------------
# 6 hours
# -------------------------

ferreiros.6h <- read.delim("/mnt/projects/ikaros/data/public/Ferreiros-Vidal_2013_SuppData1_6hrs.txt", stringsAsFactors = F)
humOrt <- getLDS(attributes = c("mgi_symbol"), filters = "mgi_symbol", values=ferreiros.6h$MGI, mart = mouse, attributesL = c("hgnc_symbol", "entrezgene"), martL = human)
colnames(humOrt) <- c("mgi_symbol", "hgnc", "entrezgene")
humOrt <- humOrt[humOrt$hgnc != "",]
ferreiros.6h <- merge(ferreiros.6h, humOrt[,c("mgi_symbol", "hgnc")], by.x="MGI", by.y="mgi_symbol")
ferreiros.6h <- ferreiros.6h[order(abs(ferreiros.6h$fc), decreasing = T),]
ferreiros.6h <- ferreiros.6h[!duplicated(ferreiros.6h$hgnc),]

# IKDvsIKCp 6h

outliers <- c("IKZF1")
IKDvsIKCp <- read.delim("/mnt/projects/ikaros/results/anduril/execute/deseqAnnotated_IKD_vs_IKCp/table.csv", stringsAsFactors = F)
m <- merge(ferreiros.6h, IKDvsIKCp[,c("Gene", "fc", "q")], by.x="hgnc", by.y="Gene", suffixes=c(".ferreiros.6h", ".conny"))
m <- m[!m$hgnc %in% outliers,] # remove outliers
fit <- lm(fc.conny~fc.ferreiros.6h, m)
p <- anova(fit)$'Pr(>F)'[1]
R <- summary(fit)$r.squared
plot(m$fc.ferreiros.6h, m$fc.conny, xlim=c(min(m$fc.ferreiros.6h)-0.2, max(m$fc.ferreiros.6h)+0.2), main=sprintf("Conny IKDvsIKCp vs. Ferreiros2013 oeIK6vsCtrl6h\np=%.2g, R=%.2g", p, R), xlab="FC Ferreiros-Vidal", ylab="FC Conny", pch=20, cex=0.3, col=ifelse(is.na(m$q.conny) | m$q.conny > 0.1, "black", "red"))
abline(h=0, lty=3)
abline(v=0, lty=3)
abline(fit, col="red")
labels <- m[order(rowMeans(scale(m[,c("fc.ferreiros.6h", "fc.conny")])), decreasing=F), c("fc.ferreiros.6h", "fc.conny", "q.conny", "hgnc", "peaked")]
with(labels[1:20,], text(fc.ferreiros.6h+0.03, fc.conny+0.05, hgnc, cex=0.3, adj=0, col=ifelse(is.na(q.conny) | q.conny > 0.1, "black", "red"), font=1+peaked))
with(labels[(nrow(labels)-20):nrow(labels),], text(fc.ferreiros.6h+0.03, fc.conny+0.05, hgnc, cex=0.3, adj=0, col=ifelse(is.na(q.conny) | q.conny > 0.1, "black", "red"), font=1+peaked))
text(0.1, -2.5, paste0("Removed outliers: ", paste0(outliers, collapse=", ")), adj=0, cex=0.5)

# IKNpvsIKCp 6h

outliers <- c("IKZF1")
IKNpvsIKCp <- read.delim("/mnt/projects/ikaros/results/anduril/execute/deseqAnnotated_IKNp_vs_IKCp/table.csv", stringsAsFactors = F)
m <- merge(ferreiros.6h, IKNpvsIKCp[,c("Gene", "fc", "q")], by.x="hgnc", by.y="Gene", suffixes=c(".ferreiros.6h", ".conny"))
m <- m[!m$hgnc %in% outliers,] # remove outliers
fit <- lm(fc.conny~fc.ferreiros.6h, m)
p <- anova(fit)$'Pr(>F)'[1]
R <- summary(fit)$r.squared
plot(m$fc.ferreiros.6h, m$fc.conny, xlim=c(min(m$fc.ferreiros.6h)-0.2, max(m$fc.ferreiros.6h)+0.2), main=sprintf("Conny IKNpvsIKCp vs. Ferreiros2013 oeIK6vsCtrl6h\np=%.2g, R=%.2g", p, R), xlab="FC Ferreiros-Vidal", ylab="FC Conny", pch=20, cex=0.3, col=ifelse(is.na(m$q.conny) | m$q.conny > 0.1, "black", "red"))
abline(h=0, lty=3)
abline(v=0, lty=3)
abline(fit, col="red")
labels <- m[order(rowMeans(scale(m[,c("fc.ferreiros.6h", "fc.conny")])), decreasing=F), c("fc.ferreiros.6h", "fc.conny", "q.conny", "hgnc", "peaked")]
with(labels[1:20,], text(fc.ferreiros.6h+0.02, fc.conny+0.05, hgnc, cex=0.3, adj=0, col=ifelse(is.na(q.conny) | q.conny > 0.1, "black", "red"), font=1+peaked))
with(labels[(nrow(labels)-20):nrow(labels),], text(fc.ferreiros.6h+0.02, fc.conny+0.05, hgnc, cex=0.3, adj=0, col=ifelse(is.na(q.conny) | q.conny > 0.1, "black", "red"), font=1+peaked))
text(0.1, -3, paste0("Removed outliers: ", paste0(outliers, collapse=", ")), adj=0, cex=0.5)

# -------------------------
# 48 hours
# -------------------------

ferreiros.48h <- read.delim("/mnt/projects/ikaros/data/public/Ferreiros-Vidal_2013_SuppData1_48hrs.txt", stringsAsFactors = F)
humOrt <- getLDS(attributes = c("mgi_symbol"), filters = "mgi_symbol", values=ferreiros.48h$MGI, mart = mouse, attributesL = c("hgnc_symbol", "entrezgene"), martL = human)
colnames(humOrt) <- c("mgi_symbol", "hgnc", "entrezgene")
humOrt <- humOrt[humOrt$hgnc != "",]
ferreiros.48h <- merge(ferreiros.48h, humOrt[,c("mgi_symbol", "hgnc")], by.x="MGI", by.y="mgi_symbol")
ferreiros.48h <- ferreiros.48h[order(abs(ferreiros.48h$fc), decreasing = T),]
ferreiros.48h <- ferreiros.48h[!duplicated(ferreiros.48h$hgnc),]
#ferreiros.48h <- ferreiros.48h[ferreiros.48h$FDR<0.01,]

# IKDvsIKCp 48h

outliers <- c("IKZF1")
IKDvsIKCp <- read.delim("/mnt/projects/ikaros/results/anduril/execute/deseqAnnotated_IKD_vs_IKCp/table.csv", stringsAsFactors = F)
m <- merge(ferreiros.48h, IKDvsIKCp[,c("Gene", "fc", "q")], by.x="hgnc", by.y="Gene", suffixes=c(".ferreiros.48h", ".conny"))
m <- m[!m$hgnc %in% outliers,] # remove outliers
fit <- lm(fc.conny~fc.ferreiros.48h, m)
p <- anova(fit)$'Pr(>F)'[1]
R <- summary(fit)$r.squared
plot(m$fc.ferreiros.48h, m$fc.conny, xlim=c(min(m$fc.ferreiros.48h)-0.2, max(m$fc.ferreiros.48h)+0.2), main=sprintf("Conny IKDvsIKCp vs. Ferreiros2013 oeIK6vsCtrl48h\np=%.2g, R=%.2g", p, R), xlab="FC Ferreiros-Vidal", ylab="FC Conny", pch=20, cex=0.3, col=ifelse(is.na(m$q.conny) | m$q.conny > 0.1, "black", "red"))
abline(h=0, lty=3)
abline(v=0, lty=3)
abline(fit, col="red")
labels <- m[order(rowMeans(scale(m[,c("fc.ferreiros.48h", "fc.conny")])), decreasing=F), c("fc.ferreiros.48h", "fc.conny", "q.conny", "hgnc", "peaked")]
with(labels[1:20,], text(fc.ferreiros.48h+0.1, fc.conny+0.05, hgnc, cex=0.3, adj=0, col=ifelse(is.na(q.conny) | q.conny > 0.1, "black", "red"), font=1+peaked))
with(labels[(nrow(labels)-20):nrow(labels),], text(fc.ferreiros.48h+0.1, fc.conny+0.05, hgnc, cex=0.3, adj=0, col=ifelse(is.na(q.conny) | q.conny > 0.1, "black", "red"), font=1+peaked))
text(1, -2.5, paste0("Removed outliers: ", paste0(outliers, collapse=", ")), adj=0, cex=0.5)

# IKNpvsIKCp 48h

outliers <- c("IKZF1")
IKNpvsIKCp <- read.delim("/mnt/projects/ikaros/results/anduril/execute/deseqAnnotated_IKNp_vs_IKCp/table.csv", stringsAsFactors = F)
m <- merge(ferreiros.48h, IKNpvsIKCp[,c("Gene", "fc", "q")], by.x="hgnc", by.y="Gene", suffixes=c(".ferreiros.48h", ".conny"))
m <- m[!m$hgnc %in% outliers,] # remove outliers
fit <- lm(fc.conny~fc.ferreiros.48h, m)
p <- anova(fit)$'Pr(>F)'[1]
R <- summary(fit)$r.squared
plot(m$fc.ferreiros.48h, m$fc.conny, xlim=c(min(m$fc.ferreiros.48h)-0.2, max(m$fc.ferreiros.48h)+0.2), main=sprintf("Conny IKNpvsIKCp vs. Ferreiros2013 oeIK6vsCtrl48h\np=%.2g, R=%.2g", p, R), xlab="FC Ferreiros-Vidal", ylab="FC Conny", pch=20, cex=0.3, col=ifelse(is.na(m$q.conny) | m$q.conny > 0.1, "black", "red"))
abline(h=0, lty=3)
abline(v=0, lty=3)
abline(fit, col="red")
labels <- m[order(rowMeans(scale(m[,c("fc.ferreiros.48h", "fc.conny")])), decreasing=F), c("fc.ferreiros.48h", "fc.conny", "q.conny", "hgnc", "peaked")]
with(labels[1:20,], text(fc.ferreiros.48h+0.1, fc.conny+0.05, hgnc, cex=0.3, adj=0, col=ifelse(is.na(q.conny) | q.conny > 0.1, "black", "red"), font=1+peaked))
with(labels[(nrow(labels)-20):nrow(labels),], text(fc.ferreiros.48h+0.1, fc.conny+0.05, hgnc, cex=0.3, adj=0, col=ifelse(is.na(q.conny) | q.conny > 0.1, "black", "red"), font=1+peaked))
text(1, -3.5, paste0("Removed outliers: ", paste0(outliers, collapse=", ")), adj=0, cex=0.5)

dev.off()