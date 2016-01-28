iacobucci.activated <- read.delim("/mnt/projects/ikaros/data/public/Iacobucci_2012_IKD_vs_WT_SuppTable3_upregulated.txt", stringsAsFactors = F)
iacobucci.repressed <- read.delim("/mnt/projects/ikaros/data/public/Iacobucci_2012_IKD_vs_WT_SuppTable2_downregulated.txt", stringsAsFactors = F)
iacobucci <- rbind(iacobucci.activated, iacobucci.repressed)
iacobucci <- iacobucci[!is.na(iacobucci$GeneSymbol) & iacobucci$GeneSymbol != "",]
iacobucci <- iacobucci[order(iacobucci$p),]
iacobucci <- iacobucci[!duplicated(iacobucci$GeneSymbol),]

# IKDvsIKCp
pdf("/mnt/projects/ikaros/results/iacobucci2012.pdf", height=12, width=6)
par(mfrow=c(2,1))

IKDvsIKCp <- read.delim("/mnt/projects/ikaros/results/anduril/execute/deseqAnnotated_IKD_vs_IKCp/table.csv", stringsAsFactors = F)
#IKDvsIKCp <- IKDvsIKCp[IKDvsIKCp$p < 0.1,]
m <- merge(iacobucci[,c("GeneSymbol", "p")], IKDvsIKCp[,c("Gene", "p", "fc")], by.x="GeneSymbol", by.y="Gene", suffixes=c(".iacobucci", ".conny"))
names(m)[names(m)=="fc"] <- "fc.conny"
m$logp.iacobucci <- -log10(m$p.iacobucci)
m$logp.conny <- -log10(m$p.conny)
fit <- lm(logp.conny~logp.iacobucci, m)
p <- anova(fit)$'Pr(>F)'[1]
R <- summary(fit)$r.squared
plot(m$logp.iacobucci, m$logp.conny, main=sprintf("Conny IKDvsIKCp vs. Iacobucci2012 IKDvsWT\np=%.2g, R=%.2g", p, R), xlab="-log10(p) Iacobucci", ylab="-log10(p) Conny", pch=20, cex=0.3, col=ifelse(is.na(m$fc.conny) | abs(m$fc.conny) < 1, "black", "red"))
abline(fit, col="red")
abline(h=-log10(0.05), lty=3)

labels <- m[order(rowMeans(scale(m[,c("logp.conny", "logp.iacobucci")])), decreasing=T), c("logp.conny", "logp.iacobucci", "GeneSymbol", "fc.conny")]
labels <- labels[1:50,]
text(labels$logp.iacobucci+0.02, labels$logp.conny+0.02, labels$GeneSymbol, cex=0.4, adj=0, col=ifelse(is.na(labels$fc.conny) | abs(labels$fc.conny) < 1, "black", "red"))

# IKNpvsIKCp
IKNpvsIKCp <- read.delim("/mnt/projects/ikaros/results/anduril/execute/deseqAnnotated_IKNp_vs_IKCp/table.csv", stringsAsFactors = F)
#IKNpvsIKCp <- IKNpvsIKCp[IKNpvsIKCp$p < 0.1,]
m <- merge(iacobucci[,c("GeneSymbol", "p")], IKNpvsIKCp[,c("Gene", "p", "fc")], by.x="GeneSymbol", by.y="Gene", suffixes=c(".iacobucci", ".conny"))
names(m)[names(m)=="fc"] <- "fc.conny"
m$logp.iacobucci <- -log10(m$p.iacobucci)
m$logp.conny <- -log10(m$p.conny)
fit <- lm(logp.conny~logp.iacobucci, m)
p <- anova(fit)$'Pr(>F)'[1]
R <- summary(fit)$r.squared
plot(m$logp.iacobucci, m$logp.conny, main=sprintf("Conny IKNpvsIKCp vs. Iacobucci2012 IKDvsWT\np=%.2g, R=%.2g", p, R), xlab="-log10(p) Iacobucci", ylab="-log10(p) Conny", pch=20, cex=0.3, col=ifelse(is.na(m$fc.conny) | abs(m$fc.conny) < 1, "black", "red"))
abline(fit, col="red")
abline(h=-log10(0.05), lty=3)

labels <- m[order(rowMeans(scale(m[,c("logp.conny", "logp.iacobucci")])), decreasing=T), c("logp.conny", "logp.iacobucci", "GeneSymbol", "fc.conny")]
labels <- labels[1:50,]
text(labels$logp.iacobucci+0.02, labels$logp.conny+0.02, labels$GeneSymbol, cex=0.4, adj=0, col=ifelse(is.na(labels$fc.conny) | abs(labels$fc.conny) < 1, "black", "red"))

dev.off()
