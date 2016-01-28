vitanza.activated <- read.delim("/mnt/projects/ikaros/data/public/Vitanza_2014_IKD_vs_WT_SuppTable1_upregulated.txt", stringsAsFactors = F)
vitanza.repressed <- read.delim("/mnt/projects/ikaros/data/public/Vitanza_2014_IKD_vs_WT_SuppTable1_downregulated.txt", stringsAsFactors = F)
vitanza <- rbind(vitanza.activated, vitanza.repressed)
vitanza <- vitanza[!is.na(vitanza$GeneSymbol) & vitanza$GeneSymbol != "",]
vitanza <- vitanza[order(abs(vitanza$fc), decreasing=T),]
vitanza <- vitanza[!duplicated(vitanza$GeneSymbol),]

pdf("/mnt/projects/ikaros/results/vitanza2014.pdf", height=12, width=6)
par(mfrow=c(2,1))

# IKDvsIKCp
IKDvsIKCp <- read.delim("/mnt/projects/ikaros/results/anduril/execute/deseqAnnotated_IKD_vs_IKCp/table.csv", stringsAsFactors = F)
#IKDvsIKCp <- IKDvsIKCp[IKDvsIKCp$p < 0.1,]
m <- merge(vitanza[,c("GeneSymbol", "fc", "p")], IKDvsIKCp[,c("Gene", "fc", "q")], by.x="GeneSymbol", by.y="Gene", suffixes=c(".vitanza", ".conny"))
names(m)[names(m)=="q"] <- "q.conny"
fit <- lm(fc.conny~fc.vitanza, m)
p <- anova(fit)$'Pr(>F)'[1]
R <- summary(fit)$r.squared
plot(m$fc.vitanza, m$fc.conny, main=sprintf("Conny IKDvsIKCp vs. Vitanza2014 IKDvsWT\np=%.2g, R=%.2g", p, R), xlab="FC Vitanza", ylab="FC Conny", pch=20, cex=0.3, col=ifelse(is.na(m$q.conny) | m$q.conny > 0.1, "black", "red"))
abline(fit, col="red")

labels <- m[order(rowMeans(scale(m[,c("fc.vitanza", "fc.conny")])), decreasing=F), c("fc.vitanza", "fc.conny", "q.conny", "GeneSymbol")]
with(labels[1:20,], text(fc.vitanza+0.05, fc.conny+0.05, GeneSymbol, cex=0.4, adj=0))
with(labels[(nrow(labels)-20):nrow(labels),], text(fc.vitanza+0.05, fc.conny+0.05, GeneSymbol, cex=0.3, adj=0, col=ifelse(is.na(q.conny) | q.conny > 0.1, "black", "red")))

# IKNpvsIKCp
IKNpvsIKCp <- read.delim("/mnt/projects/ikaros/results/anduril/execute/deseqAnnotated_IKNp_vs_IKCp/table.csv", stringsAsFactors = F)
#IKNpvsIKCp <- IKNpvsIKCp[IKNpvsIKCp$p < 0.1,]
m <- merge(vitanza[,c("GeneSymbol", "fc", "p")], IKNpvsIKCp[,c("Gene", "fc", "q")], by.x="GeneSymbol", by.y="Gene", suffixes=c(".vitanza", ".conny"))
names(m)[names(m)=="q"] <- "q.conny"
fit <- lm(fc.conny~fc.vitanza, m)
p <- anova(fit)$'Pr(>F)'[1]
R <- summary(fit)$r.squared
plot(m$fc.vitanza, m$fc.conny, main=sprintf("Conny IKNpvsIKCp vs. Vitanza2014 IKDvsWT\np=%.2g, R=%.2g", p, R), xlab="FC Vitanza", ylab="FC Conny", pch=20, cex=0.3, col=ifelse(is.na(m$q.conny) | m$q.conny > 0.1, "black", "red"))
abline(fit, col="red")

labels <- m[order(rowMeans(scale(m[,c("fc.vitanza", "fc.conny")])), decreasing=F), c("fc.vitanza", "fc.conny", "q.conny", "GeneSymbol")]
with(labels[1:20,], text(fc.vitanza+0.05, fc.conny+0.05, GeneSymbol, cex=0.4, adj=0))
with(labels[(nrow(labels)-20):nrow(labels),], text(fc.vitanza+0.05, fc.conny+0.05, GeneSymbol, cex=0.3, adj=0, col=ifelse(is.na(q.conny) | q.conny > 0.1, "black", "red")))

dev.off()
