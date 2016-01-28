m <- read.delim("/mnt/projects/ikaros/results/anduril/execute/deseqExprMatrix/table.csv", stringsAsFactors = F, row.names = 1)

crlf2 <- as.data.frame(t(m["ENSG00000205755",])) ; names(crlf2) <- "value"
crlf2$group <- NA
crlf2$group[rownames(crlf2) %in% c("IK6_6", "IKD_5")] <- "WT rel"
crlf2$group[rownames(crlf2) %in% c("IK6_2", "IKC_1", "IKC_2", "IKC_3", "IKC_4", "IKC_5", "IKC_8", "IKC_9", "IKC_10", "IKC_11", "IKC_12", "IKC_13", "IKD_2", "IKD_3")] <- "P/C dia"
crlf2$group[rownames(crlf2) %in% c("IK6_1", "IK6_4", "IK6_5", "IKC_6", "IKC_7", "IKD_1", "IKD_4")] <- "P/C rel"

pdf("/mnt/projects/ikaros/results/crlf2-expression.pdf")
boxplot(value~group, data=crlf2, ylab="Normalized expression", main="CRLF2 expression", na.action=na.exclude, outline=F, cex.axis=1)
stripchart(value~group, data=crlf2, method="jitter", na.action=na.exclude, vertical=T, pch=19, col=1:length(levels(as.factor(crlf2$group))), add=T)
dev.off()