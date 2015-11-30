e <- read.delim("/mnt/projects/ikaros/results/anduril/execute/deseqAnnotated_IKNp_vs_IKCp/table.csv")
e$mean <- pmin(e$meanA, e$meanB)
plot(log10(e$mean), -log10(e$q), cex=0.2)

library(ggplot2)
minp <- 0.1
p <- ggplot(e[!is.na(e$p) & e$p <= minp & e$mean > 0,], aes(x=log10(mean), y=-log10(p))) +
  theme_bw() +
  geom_point() + 
  stat_smooth(se=TRUE, method="loess", span = 1, size = 1.2, colour = "orange") + 
  stat_smooth(se=TRUE, method="lm", size = 1.2, colour = "blue") + 
  scale_y_continuous(limits = c(-log10(minp), 10))
  ggtitle("Read count vs. significance")
print(p)
