library(biomaRt)
human = useMart(biomart="ENSEMBL_MART_ENSEMBL", host="grch37.ensembl.org", path="/biomart/martservice", dataset="hsapiens_gene_ensembl") # GRCh37, v75
mouse = useMart(biomart="ENSEMBL_MART_ENSEMBL", host="feb2014.archive.ensembl.org", dataset="mmusculus_gene_ensembl") # GRCm38, v75

fc.cutoff <- 20

schwickert.activated <- read.delim("/mnt/projects/ikaros/data/public/schwickert_2014_suppl_table_1_in_vitro_ikaros_activated.txt", stringsAsFactors = F)
schwickert.repressed <- read.delim("/mnt/projects/ikaros/data/public/schwickert_2014_suppl_table_1_in_vitro_ikaros_repressed.txt", stringsAsFactors = F)
schwickert.activated$fc <- -schwickert.activated$fc
schwickert <- rbind(schwickert.activated, schwickert.repressed)

humOrt <- getLDS(attributes = c("mgi_symbol"),
                 filters = "mgi_symbol", values=c(schwickert.activated$gene, schwickert.repressed$gene), mart = mouse,
                 attributesL = c("hgnc_symbol", "entrezgene"), martL = human)
colnames(humOrt) <- c("mgi_symbol", "hgnc", "entrezgene")
humOrt <- humOrt[humOrt$hgnc != "",]

schwickert <- merge(schwickert, humOrt[,c("mgi_symbol", "hgnc")], by.x="gene", by.y="mgi_symbol")

IKDNPvsIKCDP <- read.delim("/mnt/projects/ikaros/results/anduril/execute/deseqAnnotated_IKDNPvsIKCDP/table.csv", stringsAsFactors = F)

m <- merge(schwickert, IKDNPvsIKCDP[,c("Gene", "fc", "q")], by.x="hgnc", by.y="Gene", suffixes=c(".schwickert", ".conny"))
m$fc.schwickert[m$fc.schwickert>fc.cutoff] <- fc.cutoff
m$fc.schwickert[m$fc.schwickert< -fc.cutoff] <- -fc.cutoff

plot(m$fc.schwickert, m$fc.conny, xlim=c(-fc.cutoff, fc.cutoff))
fit <- lm(fc.conny~fc.schwickert, m)
abline(fit, col="red")

