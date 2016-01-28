options("help_type"="html")
options(stringsAsFactors=F)

getwd()
rm(list=ls(all=TRUE))

pathA <- "/mnt/projects/ikaros/results/truus_microarray"
pathC <- "/mnt/projects/ikaros/data/truus_microarray"
pathGSEA <- "/home/STANNANET/maximilian.kauer/DATA/amnesia/AnnaMaria/Analysis/GSEA"
pathMot <- "/home/STANNANET/maximilian.kauer/DATA/amnesia/AnnaMaria/Analysis/motifSearch"
pathGenome <- "/home/STANNANET/maximilian.kauer/DATA/amnesia/genome/motifs"

pathFct <-  "/home/STANNANET/maximilian.kauer/DATA/amnesia/R.Functions.etc"
pathAnnot <- "/home/STANNANET/maximilian.kauer/DATA/amnesia/genome/AnnotationData"  
pathGSEAFct <- "/home/STANNANET/maximilian.kauer/DATA/amnesia/R.Functions.etc/PGSEA"

pathRaft1 <- "/home/STANNANET/maximilian.kauer/DATA/amnesia/Lab3.MISC/Raftlin/Analysis"
pathRaft2 <- "/home/STANNANET/maximilian.kauer/DATA/amnesia/Lab3.MISC/Raftlin/NEU/Analysis"

pathTcn <- "/home/STANNANET/maximilian.kauer/DATA/amnesia/timecourse.new/Analysis"
pathMisc <- "/home/STANNANET/maximilian.kauer/DATA/amnesia/Lab3.MISC"


"%i%" <- intersect

load(paste(pathFct, "SLmisc", sep="/"))
load(paste(pathFct, "max.func.env", sep="/"))
source(paste(pathFct, "usefulFunctions.R", sep="/"))

"%notin%" <- function(x, table) match(x, table, nomatch = 0) == 0

source(paste(pathFct, "load.packages.R", sep="/"))

source("https://bioconductor.org/biocLite.R")
biocLite("hgu133plus2.db")
biocLite("affy")
biocLite("arrayQualityMetrics")
biocLite("gcrma")
biocLite("frma")
biocLite("hgu133plus2frmavecs")
biocLite("genefilter")
library(hgu133plus2probe)
library(hgu133plus2cdf)
library(affy)
library(arrayQualityMetrics)
library(gcrma)
library(frma)
library(hgu133plus2frmavecs)
library(genefilter)

NAME <- "TruusIKZF1"
db <- "hgu133plus2.db"

##############################################################
## Affybatch
##############################################################

pD <- read.csv("/mnt/projects/ikaros/data/truus_microarray/sample_annotation.tsv", sep="\t", row.names = 1); pD

setwd(pathC)
Data <- ReadAffy(filenames=as.character(rownames(pD)), compress=TRUE)
sampleNames(Data) <- as.vector(pD$id)
rownames(pD) <- as.vector(pD$id)

pD <- as(pD, "AnnotatedDataFrame")
phenoData(Data) <- pD

save(Data, file=paste0("/mnt/projects/ikaros/results/truus_microarray/", NAME, ".affybatch.Rdata"), compress=TRUE)
#load(paste(NAME,".affybatch",sep=""))

#### arrayQualityMetrics
setwd(pathA)

# NOTE: fails with error: Error in file.info(x) : invalid filename argument
# likely because of missing/outdated dependency; run on biowaste instead (where it works)
arrayQualityMetrics(Data, spatial=F, intgroup="IKZF1_status", outdir="/mnt/projects/ikaros/results/truus_microarray/qc", force=T)

## eset
eset_gcrma <- gcrma(Data, fast=T); pData(eset_gcrma)
sampleNames(eset_gcrma) <- as.vector(Data$id)
save(eset_gcrma, file=paste0("/mnt/projects/ikaros/results/truus_microarray/", NAME, ".eset.gcrma.Rdata"), compress=TRUE)

##########################################################



##########################################################
## fRMA
##########################################################
setwd(pathA)
load(paste0("/mnt/projects/ikaros/results/truus_microarray/", NAME, ".affybatch.Rdata"))

library(frma)
esetfrma <- frma(Data, summarize="robust_weighted_average") ; sumName <- "r_w_a"
phenoData(esetfrma) <- phenoData(Data)
sampleNames(esetfrma) <- as.vector(Data$id)

save(esetfrma, file=paste0("/mnt/projects/ikaros/results/truus_microarray/", NAME, ".eset.frma.Rdata"), compress=TRUE)

bc.bin <- barcode(esetfrma, output="binary")
save(bc.bin, file=paste0("/mnt/projects/ikaros/results/truus_microarray/", NAME, ".eset.frma.barcode_bin.Rdata"), compress=TRUE)

bc.z <- barcode(esetfrma, output="z-score")
save(bc.z, file=paste0("/mnt/projects/ikaros/results/truus_microarray/", NAME, ".eset.frma.barcode_zScore.Rdata"), compress=TRUE)

par(mfrow=c(3,4))
for(i in 1:ncol(bc.z)) {
	plot(density(exprs(esetfrma)[which(bc.z[,i] > 0) ,i]))
}

##########################################################



## some dotplots #########################################
ESET <- eset_nsF ; esetName <- "frmansF"
#ESET <- esetfrma ; esetName <- "frma"

allsyms <- getSYMBOL( rownames(exprs(ESET) ), db ) 
symList <- sort( unique( c("MKL1", "MKL2", "TAGLN", "NKX2-2", "NR0B1", "EWSR1", "FLI1", "HEY1", "DKK1", "DKK2", "CXCR4", "TGFB1","TGFBR2", "IGF1", "IGFBP3", "IGFBP5", "IGFBP7", "CYR61", "PTGER3", "CALD1","NR0B1","NKX2.2", "CUL1", "STEAP1",  "MAPT", "MDM4", "NOTCH1", "TP53", "TERT", "PDGFC", "UPP1", "CAV1", "PDGFC", "PLD2", "ID2", "POLA1", "SKP2", "TYMS", "E2F1", "E2F2", "E2F3", "E2F4", "E2F5", "E2F6", "E2F7", "E2F8"))); listName <- "someTestGenes"  
q <- vector(length=length(symList)) 
for (i in 1:length(q)) {q[i] <- paste("^",symList[i], "$", sep="")}; q 

pdf(file=paste("dotcharts", esetName, listName, "pdf", sep="."), height=12, width=8.5)
for(i in 1:length(symList)) {
	wh1 <- grep( q[i], allsyms ); wh1
	if( length( wh1 ) > 0 ) {
		cols <- rainbow(length(wh1)); if(length(cols)==2) { cols <- c("red", "blue") } 
		dotchart( exprs(ESET)[wh1,], cex=1, col=cols, pch=20, main=q[i] )
	}
}
graphics.off()
##########################################################


##########################################################
## simple plots
##########################################################
sampleNames(Data) <- as.vector(Data$id)
cols <- rainbow(length(Data$id))[ as.factor(Data$id) ]

pdf(file=paste("/mnt/projects/ikaros/results/truus_microarray/", NAME, ".QC.bp.hist.pdf", sep=""), width=12, height=9)

par(mfrow=c(1,2)); par(mar=c(10,2,2,2))
boxplot(Data,col=cols, names=Data$id, las=3, cex=0.4, main="raw - affybatch")
hist(Data,col=cols,lty=1,xlab="Log (base 2) intensities", main="raw - affybatch")
legend(12,0.25, Data$id,cex=0.6,lty=1,col=cols)

#boxplot(eset_gcrma, col=cols, names=eset_gcrma$id, las=3, cex=0.4, main="eset_gcrma")
#hist(eset_gcrma, col=cols,lty=1,xlab="Log (base 2) intensities", main="eset_gcrma")
#legend(12,0.4, eset_gcrma$id, cex=0.6,lty=1,col=cols)

boxplot(esetfrma,col=cols, names=esetfrma$id, las=3, cex=0.4, main="eset_frma")
hist(esetfrma,col=cols,lty=1,xlab="Log (base 2) intensities", main="eset_frma")
legend(12,0.4, esetfrma$id,cex=0.6,lty=1,col=cols)

graphics.off()

smoothScatter( exprs(eset_gcrma)[,1], exprs(esetfrma)[,1] )
smoothScatter( exprs(esetfrma)[1:54613,1], bc.z[,1] )
##########################################################



firstesetName <- "frmaeset"; eset <- esetfrma 

##########################################################
## barcode  filter
##########################################################
ESET <- eset; esetName <- firstesetName
#ESET <- eset_gcrma; esetName <- "gcrma"


### filter out all under threshold: z_th
z_th <- 2
whlow <- as.vector(which(apply(bc.z, 1, function(x) all(x < z_th )))); head(whlow)
bc.z_f <- bc.z[-whlow,]; dim(bc.z_f)

whhigh <- which(rownames(exprs(ESET)) %in% rownames(bc.z_f) ); length(whhigh)
eset_f <- ESET[whhigh,]; dim(eset_f)
#######################


### nsF Filter ########
nsF <-  nsFilter(eset_f, require.entrez = T, require.GOBP = FALSE,
		require.GOCC = FALSE, require.GOMF = FALSE,
		remove.dupEntrez = T, var.func = sd, var.cutoff = 0.00001,
		var.filter = TRUE, feature.exclude="^AFFX")
eset_nsF <- nsF$eset; dim(exprs(eset_nsF))
rm(nsF)
##########################################################



#ESET <- eset; esetName <- firstesetName
ESET <- eset_nsF; esetName <- "frmaeset_nsF"
#ESET <- eset_nsF; esetName <- "gcrmaeset_nsF"

head( exprs(ESET) )
###################################################################



######## sd test ##################################################
Sd <- sort(apply(exprs(ESET), 1, sd), decreasing=T); head(Sd)
###################################################################



#### correlationof samples ########################################
howMany <- c(dim(ESET)[[1]], 5000, 1000, 500); 
pdf(paste("cor.heatmap", esetName,"zth",z_th,"pdf",sep="."), height=9, width=12) 
cols <- redblue(256) # colorpanel(100, "honeydew", "yellow2","tomato") 
Col <- rainbow(length(Data$IKZF1_status))[as.factor(Data$IKZF1_status)]
#layout(matrix(c(1,2,3,4),2,2,byrow=T)); #layout.show(5)
for(i in howMany) {
	ps <- names(Sd)[1:i]
	Cor <- cor(exprs(ESET)[ps,], method="pearson")   
	hmat <- Cor    
	heatmap.2(as.matrix(hmat), main=paste(esetName, "zth",z_th,i, "Genes"), distfun=function(x)  dist(x, method = "euclidean"), hclustfun=function(m) hclust(m, method="average"), 
			colsep=seq(1:dim(hmat)[[2]]), rowsep=seq(1:dim(hmat)[[1]]), sepcolor="grey92", sepwidth=c(0.005,0.005),
			ColSideColors=Col, 
			keysize=0.8, col=cols, density.info="none", symkey=FALSE, labRow=rownames(hmat), labCol=colnames(hmat), trace="none",cex.main=0.7, cexRow=1, cexCol=1, mar=c(15,15) )
	
}
graphics.off()
###################################################################



## limma to test MKLkd vs. ns in 4 conditions #####################
setwd( pathA )

## wie in UserGuide "9.5.2 Analysing as for a Single Factor"
allcond <- as.factor( paste( ESET$treat2, ESET$treat, ESET$treat1, sep="_" ) ); allcond
levels( allcond )

design <- model.matrix( ~0+allcond )
colnames( design ) <- levels( allcond ); design
fit <- lmFit( ESET, design )
heatmap( cor(fit$coefficients), mar=c(15,15) )

contrast.matrix <- makeContrasts(                  
		
		EFkd       = doxy_starved_shNS-con_starved_shNS,
		MKLkd_d_st = doxy_starved_shMKL1_2-doxy_starved_shNS,
		MKLkd_c_st = con_starved_shMKL1_2-con_starved_shNS,
		MKLkd_d_si = doxy_SI_shMKL1_2-doxy_SI_shNS,
		MKLkd_c_si = con_SI_shMKL1_2-con_SI_shNS,
		
		SI_c_ns    = con_SI_shNS-con_starved_shNS,
		SI_c_mkkd  = con_SI_shMKL1_2-con_starved_shMKL1_2,
		SI_d_ns    = doxy_SI_shNS-doxy_starved_shNS,
		SI_d_mkkd  = doxy_SI_shMKL1_2-doxy_starved_shMKL1_2,
		
		ddiff_SI_c_mkkd_VS_SI_c_ns = (con_SI_shMKL1_2-con_starved_shMKL1_2)-(con_SI_shNS-con_starved_shNS),
		ddiff_SI_d_mkkd_VS_SI_d_ns = (doxy_SI_shMKL1_2-doxy_starved_shMKL1_2)-(doxy_SI_shNS-doxy_starved_shNS),
		
		ddiff_MKLkd_d_si_VS_MKLkd_c_si = (doxy_SI_shMKL1_2-doxy_SI_shNS)-(con_SI_shMKL1_2-con_SI_shNS),
		ddiff_MKLkd_d_st_VS_MKLkd_c_st = (doxy_starved_shMKL1_2-doxy_starved_shNS)-(con_starved_shMKL1_2-con_starved_shNS),                                                                      
		levels=design )
contrast.matrix
fit2 <- contrasts.fit( fit, contrast.matrix )
fit2 <- eBayes(fit2) 
head( fit2$coefficients)
head( fit2$p.value )
############################


############################
colnames( fit2$coefficients ) <- gsub(" - ", "_vs_", colnames(fit2$coefficients) )
colnames( fit2$p.value ) <- paste( "Pvalue", gsub(" - ", "_vs_", colnames( fit2$p.value)),sep="." )

## adjust P ###
adjP <- fit2$p.value; 
colnames(adjP) <- gsub("Pvalue", "adjP",colnames(fit2$p.value))
for (j in 1:ncol(adjP)) { adjP[, j] <- p.adjust(fit2$p.value[, j], method = "BH") }

P=1; L=1.5


## venn
pdf( file=paste("vennD", "Pth", P, "logFC", L, esetName,"zth",z_th, "pdf", sep="."), height=9, width=9 )    
results <- decideTests( fit2, p.value=P, lfc=L ); 

par( mfrow=c(2,2))
#P=0.05; L=0.5; results <- decideTests( fit2, p.value=P, lfc=L ); vennDiagram( results, cex=0.65, main=(paste("p-value:",P, "logFC",L)) )
ind <- combn(length( colnames( results)), 2)
for( i in 1:ncol(ind) ) {
	vennDiagram( results[, ind[,i]], cex=0.5, main=( paste("p-value:",P, "logFC",L) ) )
}
#P=0.01; L=1.5; results <- decideTests( fit2, p.value=P, lfc=L ); vennDiagram( results, cex=0.65, main=(paste("p-value:",P, "logFC",L)) )
graphics.off()
############################


## hists of logFC ##########
pdf( file=paste("hists.of.logFCs", esetName,"zth",z_th, "pdf", sep="."), height=9, width=9 )  
par(mfrow=c(3,3))
for (j in 1:ncol(fit2$coefficients)) {
	hist( fit2$coefficients[,j], main=colnames(fit2$coefficients)[j], breaks=100 )
}
graphics.off()

pdf(file=paste("hists.of.PVals", esetName,"zth",z_th, "pdf", sep="."), height=9, width=9)
par(mfrow=c(3,3))
for (j in 1:ncol(fit2$p.value)) {
	hist(fit2$p.value[,j], main=colnames(fit2$p.value)[j], breaks=100)
}
graphics.off()

pdf(file=paste("hists.of.adjPVals", esetName,"zth",z_th, "pdf", sep="."), height=9, width=9)
par(mfrow=c(3,3))
for (j in 1:ncol(adjP)) {
	hist( adjP[,j], main=colnames(adjP)[j], breaks=100 )
}
graphics.off()
####################################################################




##### annotate/write out ##########################################
stopifnot(identical(rownames(fit$coefficients), rownames(fit2$coefficients)) ) # muss TRUE sein!!!
stopifnot(identical(rownames(fit2$p.value), rownames(adjP)) ) # muss TRUE sein!!!
stopifnot(identical(rownames(fit2$p.value), rownames( exprs(ESET) )) ) # muss TRUE sein!!!

matAnn <- as.data.frame( cbind( round(fit2$coefficients,2), round(adjP,4), round(fit$coefficients,2), round(fit2$p.value,4), round(exprs(ESET),2) ), stringsAsFactors=F); head(matAnn)

matAnn$affID <- rownames(matAnn)
matAnn$EG <- getEG( rownames(matAnn), db )
matAnn$syms <- getSYMBOL( rownames(matAnn), db )
matAnn <- matAnn[, c(ncol(matAnn), (ncol(matAnn)-1), (ncol(matAnn)-2),  1:(ncol(matAnn)-3)) ]; head(matAnn)

matAnn <- matAnn[ order(matAnn[, 12] ), ]
save(matAnn, file=paste("matAnn",NAME, esetName,"zth",z_th, "annot.RObj", sep=".") , compress=T)
write.table(matAnn, file=paste("matAnn",NAME, esetName, "zth", z_th, "annot.xls", sep="."), quote=F, row.names=F, sep="\t")
###################################################################




####################################################################
mat <- fit2$coefficients; mName <- "logFCs"; head( mat );
x <- c(  "SI_c_ns",    "SI_c_mkkd",    "MKLkd_c_st",    "MKLkd_c_si",    "SI_d_ns",    "SI_d_mkkd",    "MKLkd_d_st",    "MKLkd_d_si" , "EFkd", "EFkd" )
y <- c(  "SI_ns_u2os", "SI_mkkd_u2os", "MKLkd_st_u20s", "MKLkd_si_u2os", "SI_ns_u2os", "SI_mkkd_u2os", "MKLkd_st_u20s", "MKLkd_si_u2os","MKLkd_d_st","MKLkd_st_u20s" )

round( cor( mat[ , c( x, y) ] ) ,2)

lth <- 1.3

#gn <- "withGeneNames"
gn <- "withoutGeneNames"

pdf( file=paste("someScatterOfMostInterest", mName, esetName, gn, "zth",z_th, "pdf", sep="."), height=9, width=9 )
par( mfrow=c(2,2))
for( i in 1:( length(x) ) ) {    
	smoothScatter( mat[, x[i] ], mat[, y[i] ], xlab=x[i], ylab=y[i] , main=paste("R=", round( cor( mat[, x[i] ], mat[, y[i] ] ), 2) ) )
	abline( h=0, v=0 , col="grey")
	abline( a=0, b=1 , col="green", lty=2 )
	if( gn == "withGeneNames" ) { 
		whsig <- unique( c(which( abs(mat[, x[i] ]) > lth ), which( abs(mat[, y[i] ]) > lth )) )
		if( length(whsig)>0) { text( mat[ whsig,i], mat[whsig,j], labels=getSYMBOL( rownames(mat)[whsig], db ), cex=0.6, col="tomato" ) }
	}
}
graphics.off()      
###################################################################


