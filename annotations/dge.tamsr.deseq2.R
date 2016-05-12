## get correlations between MCF-7 samples

## Read in refseq genes.
#refGene <- read.table("refGene.bed.gz")
refGene <- read.table("tuSelecter/final_tus.txt", header=TRUE)
#refGene <- rbind(refGene, read.table("../annotations/tuSelecter/final_tus.ESR1_GREB1.txt", header=TRUE))

refGene <- refGene[grep("random|Un|hap", refGene$TXCHROM, invert=TRUE),]
refGene <- refGene[(refGene$TXEND-refGene$TXSTART)>4000,]

bodies <- refGene
bodies$TXSTART[bodies$TXSTRAND == "+"] <-bodies$TXSTART[bodies$TXSTRAND == "+"]+1000
bodies$TXEND[bodies$TXSTRAND == "-"] <- bodies$TXEND[bodies$TXSTRAND == "-"]-1000

tss <- refGene
tss$TXEND[tss$TXSTRAND == "+"] <-tss$TXSTART[tss$TXSTRAND == "+"]+1000
tss$TXSTART[tss$TXSTRAND == "-"] <- tss$TXEND[tss$TXSTRAND == "-"]-1000

## Read in bigWigs.
require(bigWig)

path="../data/"
B7_pl <- load.bigWig(paste(path, "sB7_pl.bw", sep=""))
B7_mn <- load.bigWig(paste(path, "sB7_mn.bw", sep=""))
G11_pl <- load.bigWig(paste(path, "rG11_pl.bw", sep=""))
G11_mn <- load.bigWig(paste(path, "rG11_mn.bw", sep=""))
H9_pl <- load.bigWig(paste(path, "rH9_pl.bw", sep=""))
H9_mn <- load.bigWig(paste(path, "rH9_mn.bw", sep=""))
C11_pl <- load.bigWig(paste(path, "sC11_pl.bw", sep=""))
C11_mn <- load.bigWig(paste(path, "sC11_mn.bw", sep=""))

B7_BBCA_pl <- load.bigWig(paste(path,"sB7_bbca_pl.bw", sep=""))
B7_BBCA_mn <- load.bigWig(paste(path,"sB7_bbca_mn.bw", sep=""))
G11_BBCA_pl <- load.bigWig(paste(path,"rG11_bbca_pl.bw", sep=""))
G11_BBCA_mn <- load.bigWig(paste(path,"rG11_bbca_mn.bw", sep=""))

## Count reads in each ...
B7 <- bed6.region.bpQuery.bigWig(B7_pl, B7_mn, bodies, abs.value = TRUE)
G11<- bed6.region.bpQuery.bigWig(G11_pl, G11_mn, bodies, abs.value = TRUE)
H9 <- bed6.region.bpQuery.bigWig(H9_pl, H9_mn, bodies, abs.value = TRUE)
C11<- bed6.region.bpQuery.bigWig(C11_pl, C11_mn, bodies, abs.value = TRUE)
B7b <- bed6.region.bpQuery.bigWig(B7_BBCA_pl, B7_BBCA_mn, bodies, abs.value = TRUE)
G11b<- bed6.region.bpQuery.bigWig(G11_BBCA_pl, G11_BBCA_mn, bodies, abs.value = TRUE)
gene_body_counts <- cbind(B7, C11, H9, G11, B7b, G11b)

save.image("Counts.RData")

## 
require(DESeq2)
PVAL <- 0.01

## Build experimental design matrix
sampleID <- colnames(gene_body_counts)
tamres   <- c("sen", "sen", "res", "res", "sen", "res")
bbcatrt  <- c("u", "u", "u", "u", "t", "t")
data.frame(sampleID, tamres, bbcatrt)

## Create DESeq2 object.
design <- data.frame(Condition= tamres, Treatment= bbcatrt, row.names=colnames(gene_body_counts))
dds <- DESeqDataSetFromMatrix(countData= gene_body_counts, colData= design, design= ~ Condition)
dds$Condition <- relevel(dds$Condition, ref="res") ## Set the reference condition as the primary tumor.
dds <- DESeq(dds)
res <- results(dds)

## Fit bbca.
ddsB <- DESeqDataSetFromMatrix(countData= gene_body_counts, colData= design, design= ~ Treatment)
ddsB$Treatment <- relevel(ddsB$Treatment, ref="u") ## Set the reference condition as the primary tumor.
ddsB <- DESeq(ddsB)
resB <- results(ddsB)

#ss <- fitModel()
gene_pvals <- cbind(refGene, FDR_TAM= res$padj, FC_TAM= res$log2FoldChange,
                             FDR_BBCA= resB$padj, FC_BBCA= resB$log2FoldChange)


## Add specific genes to the plot.
addlab <- function(gene_ID, deRes, genes, ...) {
 idx<-sapply(gene_ID, function(gene_ID) {
  ig <- which(genes[,7] == gene_ID)
  io <- ig[which.min(deRes$padj[ig])]
  if(NROW(io)>0) {
        text(deRes$baseMean[io], deRes$log2FoldChange[io], labels= gene_ID, cex= 1, pos= 3, ...)
        points(deRes$baseMean[io], deRes$log2FoldChange[io], col="blue", cex=1.5)
        io
  }
 })
 print(idx)
 idx <- unlist(idx); idx <- idx[!is.null(idx)]
 return(data.frame(Gene= genes[idx,7], AveExpr= deRes$baseMean[idx], logFC= deRes$log2FoldChange[idx], adj.P.Val= deRes$padj[idx]))
}

## MA-plot
pdf("MAPlot.deSeq2.pdf")
plotMA(res, ylim=c(-9,9), xlim=c(1e-1, 2e5), alpha= 0.01, cex=0.75, ylab= "Log Fold-Change, TamS/ TamR")

plotMA(res, ylim=c(-9,9), xlim=c(1e-1, 2e5), alpha= 0.01, cex=0.75, ylab= "Log Fold-Change, TamS/ TamR")
addlab(c("GREB1", "PGR", "ELOVL2", "NOS1AP"), res, bodies) ## E2 REG (Hah et. al. || classical)
addlab(c("BMP5", "C10orf112", "PCGEM1"), res, bodies) ## Non E2 REG
addlab(c("GDNF", "MPPED1", "BANK1", "PPP3CA", "RP11-369C8.1", "RP4-754E20__A.5", "FOLH1", "SIGLEC6", "KCNJ8"), res, bodies) ## Down reg

plot(0, 0, ylim=c(-9,9), xlim=c(1e-1, 2e5), log="x")
addlab(c("GREB1", "PGR", "ELOVL2", "NOS1AP"), res, bodies) ## E2 REG (Hah et. al. || classical)
addlab(c("BMP5", "C10orf112", "PCGEM1"), res, bodies) ## Non E2 REG
addlab(c("GDNF", "MPPED1", "BANK1", "PPP3CA", "RP11-369C8.1", "RP4-754E20__A.5", "FOLH1", "SIGLEC6", "KCNJ8"), res, bodies) ## Down reg

dev.off()

head(gene_pvals[order(gene_pvals$FDR_TAM),c(1:3,6:9,14:NCOL(gene_pvals))], n=30)
head(gene_pvals[order(gene_pvals$FDR_BBCA),c(1:3,6:7,10:11)], n=30)
head(gene_pvals[order(gene_pvals$FDR_BBCA_G11),c(1:3,6:7,12:13)], n=30)

## Number of genes w/ p< 0.01
NROW(unique(gene_pvals$GENEID[gene_pvals$FDR_TAM < 0.01]))
NROW(unique(gene_pvals$GENEID[gene_pvals$FDR_BBCA < 0.01]))

write.table(gene_pvals[order(gene_pvals$FDR_TAM),c(1:4,6:8,14:NCOL(gene_pvals))], "Signif.Changes.TamRes.tsv", row.names=FALSE, col.names=TRUE, quote=FALSE, sep="\t")

## Dotplots.


## SANITY CHECKS.
sanityCheck <- function(gene) {
	indx <- grep(gene, refGene$V7)
	print(data.frame(MGID= as.character(refGene[indx,7]), gene_body_counts[indx,]))
}
sanityCheck("CD36")
sanityCheck("MPPED")
sanityCheck("ESR1")



