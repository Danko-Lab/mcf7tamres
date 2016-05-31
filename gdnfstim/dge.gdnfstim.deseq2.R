## get correlations between MCF-7 samples

## Read in refseq genes.
#refGene <- read.table("refGene.bed.gz")
refGene <- read.table("../annotations/tuSelecter/final_tus.txt", header=TRUE)
#refGene <- rbind(refGene, read.table("../annotations/tuSelecter/final_tus.ESR1_GREB1.txt", header=TRUE))

refGene <- refGene[grep("random|Un|hap", refGene$TXCHROM, invert=TRUE),]
refGene <- refGene[(refGene$TXEND-refGene$TXSTART)>2000,]

bodies <- refGene
bodies$TXSTART[bodies$TXSTRAND == "+"] <-bodies$TXSTART[bodies$TXSTRAND == "+"]+1000
bodies$TXEND[bodies$TXSTRAND == "-"] <- bodies$TXEND[bodies$TXSTRAND == "-"]-1000

## Adjust for gene size.
size <- refGene$V3-refGene$V2
maxsize <- 60000
bodies$TXEND[bodies$TXSTRAND == "+" & size>maxsize] <- refGene$TXSTART[bodies$TXSTRAND == "+" & size>maxsize]+maxsize
bodies$TXSTART[bodies$TXSTRAND == "-" & size>maxsize] <- refGene$TXEND[bodies$TXSTRAND == "-" & size>maxsize]-maxsize

tss <- refGene
tss$TXEND[tss$TXSTRAND == "+"] <-tss$TXSTART[tss$TXSTRAND == "+"]+1000
tss$TXSTART[tss$TXSTRAND == "-"] <- tss$TXEND[tss$TXSTRAND == "-"]-1000

## Read in bigWigs.
require(bigWig)

countBigWig <- function(prefix, path="../data/") {
 pl <- load.bigWig(paste(path, "MCF-7_", prefix, "_plus.bw", sep=""))
 mn <- load.bigWig(paste(path, "MCF-7_", prefix, "_minus.bw", sep=""))

 counts <- bed6.region.bpQuery.bigWig(pl, mn, bodies, abs.value = TRUE)
 return(counts)
}

countTimePoints <- function(prefix) {
  count0 <- countBigWig(paste(prefix, "_GDNF_0_hr", sep=""))
  count1 <- countBigWig(paste(prefix, "_GDNF_1_hr", sep=""))
  count24 <- countBigWig(paste(prefix, "_GDNF_24_hr", sep=""))

  ret_value <- cbind(count0, count1, count24)
  colnames(ret_value) <- paste(prefix, c("_0h", "_1h", "_24h"), sep="")
  return(ret_value)
}

gene_body_counts <- cbind( countTimePoints("B7"), countTimePoints("C11"), countTimePoints("G11"), countTimePoints("H9") )
save.image("Counts.RData")

## Count reads in each ...
require(DESeq2)
PVAL <- 0.01

## Build experimental design matrix
sampleID <- colnames(gene_body_counts)
tamres   <- c(rep("sen", 6), rep("res", 6))
gdnftrt  <- rep(c("gdnf0h", "gdnf1h", "gdnf24h"), 4)
data.frame(sampleID, tamres, gdnftrt)

## Create DESeq2 object.
design <- data.frame(Condition= tamres, Treatment= gdnftrt, row.names=colnames(gene_body_counts))
dds <- DESeqDataSetFromMatrix(countData= gene_body_counts, colData= design, design= ~ Condition+Treatment)
dds$Treatment <- relevel(dds$Treatment, ref="gdnf0h") ## Set the reference condition as the primary tumor.
dds <- DESeq(dds)

tam_SR <- results(dds, contrast= c("Condition", "sen", "res"))
gdnf_1h <- results(dds, contrast= c("Treatment", "gdnf1h", "gdnf0h"))
gdnf_24h <- results(dds, contrast= c("Treatment", "gdnf24h", "gdnf0h"))

sum(gdnf_1h$padj < 0.01, na.rm=TRUE)
sum(gdnf_24h$padj < 0.01, na.rm=TRUE)

############################################################################################33
## Limit analysis to sensitive || resistant lines where the effect sizes are expected to be largest.
design <- data.frame(Condition= tamres, Treatment= gdnftrt, row.names=colnames(gene_body_counts))[1:6,]
dds <- DESeqDataSetFromMatrix(countData= gene_body_counts[,1:6], colData= design, design= ~ Treatment)
dds$Treatment <- relevel(dds$Treatment, ref="gdnf0h") ## Set the reference condition as the primary tumor.
dds <- DESeq(dds)
gdnf_1h_s <- results(dds, contrast= c("Treatment", "gdnf1h", "gdnf0h"))
gdnf_24h_s <- results(dds, contrast= c("Treatment", "gdnf24h", "gdnf0h"))

sum(gdnf_1h_s$padj < 0.01, na.rm=TRUE)
sum(gdnf_24h_s$padj < 0.01, na.rm=TRUE)

design <- data.frame(Condition= tamres, Treatment= gdnftrt, row.names=colnames(gene_body_counts))[7:12,]
dds <- DESeqDataSetFromMatrix(countData= gene_body_counts[,7:12], colData= design, design= ~ Treatment)
dds$Treatment <- relevel(dds$Treatment, ref="gdnf0h") ## Set the reference condition as the primary tumor.
dds <- DESeq(dds)
gdnf_1h_r <- results(dds, contrast= c("Treatment", "gdnf1h", "gdnf0h"))
gdnf_24h_r <- results(dds, contrast= c("Treatment", "gdnf24h", "gdnf0h"))

sum(gdnf_1h_r$padj < 0.01, na.rm=TRUE)
sum(gdnf_24h_r$padj < 0.01, na.rm=TRUE)

#ss <- fitModel()
gene_pvals <- cbind(refGene, FDR_TAM= tam_SR$padj, FC_TAM= tam_SR$log2FoldChange,
                             FDR_GDNF_1h= gdnf_1h$padj, FC_GDNF_1h= gdnf_1h$log2FoldChange,
                             FDR_GDNF_24h= gdnf_24h$padj, FC_GDNF_24h= gdnf_24h$log2FoldChange,
                             FDR_GDNF_1h_S= gdnf_1h_s$padj, FC_GDNF_1h_S= gdnf_1h_s$log2FoldChange,
                             FDR_GDNF_24h_S= gdnf_24h_s$padj, FC_GDNF_24h_S= gdnf_24h_s$log2FoldChange,
                             FDR_GDNF_1h_R= gdnf_1h_r$padj, FC_GDNF_1h_R= gdnf_1h_r$log2FoldChange,
                             FDR_GDNF_24h_R= gdnf_24h_r$padj, FC_GDNF_24h_R= gdnf_24h_r$log2FoldChange
							)

write.table(gene_pvals, "GDNF_treatment.tsv", sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)

head(gene_pvals[order(gene_pvals$FDR_TAM),], 20)
head(gene_pvals[order(gene_pvals$FDR_GDNF_24h),c(7,10:NCOL(gene_pvals))], 100)
head(gene_pvals[order(gene_pvals$FDR_GDNF_1h),c(7,10:NCOL(gene_pvals))], 100)

head(gene_pvals[gene_pvals$GENENAME == "ESR1",c(7,10:NCOL(gene_pvals))], 100)
head(gene_pvals[gene_pvals$GENENAME == "CCND1",c(7,10:NCOL(gene_pvals))], 100)
head(gene_pvals[gene_pvals$GENENAME == "GREB1",c(7,10:NCOL(gene_pvals))], 100)
head(gene_pvals[gene_pvals$GENENAME == "HOXC13",c(7,10:NCOL(gene_pvals))], 100)


head(gene_pvals[gene_pvals$V7 == "AKT",c(7,10:NCOL(gene_pvals))], 100)

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

#gene_pvals[which(gdnf_24h_s$log2FoldChange > 2.5 & gdnf_24h_s$baseMean > 0 & !is.na(gdnf_24h_s$log2FoldChange)),]
laball_1h <- c("EGR1", "EGR2", "FOSL1", "LRRC15", "EGR3", "NR4A3", "ETS2", "ESR1", "NDRG3", "GPER1")
laball_24h <- c("SHC4", "PLD1", "ERRFI1", "RASGEF1A", "LIPK", "VAV3", "ESR1", "PGR", "GATA4", "RET", "ACACB", "BRIP1")

WriteMAPlot <- function(res, laball, prefix, ylim=c(-3,3)) {
## MA-plot
pdf(paste("MAPlot.deSeq2.",prefix,".pdf", sep=""))
 plotMA(res, xlim=c(1e-1, 2e5), ylim=ylim, alpha= 0.01, cex=0.75, ylab= "Log Fold-Change, TamS/ TamR")

 plotMA(res, xlim=c(1e-1, 2e5), ylim=ylim, alpha= 0.01, cex=0.75, ylab= "Log Fold-Change, TamS/ TamR")
 addlab(laball, res, bodies)
 
 plot(0, 0, ylim=ylim, xlim=c(1e-1, 2e5), log="x")
 addlab(laball, res, bodies)
dev.off()

}

WriteMAPlot(res= gdnf_1h_s, laball_1h, "1h_s", ylim=c(-2.5,6))
WriteMAPlot(res= gdnf_24h_s, laball_24h, "24h_s", ylim=c(-3,3))

WriteMAPlot(res= gdnf_1h, laball_1h, "1h_all", ylim=c(-2,3))
WriteMAPlot(res= gdnf_24h, laball_24h, "24h_all", ylim=c(-2,2))


## Test correlation btwn. TamS and TamR.
cor.test(gene_pvals$FC_GDNF_1h_S, gene_pvals$FC_GDNF_1h_R)
cor.test(gene_pvals$FC_GDNF_24h_S, gene_pvals$FC_GDNF_24h_R)

plot(gene_pvals$FC_GDNF_1h_S, gene_pvals$FC_GDNF_1h_R)
