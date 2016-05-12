## get correlations between MCF-7 samples

## Read in refseq genes.
dreg  <-read.table("mcf7_tamSR.bedgraph.bed.gz")
dregHD<-read.table("dregHD/mcf7.dREG_HD.bed"); dregHD <- dregHD[dregHD$V3-dregHD$V2 > 0,]

tres <- dregHD

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

countTRE <- function(plus, minus, tres) {
  ## For dREG HD
  bed.region.bpQuery.bigWig(plus, tres, abs.value = TRUE)+bed.region.bpQuery.bigWig(minus, tres, abs.value = TRUE)
}

## Count reads in each ...
B7 <- countTRE(B7_pl, B7_mn, tres)
G11<- countTRE(G11_pl, G11_mn, tres)
H9 <- countTRE(H9_pl, H9_mn, tres)
C11<- countTRE(C11_pl, C11_mn, tres)
B7b <- countTRE(B7_BBCA_pl, B7_BBCA_mn, tres)
G11b<- countTRE(G11_BBCA_pl, G11_BBCA_mn, tres)
tre_counts <- cbind(B7, C11, H9, G11, B7b, G11b)

## 
require(DESeq2)
PVAL <- 0.01

## Build experimental design matrix
sampleID <- colnames(tre_counts)
tamres   <- c("sen", "sen", "res", "res", "sen", "res")
bbcatrt  <- c("u", "u", "u", "u", "t", "t")
data.frame(sampleID, tamres, bbcatrt)

## Create DESeq2 object.
design <- data.frame(Condition= tamres, Treatment= bbcatrt, row.names=colnames(tre_counts))
dds <- DESeqDataSetFromMatrix(countData= tre_counts, colData= design, design= ~ Condition)
dds$Condition <- relevel(dds$Condition, ref="res") ## Set the reference condition as the primary tumor.
dds <- DESeq(dds)
res <- results(dds)

## Fit bbca.
ddsB <- DESeqDataSetFromMatrix(countData= tre_counts, colData= design, design= ~ Treatment)
ddsB$Treatment <- relevel(ddsB$Treatment, ref="u") ## Set the reference condition as the primary tumor.
ddsB <- DESeq(ddsB)
resB <- results(ddsB)

#ss <- fitModel()
tre_pvals <- cbind(refGene, FDR_TAM= res$padj, FC_TAM= res$log2FoldChange,
                             FDR_BBCA= resB$padj, FC_BBCA= resB$log2FoldChange)

## Number of genes w/ p< 0.01
NROW(tre_pvals[tre_pvals$FDR_TAM < 0.01,])
NROW(tre_pvals[tre_pvals$FDR_BBCA < 0.01,])

write.table(tre_pvals[order(tre_pvals$FDR_TAM),], "TRE.Signif.Changes.TamRes.tsv", row.names=FALSE, col.names=TRUE, quote=FALSE, sep="\t")

