## get correlations between MCF-7 samples

## Read in refseq genes.
dregHD<-read.table("mcf7_gdnf.bedgraph.gz.bed.gz_dREG_HD_relaxed_fixed90.bed"); dregHD <- dregHD[dregHD$V3-dregHD$V2 > 0,]

tres <- dregHD

## Read in bigWigs.
require(bigWig)

countBigWig <- function(prefix, path="../data/", rpkm=FALSE) {
 pl <- load.bigWig(paste(path, "MCF-7_", prefix, "_plus.bw", sep=""))
 mn <- load.bigWig(paste(path, "MCF-7_", prefix, "_minus.bw", sep=""))

 counts <- bed.region.bpQuery.bigWig(pl, tres, abs.value = TRUE)+bed.region.bpQuery.bigWig(mn, tres, abs.value = TRUE)
        if(rpkm==TRUE) {
                counts <- counts * (1000/(tres[,3]-tres[,2])) * (1e6/(abs(pl$mean)*pl$basesCovered+abs(mn$mean)*mn$basesCovered))
        }

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

tre_counts <- cbind( countTimePoints("B7"), countTimePoints("C11") )#, countTimePoints("G11"), countTimePoints("H9") ) ## Focus only on sensitive lines.
save.image("gdnf.tres.Counts.RData")

## 
require(DESeq2)
PVAL <- 0.01

## Build experimental design matrix
sampleID <- colnames(tre_counts)
tamres   <- c(rep("sen", 6))#, rep("res", 6))
gdnftrt  <- rep(c("gdnf0h", "gdnf1h", "gdnf24h"), 2)#4)
data.frame(sampleID, tamres, gdnftrt)

## Create DESeq2 object.  Note that Condition (tamres) is excluded from the current design focusing on the TamS lines.
design <- data.frame(Treatment= gdnftrt, row.names=colnames(tre_counts))
dds <- DESeqDataSetFromMatrix(countData= tre_counts, colData= design, design= ~ Treatment)
dds$Treatment <- relevel(dds$Treatment, ref="gdnf0h") ## Set the reference condition as the primary tumor.
dds <- DESeq(dds)

gdnf_1h <- results(dds, contrast= c("Treatment", "gdnf1h", "gdnf0h"))
gdnf_24h <- results(dds, contrast= c("Treatment", "gdnf24h", "gdnf0h"))

sum(gdnf_1h$padj < 0.01, na.rm=TRUE)
sum(gdnf_24h$padj < 0.01, na.rm=TRUE)

## Summarize data.
tre_pvals <- cbind(tres,  FDR_GDNF_1h= gdnf_1h$padj, FC_GDNF_1h= gdnf_1h$log2FoldChange,
                          FDR_GDNF_24h= gdnf_24h$padj, FC_GDNF_24h= gdnf_24h$log2FoldChange)

## Number of genes w/ p< 0.01
NROW(tre_pvals[tre_pvals$FDR_GDNF_1h < 0.01,])
NROW(tre_pvals[tre_pvals$FDR_GDNF_24h < 0.01,])

write.table(tre_pvals, "GDNF_treatment.dREG-HD.tsv", row.names=FALSE, col.names=TRUE, quote=FALSE, sep="\t")

