## get correlations between MCF-7 samples

## Read in refseq genes.
refGene <- read.table("../dreg/mcf7.gdnf.dREG_HD.bed")
refGene <- cbind(refGene, name=paste(refGene[,1], ":",refGene[,2], "-", refGene[,3], sep=""), refGene) 

refGene[,2] <- refGene[,2]-200
refGene[,3] <- refGene[,3]+200

## Read in bigWigs.
require(bigWig)

countBigWig <- function(prefix, path="../data/") {
 pl <- load.bigWig(paste(path, "MCF-7_", prefix, "_plus.bw", sep=""))
 mn <- load.bigWig(paste(path, "MCF-7_", prefix, "_minus.bw", sep=""))

 counts <- bed.region.bpQuery.bigWig(pl, refGene[1:3], abs.value = TRUE) + abs(bed.region.bpQuery.bigWig(mn, refGene[,1:3], abs.value = TRUE))
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

## Count reads in each ...

require(edgeR)
PVAL <- 0.01

## Build experimental design matrix
sampleID <- colnames(gene_body_counts)
tamres   <- c(rep("sen", 6), rep("res", 6))
gdnftrt  <- rep(c("gdnf0h", "gdnf1h", "gdnf24h"), 4)
data.frame(sampleID, tamres, gdnftrt)

  design <- model.matrix(~tamres+gdnftrt)
  rownames(design) <- colnames(gene_body_counts)

  ## Estimate dispersions.
  dge <- estimateGLMCommonDisp(gene_body_counts, design, verbose=TRUE)
  #dge <- estimateGLMTrendedDisp(dge,design)

  ## Fit neg. binom. GLM to sensitive/ resistent.  Compute p-values using LRT.
  fit <- glmFit(gene_body_counts,design,dge)

  sen_cont <- makeContrasts(tamressen, levels=design)
  ss <- glmLRT(fit,contrast=sen_cont)

  ## Fit 1h GDNF.
  gdnf_1h_cont <- makeContrasts(gdnftrtgdnf1h, levels=design)
  gdnf_1h_lrt <- glmLRT(fit,contrast=gdnf_1h_cont)

  ## Fit 24h GDNF.
  gdnf_24h_cont <- makeContrasts(gdnftrtgdnf24h, levels=design)
  gdnf_24h_lrt <- glmLRT(fit,contrast=gdnf_24h_cont)

############################################################################################33
## Limit analysis to sensitive lines where the effect sizes are expected to be largest.
  gdnftrt_s <- gdnftrt[1:6]
  design <- model.matrix(~gdnftrt_s)
  rownames(design) <- colnames(gene_body_counts[,1:6])

  ## Estimate dispersions.
  dge <- estimateGLMCommonDisp(gene_body_counts[,1:6], design, verbose=TRUE)
  fit <- glmFit(gene_body_counts[,1:6],design,dge)

  ## Fit 1h GDNF.
  gdnf_1h_cont <- makeContrasts(gdnftrt_sgdnf1h, levels=design)
  gdnf_1h_lrt_s <- glmLRT(fit,contrast=gdnf_1h_cont)

  ## Fit 24h GDNF.
  gdnf_24h_cont <- makeContrasts(gdnftrt_sgdnf24h, levels=design)
  gdnf_24h_lrt_s <- glmLRT(fit,contrast=gdnf_24h_cont)

############################################################################################33
## Limit analysis to sensitive lines, 1h ... does the strange bias for down-regulation hold?!
  indx <- c(1:2,4:5)
  gdnftrt_s <- gdnftrt[indx]
  design <- model.matrix(~gdnftrt_s)
  rownames(design) <- colnames(gene_body_counts[,indx])

  ## Estimate dispersions.
  dge <- estimateGLMCommonDisp(gene_body_counts[,indx], design, verbose=TRUE)
  fit <- glmFit(gene_body_counts[,indx],design,dge)

  ## Fit 1h GDNF.
  gdnf_1h_cont <- makeContrasts(gdnftrt_sgdnf1h, levels=design)
  gdnf_1h_lrt_s <- glmLRT(fit,contrast=gdnf_1h_cont)

  ## And now 24h...
  indx <- c(1,3:4,6)
  gdnftrt_s <- gdnftrt[indx]
  design <- model.matrix(~gdnftrt_s)
  rownames(design) <- colnames(gene_body_counts[,indx])

  ## Estimate dispersions.
  dge <- estimateGLMCommonDisp(gene_body_counts[,indx], design, verbose=TRUE)
  fit <- glmFit(gene_body_counts[,indx],design,dge)

  ## Fit 24h GDNF.
  gdnf_24h_cont <- makeContrasts(gdnftrt_sgdnf24h, levels=design)
  gdnf_24h_lrt_s <- glmLRT(fit,contrast=gdnf_24h_cont)

#ss <- fitModel()
gene_pvals <- cbind(refGene, FDR_TAM= p.adjust(ss$table$PValue), FC_TAM= ss$table$logFC,
                                                                FDR_GDNF_1h= p.adjust(gdnf_1h_lrt$table$PValue), FC_GDNF_1h= gdnf_1h_lrt$table$logFC,
                                                                FDR_GDNF_24h= p.adjust(gdnf_24h_lrt$table$PValue), FC_GDNF_24h= gdnf_24h_lrt$table$logFC,
                                                                FDR_GDNF_1h_S= p.adjust(gdnf_1h_lrt_s$table$PValue), FC_GDNF_1h_S= gdnf_1h_lrt_s$table$logFC,
                                                                FDR_GDNF_24h_S= p.adjust(gdnf_24h_lrt_s$table$PValue), FC_GDNF_24h_S= gdnf_24h_lrt_s$table$logFC)

sum(gene_pvals$FDR_GDNF_24h_S < 0.05 & gene_pvals$FC_GDNF_24h_S > 0)
sum(gene_pvals$FDR_GDNF_1h_S < 0.05 & gene_pvals$FC_GDNF_1h_S > 0)

sum(gene_pvals$FDR_GDNF_24h_S < 0.05 & gene_pvals$FC_GDNF_24h_S < 0)
sum(gene_pvals$FDR_GDNF_1h_S < 0.05 & gene_pvals$FC_GDNF_1h_S < 0)



write.table(gene_pvals, "GDNF_treatment.dREG-HD.tsv", sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)

head(gene_pvals[order(gene_pvals$FDR_TAM),], 20)
head(gene_pvals[order(gene_pvals$FDR_GDNF_24h),c(7,10:NCOL(gene_pvals))], 100)
head(gene_pvals[order(gene_pvals$FDR_GDNF_1h),c(7,10:NCOL(gene_pvals))], 100)

head(gene_pvals[gene_pvals$V7 == "ESR1",c(7,10:NCOL(gene_pvals))], 100)
head(gene_pvals[gene_pvals$V7 == "CCND1",c(7,10:NCOL(gene_pvals))], 100)
head(gene_pvals[gene_pvals$V7 == "GREB1",c(7,10:NCOL(gene_pvals))], 100)
head(gene_pvals[gene_pvals$V7 == "HOXC13",c(7,10:NCOL(gene_pvals))], 100)


head(gene_pvals[gene_pvals$V7 == "AKT",c(7,10:NCOL(gene_pvals))], 100)


## MA-plot
pdf("MAPlot.pdf")
sign <- which(ss$table$PValue < 0.01)
#samp <- sample(NROW(gene_pvals), 0.3*NROW(gene_pvals))
#maPlot(logAbundance= ss$table$logCPM[c(sign, samp)], logFC= ss$table$logFC[c(sign, samp)], de.tags= 1:NROW(sign), pch=19)
maPlot(logAbundance= ss$table$logCPM, logFC= ss$table$logFC, de.tags= sign, pch=19)

addlab <- function(gene_ID, ...) {
  ig <- which(gene_pvals$V7 == gene_ID)
  io <- ig[which.min(gene_pvals$FDR_TAM[ig])]
  text(ss$table$logCPM[io], ss$table$logFC[io], labels= gene_ID, cex= 0.7, pos= 3, ...)
} 

laball <- function() {
addlab("GREB1")
addlab("GDNF")
addlab("PGR")
addlab("BMP5")

addlab("CD36")
addlab("MPPED1")

#addlab("ESR1", col="white")
#addlab("ERBB2", col="white")
#addlab("REL", col="white")
#addlab("PIK3CA", col="white")
}
laball()

indx <- which(gene_pvals$FDR_TAM < 1e-10 | ss$table$logCPM > 10)
maPlot(logAbundance= ss$table$logCPM[indx], logFC= ss$table$logFC[indx], pch=19)
laball()
dev.off()

head(gene_pvals[order(gene_pvals$FDR_TAM),c(1:3,6:9)], n=30)
head(gene_pvals[order(gene_pvals$FDR_BBCA),c(1:3,6:7,10:11)], n=30)
head(gene_pvals[order(gene_pvals$FDR_BBCA_G11),c(1:3,6:7,12:13)], n=30)

## Number of genes w/ p< 0.01
NROW(unique(gene_pvals$V7[gene_pvals$FDR_TAM < 0.01]))
NROW(unique(gene_pvals$V7[gene_pvals$FDR_BBCA < 0.01]))

write.table(gene_pvals[order(gene_pvals$FDR_TAM),c(1:3,6:11)], "Signif.Changes.TamRes.tsv", row.names=FALSE, col.names=TRUE, quote=FALSE, sep="\t")


## SANITY CHECKS.
sanityCheck <- function(gene) {
	indx <- grep(gene, refGene$V7)
	print(data.frame(MGID= as.character(refGene[indx,7]), gene_body_counts[indx,]))
}
sanityCheck("CD36")
sanityCheck("MPPED")
sanityCheck("ESR1")



