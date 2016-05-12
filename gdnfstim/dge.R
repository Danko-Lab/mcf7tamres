## get correlations between MCF-7 samples

## Read in refseq genes.
#refGene <- read.table("refGene.bed.gz")
refGene <- read.table("../annotations/tuSelecter/final_tus.txt", header=TRUE)
#refGene <- rbind(refGene, read.table("../annotations/tuSelecter/final_tus.ESR1_GREB1.txt", header=TRUE))

refGene <- refGene[grep("random|Un|hap", refGene$TXCHROM, invert=TRUE),]
refGene <- refGene[(refGene$TXEND-refGene$TXSTART)>4000,]

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
## Limit analysis to resistent lines to compare FC with sensitive.
  indexres <- c(7:12)
  gdnftrt_r <- gdnftrt[indexres]
  design <- model.matrix(~gdnftrt_r)
  rownames(design) <- colnames(gene_body_counts[,indexres])

  ## Estimate dispersions.
  dge <- estimateGLMCommonDisp(gene_body_counts[,indexres], design, verbose=TRUE)
  fit <- glmFit(gene_body_counts[,indexres],design,dge)

  ## Fit 1h GDNF.
  gdnf_1h_cont <- makeContrasts(gdnftrt_rgdnf1h, levels=design)
  gdnf_1h_lrt_r <- glmLRT(fit,contrast=gdnf_1h_cont)

  ## Fit 24h GDNF.
  gdnf_24h_cont <- makeContrasts(gdnftrt_rgdnf24h, levels=design)
  gdnf_24h_lrt_r <- glmLRT(fit,contrast=gdnf_24h_cont)

#ss <- fitModel()
gene_pvals <- cbind(refGene, FDR_TAM= p.adjust(ss$table$PValue), FC_TAM= ss$table$logFC,
                                                                FDR_GDNF_1h= p.adjust(gdnf_1h_lrt$table$PValue), FC_GDNF_1h= gdnf_1h_lrt$table$logFC,
                                                                FDR_GDNF_24h= p.adjust(gdnf_24h_lrt$table$PValue), FC_GDNF_24h= gdnf_24h_lrt$table$logFC,
                                                                FDR_GDNF_1h_S= p.adjust(gdnf_1h_lrt_s$table$PValue), FC_GDNF_1h_S= gdnf_1h_lrt_s$table$logFC,
                                                                FDR_GDNF_24h_S= p.adjust(gdnf_24h_lrt_s$table$PValue), FC_GDNF_24h_S= gdnf_24h_lrt_s$table$logFC,
                                                                FDR_GDNF_1h_R= p.adjust(gdnf_1h_lrt_r$table$PValue), FC_GDNF_1h_R= gdnf_1h_lrt_r$table$logFC,
                                                                FDR_GDNF_24h_R= p.adjust(gdnf_24h_lrt_r$table$PValue), FC_GDNF_24h_R= gdnf_24h_lrt_r$table$logFC
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

addlab <- function(ss, gene_ID, ...) {
  ig <- which(gene_pvals$GENENAME == gene_ID)
  io <- ig[which.min(gene_pvals$FDR_TAM[ig])]
  text(ss$table$logCPM[io], ss$table$logFC[io], labels= gene_ID, cex= 0.7, pos= 3, ...)
} 

laball_1h <- function(ss) {
addlab(ss, "EGR2")
addlab(ss, "EGR3")
addlab(ss, "FOSL1")
addlab(ss, "NR4A3")
addlab(ss, "LRRC15")
addlab(ss, "PMF1")
addlab(ss, "ETS2")
}

laball_24h <- function(ss) {
addlab(ss, "SLC45A1")
addlab(ss, "RASGEF1A")
addlab(ss, "LONRF3")
addlab(ss, "PLD1")
addlab(ss, "ERRFI1")
addlab(ss, "SHC4")
addlab(ss, "VAV3")
addlab(ss, "NRP1")

addlab(ss, "ESR1")
addlab(ss, "PGR")

addlab(ss, "GATA4")
addlab(ss, "KIF12")
addlab(ss, "CTB-164N12.1")
addlab(ss, "RET")
addlab(ss, "ACACB")
addlab(ss, "GFRA1")
}

WriteMAPlot <- function(ss, laball, prefix) {

## MA-plot
pdf(paste("MAPlot.",prefix,".pdf", sep=""))

indx <- rep(TRUE, NROW(ss)) #ss$table$logCPM>0
sign <- (p.adjust(ss$table$PValue) < 0.01)
maPlot(logAbundance= ss$table$logCPM[indx], logFC= ss$table$logFC[indx], de.tags= sign[indx], pch=19)

maPlot(logAbundance= ss$table$logCPM[indx], logFC= ss$table$logFC[indx], de.tags= sign[indx], pch=19)
laball(ss)

indx <- which(p.adjust(ss$table$PValue) < 1e-5 | abs(ss$table$logFC) > 1)
maPlot(logAbundance= ss$table$logCPM[indx], logFC= ss$table$logFC[indx], pch=19)
laball(ss)

dev.off()

}

WriteMAPlot(ss= gdnf_1h_lrt_s, laball_1h, "1h_s")
WriteMAPlot(ss= gdnf_24h_lrt_s, laball_24h, "24h_s")

WriteMAPlot(ss= gdnf_1h_lrt, laball_1h, "1h_all")
WriteMAPlot(ss= gdnf_24h_lrt, laball_24h, "24h_all")


