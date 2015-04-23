## get correlations between MCF-7 samples

## Read in refseq genes.
tres <- read.table("mcf7_tamSR.bedgraph.bed.gz")

## Read in bigWigs.
require(bigWig)

path="/local/storage/projects/RawSequenceFiles/mcf7tamRes/"

B7_pl <- load.bigWig(paste(path, "5089_5598_18732_H73KFGXX_MCF7_B7_TTAGGC_R1_plus.bw",sep=""))
B7_mn <- load.bigWig(paste(path, "5089_5598_18732_H73KFGXX_MCF7_B7_TTAGGC_R1_minus.bw",sep=""))
G11_pl <- load.bigWig(paste(path, "5089_5598_18733_H73KFGXX_MCF7_G11_TGACCA_R1_plus.bw", sep=""))
G11_mn <- load.bigWig(paste(path, "5089_5598_18733_H73KFGXX_MCF7_G11_TGACCA_R1_minus.bw", sep=""))
H9_pl <- load.bigWig(paste(path, "5089_5598_18735_H73KFGXX_MCF7_H9_GCCAAT_R1_plus.bw",sep=""))
H9_mn <- load.bigWig(paste(path, "5089_5598_18735_H73KFGXX_MCF7_H9_GCCAAT_R1_minus.bw",sep=""))
C11_pl <- load.bigWig(paste(path, "5089_5598_18733_H73KFGXX_MCF7_G11_TGACCA_R1_plus.bw",sep=""))
C11_mn <- load.bigWig(paste(path, "5089_5598_18734_H73KFGXX_MCF7_C11_ACAGTG_R1_minus.bw",sep=""))

B7_BBCA_pl <- load.bigWig(paste(path, "5089_5598_18730_H73KFGXX_MCF7_B7_BBCA_ATCACG_R1_plus.bw",sep=""))
B7_BBCA_mn <- load.bigWig(paste(path, "5089_5598_18730_H73KFGXX_MCF7_B7_BBCA_ATCACG_R1_minus.bw",sep=""))
G11_BBCA_pl <- load.bigWig(paste(path, "5089_5598_18731_H73KFGXX_MCF7_G11_BBCA_CGATGT_R1_plus.bw", sep=""))
G11_BBCA_mn <- load.bigWig(paste(path, "5089_5598_18731_H73KFGXX_MCF7_G11_BBCA_CGATGT_R1_minus.bw", sep=""))

countTRE <- function(plus, minus, tres) {
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
require(edgeR)
PVAL <- 0.01

## Build experimental design matrix
sampleID <- colnames(tre_counts)
tamres   <- c("sen", "sen", "res", "res", "sen", "res")
bbcatrt  <- c("u", "u", "u", "u", "t", "t")
data.frame(sampleID, tamres, bbcatrt)

#fitModel <- function() {
  design <- model.matrix(~tamres+bbcatrt)
  rownames(design) <- colnames(tre_counts)

  ## Estimate dispersions.
  dge <- estimateGLMCommonDisp(tre_counts, design, verbose=TRUE)
  #dge <- estimateGLMTrendedDisp(dge,design)

  ## Fit neg. binom. GLM to sensitive/ resistent.  Compute p-values using LRT.
  fit <- glmFit(tre_counts,design,dge)

  sen_cont <- makeContrasts(tamressen, levels=design)
  ss <- glmLRT(fit,contrast=sen_cont)

  ## Fit bbca.
  bbca_cont <- makeContrasts(bbcatrtu, levels=design)
  sb <- glmLRT(fit,contrast=bbca_cont)

#  return(ss, )
#}

#ss <- fitModel()
tre_pvals <- cbind(tres, FDR_TAM= p.adjust(ss$table$PValue), FC_TAM= ss$table$logFC, 
								FDR_BBCA= p.adjust(sb$table$PValue), FC_BBCA= sb$table$logFC) 
								#FDR_BBCA_G11= p.adjust(sbg11$table$PValue), FC_BBCA_G11= sbg11$table$logFC)

head(tre_pvals[order(tre_pvals$FDR_TAM),c(1:3,6:9)], n=30)
head(tre_pvals[order(tre_pvals$FDR_BBCA),c(1:3,6:7,10:11)], n=30)
head(tre_pvals[order(tre_pvals$FDR_BBCA_G11),c(1:3,6:7,12:13)], n=30)

## Number of genes w/ p< 0.01
NROW(tre_pvals[tre_pvals$FDR_TAM < 0.01,])
NROW(tre_pvals[tre_pvals$FDR_BBCA < 0.01,])

write.table(tre_pvals[order(tre_pvals$FDR_TAM),], "TRE.Signif.Changes.TamRes.tsv", row.names=FALSE, col.names=TRUE, quote=FALSE, sep="\t")

