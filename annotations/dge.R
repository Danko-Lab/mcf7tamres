## get correlations between MCF-7 samples

## Read in refseq genes.
refGene <- read.table("refGene.bed.gz")
refGene <- refGene[grep("random|Un|hap", refGene$V1, invert=TRUE),]
refGene <- refGene[(refGene$V3-refGene$V2)>5000,]

bodies <- refGene
bodies$V2[bodies$V6 == "+"] <-bodies$V2[bodies$V6 == "+"]+1000
bodies$V3[bodies$V6 == "-"] <- bodies$V3[bodies$V6 == "-"]-1000

tss <- refGene
tss$V3[tss$V6 == "+"] <-tss$V2[tss$V6 == "+"]+1000
tss$V2[tss$V6 == "-"] <- tss$V3[tss$V6 == "-"]-1000

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

## Count reads in each ...
B7 <- bed6.region.bpQuery.bigWig(B7_pl, B7_mn, bodies, abs.value = TRUE)
G11<- bed6.region.bpQuery.bigWig(G11_pl, G11_mn, bodies, abs.value = TRUE)
H9 <- bed6.region.bpQuery.bigWig(H9_pl, H9_mn, bodies, abs.value = TRUE)
C11<- bed6.region.bpQuery.bigWig(C11_pl, C11_mn, bodies, abs.value = TRUE)
B7b <- bed6.region.bpQuery.bigWig(B7_BBCA_pl, B7_BBCA_mn, bodies, abs.value = TRUE)
G11b<- bed6.region.bpQuery.bigWig(G11_BBCA_pl, G11_BBCA_mn, bodies, abs.value = TRUE)
gene_body_counts <- cbind(B7, C11, H9, G11, B7b, G11b)

## 
require(edgeR)
PVAL <- 0.01

## Build experimental design matrix
sampleID <- colnames(gene_body_counts)
tamres   <- c("sen", "sen", "res", "res", "sen", "res")
bbcatrt  <- c("u", "u", "u", "u", "t", "t")
data.frame(sampleID, tamres, bbcatrt)

#fitModel <- function() {
  design <- model.matrix(~tamres+bbcatrt)
  rownames(design) <- colnames(gene_body_counts)

  ## Estimate dispersions.
  dge <- estimateGLMCommonDisp(gene_body_counts, design, verbose=TRUE)
  #dge <- estimateGLMTrendedDisp(dge,design)

  ## Fit neg. binom. GLM to sensitive/ resistent.  Compute p-values using LRT.
  fit <- glmFit(gene_body_counts,design,dge)

  sen_cont <- makeContrasts(tamressen, levels=design)
  ss <- glmLRT(fit,contrast=sen_cont)

  ## Fit bbca.
  bbca_cont <- makeContrasts(bbcatrtu, levels=design)
  sb <- glmLRT(fit,contrast=bbca_cont)

#  return(ss, )
#}

#ss <- fitModel()
gene_pvals <- cbind(refGene, FDR_TAM= p.adjust(ss$table$PValue), FC_TAM= ss$table$logFC, FDR_BBCA= p.adjust(sb$table$PValue), FC_BBCA= sb$table$logFC)

head(gene_pvals[order(gene_pvals$FDR_TAM),c(1:3,6:9)], n=30)
head(gene_pvals[order(gene_pvals$FDR_BBCA),c(1:3,6:7,10:11)], n=30)

## Number of genes w/ p< 0.01
NROW(unique(gene_pvals$V7[gene_pvals$FDR_TAM < 0.01]))
NROW(unique(gene_pvals$V7[gene_pvals$FDR_BBCA < 0.01]))

write.table(gene_pvals[order(gene_pvals$FDR_TAM),c(1:3,6:11)], "Signif.Changes.TamRes.tsv", row.names=FALSE, col.names=TRUE, quote=FALSE, sep="\t")


