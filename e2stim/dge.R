## get correlations between MCF-7 samples

## Read in refseq genes.
refGene <- read.table("../annotations/refGene.bed.gz")
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

GC <- function(prefix, path="../data/") {
 pl <- load.bigWig(paste(path, prefix, "_plus.bw", sep=""))
 mn <- load.bigWig(paste(path, prefix, "_minus.bw", sep=""))

 ## Count reads in each ...
 cnts <- bed6.region.bpQuery.bigWig(pl, mn, bodies, abs.value = TRUE)
 cnts
}

gene_body_counts <- cbind(B7=GC("B7_CTRL"), 
                          B7.e2=GC("B7_E2"),
                          B7.tam=GC("B7_TAM"),
                          B7.Be2=GC("B7_E2_BBCA"),
                          B7.Btam=GC("B7_TAM_BBCA"),

                          G11=GC("G11_CTRL"),
                          G11.e2=GC("G11_E2"),
                          G11.tam=GC("G11_TAM"),
                          G11.Be2=GC("G11_E2_BBCA"),
                          G11.Btam=GC("G11_TAM_BBCA"))

## 
require(edgeR)
PVAL <- 0.05

## Build experimental design matrix
sampleID <- colnames(gene_body_counts)
trt  <- rep(c("u", "ze", "zt", "ze", "zt"),2)
bbca <- rep(c("u", "u", "u", "zb", "zb"), 2)
data.frame(sampleID, trt, bbca)

#fitModel <- function() {

  ####################
  ## Fit LRT model for BBCA.
  design <- model.matrix(~trt+bbca)
  rownames(design) <- colnames(gene_body_counts)

  ## Estimate dispersions.
  dge <- estimateGLMCommonDisp(gene_body_counts, design, verbose=TRUE)
  #dge <- estimateGLMTrendedDisp(dge,design)

  ## Fit neg. binom. GLM to sensitive/ resistent.  Compute p-values using LRT.
  fit <- glmFit(gene_body_counts,design,dge)

  cont <- makeContrasts(bbcazb, levels=design)
  bbca_ss <- glmLRT(fit,contrast=cont)

  ###################
  ## Fit E2 and Tam stim. separately in B7 and G11.
  design <- model.matrix(~rep(c("u", "ze", "u", "b", "b"),2)+c(rep("b7",5), rep("g11",5)))
  rownames(design) <- colnames(gene_body_counts)
  dge <- estimateGLMCommonDisp(gene_body_counts, design, verbose=TRUE)

  y <- DGEList(counts= gene_body_counts[,c(1,2)], group=1:2)
  b7e2 <- exactTest(y, dispersion= dge)
  topTags(b7e2)

  y <- DGEList(counts= gene_body_counts[,c(6,7)], group=1:2)
  g11e2 <- exactTest(y, dispersion= dge)
  topTags(g11e2)

##
gene_pvals <- cbind(refGene, FDR_bbca= p.adjust(bbca_ss$table$PValue), FC_BBCA= bbca_ss$table$logFC, 
				FDR_b7_e2= p.adjust(b7e2$table$PValue), FC_b7_e2= b7e2$table$logFC,
				FDR_g11_e2= p.adjust(g11e2$table$PValue), FC_g11_e2= g11e2$table$logFC)
## Questions...
## (1) Is E2 stimulation functional in resistant lines.
## YES!
sum(gene_pvals$FDR_g11_e2 < PVAL) ## Plus, the genes check out.
sum(gene_pvals$FDR_b7_e2 < PVAL) 


## (2) Does the E2 target gene set change.
sig_b7 <- (gene_pvals$FDR_b7_e2 < PVAL)
sig_g11<- (gene_pvals$FDR_g11_e2 < PVAL)

fc <- gene_pvals[,c("FC_b7_e2", "FC_g11_e2")]

plot(fc[sig_b7,]); abline(h=0); abline(v=0) ## Generally consistent direction, smaller effect sizes.
plot(fc[sig_g11,]); abline(h=0); abline(v=0)

sum(sig_b7&sig_g11)/sum(sig_g11) ## 60% are significant by both.


