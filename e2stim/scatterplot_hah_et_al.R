## Creates a scatterplot comparing fold-changes in Hah et. al. 2011 to 
## new data.

# read data.
hah <- read.table("GSE27463_RefSeq.reg.tsv.gz", header=TRUE)
hg19 <- read.table("../annotations/refGene.bed.gz")

# Come up with scheme for indexing Hah et. al. in hg19
hah <- hah[!is.na(match(hah[,5], hg19[,4])),]
indx <- match(hah[,5], hg19[,4])
hah[,1:4] <- hg19[indx,c(1:3,6)]


# Count reads in our new data.
require(bigWig)
data_path<- "/local/storage/projects/mcf7tamres/data"

getCounts <- function(prefix) {
  B7c_pl <- load.bigWig(paste(data_path,"/",prefix,"_CTRL_plus.bw", sep=""))
  B7c_mn <- load.bigWig(paste(data_path,"/",prefix,"_CTRL_minus.bw", sep=""))
  B7e_pl <- load.bigWig(paste(data_path,"/",prefix,"_E2_plus.bw", sep=""))
  B7e_mn <- load.bigWig(paste(data_path,"/",prefix,"_E2_minus.bw", sep=""))
  B7t_pl <- load.bigWig(paste(data_path,"/",prefix,"_TAM_plus.bw", sep=""))
  B7t_mn <- load.bigWig(paste(data_path,"/",prefix,"_TAM_minus.bw", sep=""))
  B7eb_pl <- load.bigWig(paste(data_path,"/",prefix,"_E2_BBCA_plus.bw", sep=""))
  B7eb_mn <- load.bigWig(paste(data_path,"/",prefix,"_E2_BBCA_minus.bw", sep=""))

  B7c_counts <- (B7c_pl$basesCovered*B7c_pl$mean)+abs(B7c_mn$basesCovered*B7c_mn$mean)
  B7e_counts <- (B7e_pl$basesCovered*B7e_pl$mean)+abs(B7e_mn$basesCovered*B7e_mn$mean)
  B7t_counts <- (B7t_pl$basesCovered*B7t_pl$mean)+abs(B7t_mn$basesCovered*B7t_mn$mean)
  B7eb_counts <- (B7eb_pl$basesCovered*B7eb_pl$mean)+abs(B7eb_mn$basesCovered*B7eb_mn$mean)

  B7c <- bed6.region.bpQuery.bigWig(B7c_pl, B7c_mn, hah[,c(1:3,5:6,4)], abs.value=TRUE)
  B7e2<- bed6.region.bpQuery.bigWig(B7e_pl, B7e_mn, hah[,c(1:3,5:6,4)], abs.value=TRUE)
  B7t <- bed6.region.bpQuery.bigWig(B7t_pl, B7t_mn, hah[,c(1:3,5:6,4)], abs.value=TRUE)
  B7eb<- bed6.region.bpQuery.bigWig(B7eb_pl, B7eb_mn, hah[,c(1:3,5:6,4)], abs.value=TRUE)

  rL <- list()
  rL$Ce  <- log((B7e2+1)/B7e_counts/((B7c+1)/B7c_counts), 2)
  rL$Cet <- log((B7t+1)/B7t_counts/((B7c+1)/B7c_counts), 2)
  rL$Ceb <- log((B7eb+1)/B7eb_counts/((B7c+1)/B7c_counts), 2)

  rL$all_hah <- log((hah$E2_40m+1)/ (hah$VEH+1), 2)
  rL
}

rL <- getCounts("B7")

use <- hah$E2.40m.qVal < 0.01 #& abs(b7) > 1
up  <- all_hah > 0
down<- all_hah < 0
plot(rL$Ce[use], rL$all_hah[use]);abline(h=0);abline(v=0)
plot(rL$Ce[use], rL$Ceb[use]);abline(h=0);abline(v=0)
plot(rL$Ce[use], rL$Cet[use]);abline(h=0);abline(v=0)

source("~/NHP/lib/densScatterplot.R")
densScatterplot(rL$Ce[use], rL$all_hah[use]);abline(h=0);abline(v=0)
densScatterplot(rL$Ceb[use], rL$all_hah[use]);abline(h=0);abline(v=0)

#########################
## Add violinplots.
use <- hah$E2.40m.qVal < 0.01 & abs(b7) > 1 ## Looking at the effects of  BBCA, require changed in this experiment.

gbc <- read.csv("../annotations/gene_body_counts.csv")

pdf("BBCAonE2.violinplot.pdf")
require(vioplot)
vioplot(b7[use&up], b7t[use&up], b7b[use&up], b7[use&down], b7t[use&down], b7b[use&down], names=c("Up E2", "Up TAM", "Up E2+BBCA", "DN E2", "DN Tam", "DN E2+BBCA")); abline(h=0)
dev.off()

##############################
## Compare to G11
rLg11 <- getCounts("G11")

plot(rL$Ce[use], rLg11$Ce[use]);abline(h=0);abline(v=0); abline(0,1)

densScatterplot(b7[use], g11[use]);abline(h=0);abline(v=0)

## Compare to BBCA
B7eb_pl <- load.bigWig(paste(data_path,"/B7_E2_plus.bw", sep=""))
B7eb_mn <- load.bigWig(paste(data_path,"/B7_E2_minus.bw", sep=""))



