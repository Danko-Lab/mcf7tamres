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
rLg11 <- getCounts("G11")

use <- hah$E2.40m.qVal < 0.05 #& abs(b7) > 1
up  <- rL$all_hah > 0
down<- rL$all_hah < 0
cor.test(rL$Ce[use], rL$all_hah[use], method="spearman")
plot(rL$Ce[use], rL$all_hah[use]);abline(h=0);abline(v=0)
plot(rL$Ce[use], rL$Ceb[use]);abline(h=0);abline(v=0)
plot(rL$Ce[use], rL$Cet[use]);abline(h=0);abline(v=0)

source("~/NHP/lib/densScatterplot.R")
densScatterplot(rL$Ce[use], rL$all_hah[use]);abline(h=0);abline(v=0)
densScatterplot(rL$Ceb[use], rL$all_hah[use]);abline(h=0);abline(v=0)

#########################
## Write out conditioned B7 dataset.
write <- use #& abs(rL$Ce) > 0 & (rL$Ce*rL$all_hah)>0 #& abs(rL$Ce) > 0.75 
plot(rL$Ce[write], rL$all_hah[write])
cor.test(rL$Ce[write], rL$all_hah[write])
write.table(data.frame(hah[write,1:6], fold.change.B7= rL$Ce[write], fold.change.hah= log(hah$E2_40m/ hah$VEH, 2)[write]), "Hah.B7.changes.tsv", row.names=FALSE, quote=FALSE, sep="\t")


#########################
## Add violinplots.
#use <- hah$E2.40m.qVal < 0.01 & abs(rL$Ce) > 0.25 & abs(rLg11$Ce) > 0.25 ## Looking at the effects of  BBCA, require changed >2-fold in E2 B7.  
use <- hah$E2.40m.qVal < 0.01 & abs(rL$Ce) > 1 #& abs(rLg11$Ce) > 1 ## Looking at the effects of  BBCA, require changed >2-fold in E2 B7.  

gbc <- read.csv("../annotations/gene_body_counts.csv")

pdf("BBCAonE2.violinplot.pdf")
require(vioplot)
vioplot(rL$Ce[use&up], rL$Cet[use&up], rL$Ceb[use&up], rL$Ce[use&down], rL$Cet[use&down], rL$Ceb[use&down], names=c("Up E2", "Up TAM", "Up E2+BBCA", "DN E2", "DN Tam", "DN E2+BBCA")); abline(h=0)
dev.off()

##############################
## Compare to G11

plot(rL$Ce[use], rLg11$Ce[use]);abline(h=0);abline(v=0); abline(0,1)

pdf("B7vG11.E2Response.pdf")
 densScatterplot(rL$Ce[use], rLg11$Ce[use]);abline(h=0);abline(v=0)
 vioplot(rLg11$Ce[use&up], rLg11$Cet[use&up], rLg11$Ceb[use&up], rLg11$Ce[use&down], rLg11$Cet[use&down], rLg11$Ceb[use&down], names=c("Up E2", "Up TAM", "Up E2+BBCA", "DN E2", "DN Tam", "DN E2+BBCA")); abline(h=0)
 vioplot(rL$Cet[use&up], rLg11$Cet[use&up], rL$Cet[use&down], rLg11$Cet[use&down], names=c("Up B7 t", "Up G11 t", "Dn B7 t", "Dn G11 t")); abline(h=0)
 vioplot(rL$Ce[use&up], rLg11$Ce[use&up], rL$Ce[use&down], rLg11$Ce[use&down], names=c("Up B7 e", "Up G11 e", "Dn B7 e", "Dn G11 e")); abline(h=0)
dev.off()

##  

2^summary(rLg11$Ce[use&up])
2^summary(rL$Ce[use&up])

wilcox.test(rLg11$Ce[use&up], rL$Ce[use&up])

