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
data_path<- "/local/storage/projects/RawSequenceFiles/mcf7tamRes_e2stim"

B7c_pl <- load.bigWig(paste(data_path,"/B7_CTRL_plus.bw", sep=""))
B7c_mn <- load.bigWig(paste(data_path,"/B7_CTRL_minus.bw", sep=""))
B7e_pl <- load.bigWig(paste(data_path,"/B7_E2_plus.bw", sep=""))
B7e_mn <- load.bigWig(paste(data_path,"/B7_E2_minus.bw", sep=""))

B7c_counts <- (B7c_pl$basesCovered*B7c_pl$mean)+abs(B7c_mn$basesCovered*B7c_mn$mean)
B7e_counts <- (B7e_pl$basesCovered*B7e_pl$mean)+abs(B7e_mn$basesCovered*B7e_mn$mean)


B7c <- bed6.region.bpQuery.bigWig(B7c_pl, B7c_mn, hah[,c(1:3,5:6,4)], abs.value=TRUE)
B7e2<- bed6.region.bpQuery.bigWig(B7e_pl, B7e_mn, hah[,c(1:3,5:6,4)], abs.value=TRUE)


b7 <- log((B7e2+1)/B7e_counts/((B7c+1)/B7c_counts), 2)
all_hah <- log((hah$E2_40m+1)/ (hah$VEH+1), 2)

use <- hah$E2.40m.qVal < 0.0001
plot(b7[use], all_hah[use]);abline(h=0);abline(v=0)

source("~/NHP/lib/densScatterplot.R")
densScatterplot(b7[use], all_hah[use]);abline(h=0);abline(v=0)

## Compare to G11
G11c_pl <- load.bigWig(paste(data_path,"/G11_CTRL_plus.bw", sep=""))
G11c_mn <- load.bigWig(paste(data_path,"/G11_CTRL_minus.bw", sep=""))
G11e_pl <- load.bigWig(paste(data_path,"/G11_E2_plus.bw", sep=""))
G11e_mn <- load.bigWig(paste(data_path,"/G11_E2_minus.bw", sep=""))

G11c_counts <- (G11c_pl$basesCovered*G11c_pl$mean)+abs(G11c_mn$basesCovered*G11c_mn$mean)
G11e_counts <- (G11e_pl$basesCovered*G11e_pl$mean)+abs(G11e_mn$basesCovered*G11e_mn$mean)

G11c <- bed6.region.bpQuery.bigWig(G11c_pl, G11c_mn, hah[,c(1:3,5:6,4)], abs.value=TRUE)
G11e2<- bed6.region.bpQuery.bigWig(G11e_pl, G11e_mn, hah[,c(1:3,5:6,4)], abs.value=TRUE)

g11 <- log((G11e2+1)/G11e_counts/((G11c+1)/G11c_counts), 2)
plot(b7[use], g11[use]);abline(h=0);abline(v=0); abline(0,1)

densScatterplot(b7[use], g11[use]);abline(h=0);abline(v=0)


