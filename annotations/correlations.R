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

path="../data/"
B7_pl <- load.bigWig(paste(path, "sB7_pl.bw", sep=""))
B7_mn <- load.bigWig(paste(path, "sB7_mn.bw", sep=""))
G11_pl <- load.bigWig(paste(path, "rG11_pl.bw", sep=""))
G11_mn <- load.bigWig(paste(path, "rG11_mn.bw", sep=""))
H9_pl <- load.bigWig(paste(path, "rH9_pl.bw", sep=""))
H9_mn <- load.bigWig(paste(path, "rH9_mn.bw", sep=""))
C11_pl <- load.bigWig(paste(path, "sC11_pl.bw", sep=""))
C11_mn <- load.bigWig(paste(path, "sC11_mn.bw", sep=""))
#hah_pl <- load.big

B7_BBCA_pl <- load.bigWig(paste(path,"sB7_bbca_pl.bw", sep=""))
B7_BBCA_mn <- load.bigWig(paste(path,"sB7_bbca_mn.bw", sep=""))
G11_BBCA_pl <- load.bigWig(paste(path,"rG11_bbca_pl.bw", sep=""))
G11_BBCA_mn <- load.bigWig(paste(path,"rG11_bbca_mn.bw", sep=""))

## E2 and TAM treated.
B7c_pl <- load.bigWig(paste(path, "B7_CTRL_plus.bw", sep=""))
B7c_mn <- load.bigWig(paste(path, "B7_CTRL_minus.bw", sep=""))
B7e_pl <- load.bigWig(paste(path, "B7_E2_plus.bw", sep=""))
B7e_mn <- load.bigWig(paste(path, "B7_E2_minus.bw", sep=""))
B7t_pl <- load.bigWig(paste(path, "B7_TAM_plus.bw", sep=""))
B7t_mn <- load.bigWig(paste(path, "B7_TAM_minus.bw", sep=""))
B7eb_pl <- load.bigWig(paste(path, "B7_E2_BBCA_plus.bw", sep=""))
B7eb_mn <- load.bigWig(paste(path, "B7_E2_BBCA_minus.bw", sep=""))
B7tb_pl<- load.bigWig(paste(path, "B7_TAM_BBCA_plus.bw", sep=""))
B7tb_mn<- load.bigWig(paste(path, "B7_TAM_BBCA_minus.bw", sep=""))

G11c_pl  <- load.bigWig(paste(path, "G11_CTRL_plus.bw", sep=""))
G11c_mn  <- load.bigWig(paste(path, "G11_CTRL_minus.bw", sep=""))
G11e_pl <- load.bigWig(paste(path, "G11_E2_plus.bw", sep=""))
G11e_mn <- load.bigWig(paste(path, "G11_E2_minus.bw", sep=""))
G11t_pl <- load.bigWig(paste(path, "G11_TAM_plus.bw", sep=""))
G11t_mn <- load.bigWig(paste(path, "G11_TAM_minus.bw", sep=""))
G11eb_pl<- load.bigWig(paste(path, "G11_E2_BBCA_plus.bw", sep=""))
G11eb_mn<- load.bigWig(paste(path, "G11_E2_BBCA_minus.bw", sep=""))
G11tb_pl<- load.bigWig(paste(path, "G11_TAM_BBCA_plus.bw", sep=""))
G11tb_mn<- load.bigWig(paste(path, "G11_TAM_BBCA_minus.bw", sep=""))

## Count reads in each ...
B7 <- bed6.region.bpQuery.bigWig(B7_pl, B7_mn, bodies, abs.value = TRUE)/ (bodies$V3-bodies$V2) * 1000/ (B7_pl$basesCovered * B7_pl$mean + abs(B7_mn$basesCovered * B7_mn$mean)) * 1e6
G11<- bed6.region.bpQuery.bigWig(G11_pl, G11_mn, bodies, abs.value = TRUE)/ (bodies$V3-bodies$V2) * 1000/ (G11_pl$basesCovered * G11_pl$mean + abs(G11_mn$basesCovered * G11_mn$mean)) * 1e6
H9 <- bed6.region.bpQuery.bigWig(H9_pl, H9_mn, bodies, abs.value = TRUE)/ (bodies$V3-bodies$V2) * 1000/ (H9_pl$basesCovered * H9_pl$mean + abs(H9_mn$basesCovered * H9_mn$mean)) * 1e6
C11<- bed6.region.bpQuery.bigWig(C11_pl, C11_mn, bodies, abs.value = TRUE)/ (bodies$V3-bodies$V2) * 1000/ (C11_pl$basesCovered * C11_pl$mean + abs(C11_mn$basesCovered * C11_mn$mean)) * 1e6
B7b <- bed6.region.bpQuery.bigWig(B7_BBCA_pl, B7_BBCA_mn, bodies, abs.value = TRUE)/ (bodies$V3-bodies$V2) * 1000/ (B7_BBCA_pl$basesCovered * B7_BBCA_pl$mean + abs(B7_BBCA_mn$basesCovered * B7_BBCA_mn$mean)) * 1e6
G11b<- bed6.region.bpQuery.bigWig(G11_BBCA_pl, G11_BBCA_mn, bodies, abs.value = TRUE)/ (bodies$V3-bodies$V2) * 1000/ (G11_BBCA_pl$basesCovered * G11_BBCA_pl$mean + abs(G11_BBCA_mn$basesCovered * G11_BBCA_mn$mean)) * 1e6

## Count E2TAM
B7c <- bed6.region.bpQuery.bigWig(B7c_pl, B7c_mn, bodies, abs.value = TRUE)/ (bodies$V3-bodies$V2) * 1000/ (B7c_pl$basesCovered * B7c_pl$mean + abs(B7c_mn$basesCovered * B7c_mn$mean)) * 1e6
B7e <- bed6.region.bpQuery.bigWig(B7e_pl, B7e_mn, bodies, abs.value = TRUE)/ (bodies$V3-bodies$V2) * 1000/ (B7e_pl$basesCovered * B7e_pl$mean + abs(B7e_mn$basesCovered * B7e_mn$mean)) * 1e6
B7t <- bed6.region.bpQuery.bigWig(B7t_pl, B7t_mn, bodies, abs.value = TRUE)/ (bodies$V3-bodies$V2) * 1000/ (B7t_pl$basesCovered * B7t_pl$mean + abs(B7t_mn$basesCovered * B7t_mn$mean)) * 1e6
B7eb<- bed6.region.bpQuery.bigWig(B7eb_pl, B7eb_mn, bodies, abs.value = TRUE)/ (bodies$V3-bodies$V2) * 1000/ (B7eb_pl$basesCovered * B7eb_pl$mean + abs(B7eb_mn$basesCovered * B7eb_mn$mean)) * 1e6
B7tb<- bed6.region.bpQuery.bigWig(B7tb_pl, B7tb_mn, bodies, abs.value = TRUE)/ (bodies$V3-bodies$V2) * 1000/ (B7tb_pl$basesCovered * B7tb_pl$mean + abs(B7tb_mn$basesCovered * B7tb_mn$mean)) * 1e6

G11c <- bed6.region.bpQuery.bigWig(G11c_pl, G11c_mn, bodies, abs.value = TRUE)/ (bodies$V3-bodies$V2) * 1000/ (G11c_pl$basesCovered * G11c_pl$mean + abs(G11c_mn$basesCovered * G11c_mn$mean)) * 1e6
G11e <- bed6.region.bpQuery.bigWig(G11e_pl, G11e_mn, bodies, abs.value = TRUE)/ (bodies$V3-bodies$V2) * 1000/ (G11e_pl$basesCovered * G11e_pl$mean + abs(G11e_mn$basesCovered * G11e_mn$mean)) * 1e6
G11t <- bed6.region.bpQuery.bigWig(G11t_pl, G11t_mn, bodies, abs.value = TRUE)/ (bodies$V3-bodies$V2) * 1000/ (G11t_pl$basesCovered * G11t_pl$mean + abs(G11t_mn$basesCovered * G11t_mn$mean)) * 1e6
G11eb<- bed6.region.bpQuery.bigWig(G11eb_pl, G11eb_mn, bodies, abs.value = TRUE)/ (bodies$V3-bodies$V2) * 1000/ (G11eb_pl$basesCovered * G11eb_pl$mean + abs(G11eb_mn$basesCovered * G11eb_mn$mean)) * 1e6
G11tb<- bed6.region.bpQuery.bigWig(G11tb_pl, G11tb_mn, bodies, abs.value = TRUE)/ (bodies$V3-bodies$V2) * 1000/ (G11tb_pl$basesCovered * G11tb_pl$mean + abs(G11tb_mn$basesCovered * G11tb_mn$mean)) * 1e6


gene_body_counts <- cbind(B7, C11, H9, G11, B7b, G11b, B7c, B7e, B7t, B7eb, B7tb, G11c, G11e, G11t, G11eb, G11tb)

## Make some plots
plot(B7, C11)
pairs(gene_body_counts, log="xy")
cor(gene_body_counts, method="spearman")

cor(gene_body_counts, method="spearman")

write.csv(cbind(refGene, gene_body_counts), "gene_body_counts.csv")

## Create a colorscatterplot
source("../lib/densScatterplot.R")
pdf("DensityScatterplots.pdf")

indx <- B7 > 0 & G11 > 0
densScatterplot(B7[indx], G11[indx], uselog=TRUE, xlab="B7", ylab="G11")

indx <- (B7c > 0 & B7e > 0 & B7eb > 0 & B7 > 0 & B7b > 0)
densScatterplot(B7[indx], B7b[indx], uselog=TRUE, xlab="B7", ylab="B7+BBCA")
densScatterplot(B7c[indx], B7e[indx], uselog=TRUE, xlab="B7 -- Control", ylab="B7 + E2")
densScatterplot(B7c[indx], B7eb[indx], uselog=TRUE, xlab="B7 -- Control", ylab="B7 + E2 + BBCA")

dev.off()

## Cluster changed genes...
yb.sig.pal <- function(n, scale=10) {
 ints<- c(0:(n-1))/(n-1)   ## Linear scale from 0:1 x N values.
 ints<- 1/(1+exp(scale*(0.5-ints)))## Transfer to sigmoidal scale.
 b<- min(ints)
 m<- 2*b/(n-1)
 ints<- ints+(m*c(0:(n-1)) -b)## Transform by linear function to fill colors out to maxes.
 
 ## Transfer to colorspace.
 # Yellow: 255, 255, 0
 # White:  255, 255, 255
 # Blue:   0, 0, 255
 YW <- ints[ints < 0.5] *2
 WB <- (ints[ints >= 0.5]-0.5) *2
 YW[YW<0] <- 0; WB[WB>1] <- 1
 c(rgb(1, 1, YW), rgb(1-WB, 1-WB, 1))
}

drawCor <- function(indx) {
	rpkm_df <- as.matrix(gene_body_counts[,indx]) # as.matrix(ca[,indx])#/(ca[,"mapSize"]) ## "Good?!"  Remove H2-U, H3-PI, C2-U+PI, M1-PI

	cond <- c(1,1,1,1,2,2,3,4,5,6,7,3,4,5,6,7)[indx]#"", "", "", "", "", "")# Condition[indx]
	spec <- c(1,1,2,2,1,2,1,1,1,1,1,2,2,2,2,2)[indx]#"", "", "", "", "", "")
	labs <- colnames(gene_body_counts)[indx]

	cc <- cor(rpkm_df, method="spearman")
	#clu <- agnes(t(rpkm_df))

	pal3 <- c("#E03CE9", "#17B92B", "#E6350D", "#6FD2F0", "#F9F77F", "#5B6C0C", "#68003D", "#310F08")
	
	## Print dendrogram and heatmap with latticeExtra.
	 library(latticeExtra)
	# hc1 <- agnes(1-cc, diss=TRUE, method="ward")
	# hc1 <- hclust(dist(t(rpkm_df), method = "canberra"))
	 hc1 <- hclust(dist(cc, method = "euclidean"),  method="single")
	 hc1 <- as.dendrogram(hc1)
	 ord.hc1 <- order.dendrogram(hc1)
	 hc2 <- reorder(hc1, cond[ord.hc1])
	 ord.hc2 <- order.dendrogram(hc2)
	 #region.colors <- trellis.par.get("superpose.polygon")$col

	 pl <- levelplot((cc)[ord.hc2, ord.hc2], col.regions= yb.sig.pal(100, scale=3), xlab="", ylab="", #rev(cm.colors(100)),  # #c("white", "yellow", "blue") # c("#E9F231", "#B1EC2C", "#5DBDEF")
		 colorkey = list(space="left", labels=list(cex=1.5)), 
		 scales = list(x= list(rot=90, cex=1.5, labels=labs[ord.hc2]), y=list(draw=FALSE)), #scales = list(x = list(rot = 90)), 
		 legend = list(
			right = list(fun = dendrogramGrob,
				 args = list(x = hc2, ord = ord.hc2, side = "right", #lwd=2,
				 size = 7, size.add = 0.5, 
				 add = list(rect = list(col = "transparent", fill = pal3[c(1:8)][cond])),
				 type = "rectangle")), 
			top = list(fun = dendrogramGrob,
				 args = list(x = hc2, ord = ord.hc2, side = "top", #lwd=2,
				 size = 1, size.add = 0.5, 
				 add = list(rect = list(col = "transparent", fill = pal3[2:8][spec])),
				 type = "rectangle"))
				 ))
	 print(pl)
}

pdf("correlationMatrix.pdf")
	drawCor(1:6)
	drawCor(1:4)
	drawCor(1:16)
dev.off()

