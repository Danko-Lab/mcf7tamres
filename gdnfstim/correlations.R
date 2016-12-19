## get correlations between MCF-7 samples

load("Counts.RData")

## Make some plots
plot(B7, C11)
pairs(gene_body_density, log="xy")
cor(gene_body_density, method="spearman")
cor(gene_body_density, method="spearman")

## Create a colorscatterplot
source("../lib/densScatterplot.R")
pdf("S5A.DensityScatterplots.pdf")

rpkm <- gene_body_density

densScatterplot(rpkm[,"B7_0h"], rpkm[,"C11_0h"], uselog=TRUE, xlab="PRO-seq, B7 0h [RPKM]", ylab="PRO-seq, C11 0h [rpkm]", main="0h GDNF")
densScatterplot(rpkm[,"B7_1h"], rpkm[,"C11_1h"], uselog=TRUE, xlab="PRO-seq, B7 1h [RPKM]", ylab="PRO-seq, C11 1h [rpkm]", main="1h GDNF")
densScatterplot(rpkm[,"B7_24h"], rpkm[,"C11_24h"], uselog=TRUE, xlab="PRO-seq, B7 24h [RPKM]", ylab="PRO-seq, C11 24h [rpkm]", main="24h GDNF")

densScatterplot(rpkm[,"G11_0h"], rpkm[,"H9_0h"], uselog=TRUE, xlab="PRO-seq, G11 0h [RPKM]", ylab="PRO-seq, H9 0h [rpkm]", main="0h GDNF")
densScatterplot(rpkm[,"G11_1h"], rpkm[,"H9_1h"], uselog=TRUE, xlab="PRO-seq, G11 1h [RPKM]", ylab="PRO-seq, H9 1h [rpkm]", main="1h GDNF")
densScatterplot(rpkm[,"G11_24h"], rpkm[,"H9_24h"], uselog=TRUE, xlab="PRO-seq, G11 24h [RPKM]", ylab="PRO-seq, H9 24h [rpkm]", main="24h GDNF")

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
 YW <- ints[ints < 0.25] *4
 WB <- (ints[ints >= 0.5]-0.5) *2
 YW[YW<0] <- 0; WB[WB>1] <- 1
 c(rgb(1, 1, YW), rgb(1-WB, 1-WB, 1))
}

drawCor <- function(indx) {
	rpkm_df <- as.matrix(gene_body_density[,indx]) # as.matrix(ca[,indx])#/(ca[,"mapSize"]) ## "Good?!"  Remove H2-U, H3-PI, C2-U+PI, M1-PI

	cond <- c(1,1,1,1,1,1,2,2,2,2,2,2)[indx]#"", "", "", "", "", "")# Condition[indx]
	spec <- c(1,2,3,1,2,3,1,2,3,1,2,3)[indx]#"", "", "", "", "", "")
	labs <- colnames(gene_body_density)[indx]

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

pdf("S5B.correlationMatrix.pdf")
	drawCor(1:12)
dev.off()

