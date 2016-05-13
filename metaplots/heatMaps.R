##
## Creates heatmaps on ER-alpha binding sites.
## Ultimately, divides heatmaps: E2/ CTRL, TAM/ CTRL, E2+BBCA/ CTRL.

require(bigWig)
library(pheatmap)

pth= "/local/storage/projects/mcf7tamres/data/"
dist <- 4000

## Load dREG-HD sites.
hd <- read.table("ERaBS.bed"); hd <- hd[hd[,2]-dist > 0,]

calcHeatmap<- function(hMarkFile1p, hMarkFile2p, hMarkFile1m, hMarkFile2m, name, hm_order= NULL, step=50, path=pth, cols= c("white","#00A63E"), hmcols= NULL, bk= NULL) {
	print("Counting Reads")

	## Load mark.
	hMark1p <- load.bigWig(paste(path, hMarkFile1p, sep=""))  #"/local/storage/data/hg19/cd4/epiRoadmap_histone/H3K27ac.bw")
        hMark2p <- load.bigWig(paste(path, hMarkFile2p, sep=""))  #"/local/storage/data/hg19/cd4/epiRoadmap_histone/H3K27ac.bw")
        hMark1m <- load.bigWig(paste(path, hMarkFile1m, sep=""))  #"/local/storage/data/hg19/cd4/epiRoadmap_histone/H3K27ac.bw")
        hMark2m <- load.bigWig(paste(path, hMarkFile2m, sep=""))  #"/local/storage/data/hg19/cd4/epiRoadmap_histone/H3K27ac.bw")

	## Get a matrix of counts.
	hCountMatrix1p <- bed.step.bpQuery.bigWig(hMark1p, center.bed(hd, dist, dist), step=step, abs.value=TRUE)
        hCountMatrix2p <- bed.step.bpQuery.bigWig(hMark2p, center.bed(hd, dist, dist), step=step, abs.value=TRUE)

	## Sizes.
	size1 <- (abs(hMark1p$mean)*hMark1p$basesCovered+abs(hMark1m$mean)*hMark1m$basesCovered)
	size2 <- (abs(hMark2p$mean)*hMark2p$basesCovered+abs(hMark2m$mean)*hMark2m$basesCovered)

	hmat1 <- matrix(unlist(hCountMatrix1p), nrow= NROW(hd), byrow=TRUE)
	hmat2 <- matrix(unlist(hCountMatrix2p), nrow= NROW(hd), byrow=TRUE)
	hmat  <- log((hmat1+1)/(hmat2+1)*(size2/size1))

	if(is.null(hm_order)) {
		print("Computing Heatmap Order")
		hm_order <- order(rowSums(hmat[,(NCOL(hmat)/2 -20):(NCOL(hmat)/2 +20)]), decreasing=TRUE)
	}
	hmat <- hmat[hm_order,]

	## Average by rows of 150.
	navg <- 100 ## Average every navg rows
	avgMat <- t(sapply(1:floor(NROW(hmat)/navg), function(x) {colMeans(hmat[((x-1)*navg+1):min(NROW(hmat),(x*navg)),])}))
	hmat <- avgMat

	## Write out a heatmap.
	library(pheatmap)
	if(is.null(hmcols) | is.null(bk)) {
		print("Computing Heatmap Colors")
		bk <- c(log(seq(1, 3, 0.01))) #seq(min(hmat), max(hmat), 0.01)
		hmcols <- colorRampPalette(cols, bias=1)(length(bk)-1) # red
		png(paste("heatmaps/", name,".scale.png",sep=""), width=400, height=500)
			pheatmap(rev(bk), cluster_rows = FALSE, cluster_cols = FALSE, col= hmcols, breaks = bk, legend=TRUE, legend_breaks= quantile(bk), legend_labels= signif(exp(quantile(bk)),3), show_rownames=FALSE, show_colnames=FALSE)
		dev.off()
		print(summary(exp(bk)))
	}

	png(paste("heatmaps/",name,".png",sep=""), width=400, height = 600)     # width and height are in pixels
		pheatmap( hmat, cluster_rows = FALSE, cluster_cols = FALSE, col= hmcols, breaks = bk, legend=FALSE, show_rownames=FALSE, show_colnames=FALSE )
	dev.off()
	
        return(list(matrix=hmat, order=hm_order, hmcols= hmcols, bk= bk))
}

#get order: 

getOrder <- function(b1p, b1m, b2p, b2m) {
 pbw1 <- load.bigWig(paste(pth, b1p, sep=""))
 mbw1 <- load.bigWig(paste(pth, b1m, sep=""))
 pbw2 <- load.bigWig(paste(pth, b2p, sep=""))
 mbw2 <- load.bigWig(paste(pth, b2m, sep=""))

 size1 <- (pbw1$mean*pbw1$basesCovered+abs(mbw1$mean)*mbw1$basesCovered)
 size2 <- (pbw2$mean*pbw2$basesCovered+abs(mbw2$mean)*mbw2$basesCovered)

## Get order based on GRO-cap plus and minus strand.
 p_counts_1 <- matrix(unlist(bed.step.bpQuery.bigWig(pbw1, center.bed(hd, dist, dist), step=100, abs.value=TRUE)), nrow= NROW(hd), byrow=TRUE)
 m_counts_1 <- matrix(unlist(bed.step.bpQuery.bigWig(mbw1, center.bed(hd, dist, dist), step=100, abs.value=TRUE)), nrow= NROW(hd), byrow=TRUE)
 p_counts_2 <- matrix(unlist(bed.step.bpQuery.bigWig(pbw2, center.bed(hd, dist, dist), step=100, abs.value=TRUE)), nrow= NROW(hd), byrow=TRUE)
 m_counts_2 <- matrix(unlist(bed.step.bpQuery.bigWig(mbw2, center.bed(hd, dist, dist), step=100, abs.value=TRUE)), nrow= NROW(hd), byrow=TRUE)

 hmat <- log((p_counts_1+m_counts_1+1)/ (p_counts_2+m_counts_2+1) * (size2/size1))
 ord <- order(rowSums(hmat[,(NCOL(hmat)/2 -20):(NCOL(hmat)/2 +20)]), decreasing=TRUE)
}

ord_b7_E2 <- getOrder("B7_E2_plus.bw", "B7_E2_minus.bw", "B7_CTRL_plus.bw", "B7_CTRL_minus.bw")
ord_b7_TAM<- getOrder("B7_TAM_plus.bw", "B7_TAM_minus.bw", "B7_CTRL_plus.bw", "B7_CTRL_minus.bw")
ord_b7_BB <- getOrder("B7_E2_BBCA_plus.bw", "B7_E2_BBCA_minus.bw", "B7_E2_BBCA_plus.bw", "B7_E2_BBCA_minus.bw")

ord_g11_E2<- getOrder("G11_E2_plus.bw", "G11_E2_minus.bw", "G11_CTRL_plus.bw", "G11_CTRL_minus.bw")
ord_g11_TAM<-getOrder("G11_TAM_plus.bw", "G11_TAM_minus.bw", "G11_CTRL_plus.bw", "G11_CTRL_minus.bw")
ord_g11_BB<- getOrder("G11_E2_BBCA_plus.bw", "G11_E2_BBCA_minus.bw", "G11_E2_BBCA_plus.bw", "G11_E2_BBCA_minus.bw")

## Heatmaps

#B7
plus_E2 <- calcHeatmap("B7_E2_plus.bw", "B7_CTRL_plus.bw", "B7_E2_minus.bw", "B7_CTRL_minus.bw", "B7_e2_pl", col=c("white","#fe0000"), hm_order=ord_b7_E2)
minus_E2<- calcHeatmap("B7_E2_minus.bw", "B7_CTRL_minus.bw", "B7_E2_plus.bw", "B7_CTRL_plus.bw", "B7_e2_mn", col=c("white","#0000fe"), hm_order=ord_b7_E2)

tam_plus_E2 <- calcHeatmap("B7_TAM_plus.bw", "B7_CTRL_plus.bw", "B7_TAM_minus.bw", "B7_CTRL_minus.bw", "B7_tam_pl", hmcols=plus_E2$hmcols, bk= plus_E2$bk, hm_order=ord_b7_TAM)
tam_minus_E2<- calcHeatmap("B7_TAM_minus.bw", "B7_CTRL_minus.bw", "B7_TAM_plus.bw", "B7_CTRL_plus.bw", "B7_tam_mn", hmcols=minus_E2$hmcols, bk= minus_E2$bk,hm_order=ord_b7_TAM)

bbca_plus_E2 <- calcHeatmap("B7_E2_BBCA_plus.bw", "B7_CTRL_plus.bw", "B7_E2_BBCA_minus.bw", "B7_CTRL_minus.bw", "B7_bb_pl", hmcols=plus_E2$hmcols, bk= plus_E2$bk, hm_order=ord_b7_BB)
bbca_minus_E2<- calcHeatmap("B7_E2_BBCA_minus.bw", "B7_CTRL_minus.bw", "B7_E2_BBCA_plus.bw", "B7_CTRL_plus.bw", "B7_bb_mn", hmcols=minus_E2$hmcols, bk= minus_E2$bk,hm_order=ord_b7_BB)


#G11
g11_plus_E2 <- calcHeatmap("G11_E2_plus.bw", "G11_CTRL_plus.bw", "G11_E2_minus.bw", "G11_CTRL_minus.bw", "G11_e2_pl", hm_order=ord_g11_E2, hmcols=plus_E2$hmcols, bk= plus_E2$bk)
g11_minus_E2<- calcHeatmap("G11_E2_minus.bw", "G11_CTRL_minus.bw", "G11_E2_plus.bw", "G11_CTRL_plus.bw", "G11_e2_mn", hm_order=ord_g11_E2, hmcols=minus_E2$hmcols, bk= minus_E2$bk)

g11_tam_plus_E2 <- calcHeatmap("G11_TAM_plus.bw", "G11_CTRL_plus.bw", "G11_TAM_minus.bw", "G11_CTRL_minus.bw", "G11_tam_pl", hmcols=plus_E2$hmcols, bk= plus_E2$bk, hm_order=ord_g11_TAM)
g11_tam_minus_E2<- calcHeatmap("G11_TAM_minus.bw", "G11_CTRL_minus.bw", "G11_TAM_plus.bw", "G11_CTRL_plus.bw", "G11_tam_mn", hmcols=minus_E2$hmcols, bk= minus_E2$bk,hm_order=ord_g11_TAM)

g11_bbca_plus_E2 <- calcHeatmap("G11_E2_BBCA_plus.bw", "G11_CTRL_plus.bw", "G11_E2_BBCA_minus.bw", "G11_CTRL_minus.bw", "G11_bb_pl", hmcols=plus_E2$hmcols, bk= plus_E2$bk, hm_order=ord_g11_BB)
g11_bbca_minus_E2<- calcHeatmap("G11_E2_BBCA_minus.bw", "G11_CTRL_minus.bw", "G11_E2_BBCA_plus.bw", "G11_CTRL_plus.bw", "G11_bb_mn", hmcols=minus_E2$hmcols, bk= minus_E2$bk,hm_order=ord_g11_BB)


