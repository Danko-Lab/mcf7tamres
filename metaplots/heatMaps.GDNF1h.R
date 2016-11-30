##
## Creates heatmaps on ER-alpha binding sites.
## Ultimately, divides heatmaps: E2/ CTRL, TAM/ CTRL, E2+BBCA/ CTRL.

require(bigWig)
library(pheatmap)

pth= "/local/storage/projects/mcf7tamres/data/"
dist <- 10000

## Load dREG-HD sites of changed promoters.
genes <- read.table("../gdnfstim/GDNF_treatment.tsv", header=TRUE)
genes[genes$TXSTRAND=="+","TXEND"] <- genes[genes$TXSTRAND=="+","TXSTART"]+1
genes[genes$TXSTRAND=="-","TXSTART"] <- genes[genes$TXSTRAND=="-","TXEND"]-1

hdu <- genes[!is.na(genes$FDR_GDNF_1h_S) & genes$FDR_GDNF_1h_S < 0.01 & genes$FC_GDNF_1h_S > 0,]
hdd <- genes[!is.na(genes$FDR_GDNF_1h_S) & genes$FDR_GDNF_1h_S < 0.01 & genes$FC_GDNF_1h_S < 0,]
nrg <- genes[!is.na(genes$FDR_GDNF_1h_S) & genes$FDR_GDNF_1h_S > 0.3  & abs(genes$FC_GDNF_1h_S) < 0.25,]

hdu  <- hdu[,1:6]
hdd  <- hdd[,1:6]
nrg  <- nrg[,1:6]

calcHeatmap<- function(hd, hMarkFile1p, hMarkFile2p, hMarkFile1m, hMarkFile2m, name, hm_order= NULL, step=100, path=pth, cols= c("#0000FE", "white","#FE0000"), hmcols= NULL, bk= NULL) {
	print("Counting Reads")

	## Load mark.
	hMark1p <- load.bigWig(paste(path, hMarkFile1p, sep=""))  #"/local/storage/data/hg19/cd4/epiRoadmap_histone/H3K27ac.bw")
        hMark2p <- load.bigWig(paste(path, hMarkFile2p, sep=""))  #"/local/storage/data/hg19/cd4/epiRoadmap_histone/H3K27ac.bw")
        hMark1m <- load.bigWig(paste(path, hMarkFile1m, sep=""))  #"/local/storage/data/hg19/cd4/epiRoadmap_histone/H3K27ac.bw")
        hMark2m <- load.bigWig(paste(path, hMarkFile2m, sep=""))  #"/local/storage/data/hg19/cd4/epiRoadmap_histone/H3K27ac.bw")

	## Get a matrix of counts.
	hCountMatrix1p <- bed6.step.bpQuery.bigWig(hMark1p, hMark1m, center.bed(hd, dist/5, dist), step=step, abs.value=TRUE, follow.strand=TRUE)
        hCountMatrix2p <- bed6.step.bpQuery.bigWig(hMark2p, hMark2m, center.bed(hd, dist/5, dist), step=step, abs.value=TRUE, follow.strand=TRUE)

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
	navg <- 5 ## Average every navg rows
	avgMat <- t(sapply(1:floor(NROW(hmat)/navg), function(x) {colMeans(hmat[((x-1)*navg+1):min(NROW(hmat),(x*navg)),])}))
	hmat <- avgMat

	## Write out a heatmap.
	library(pheatmap)
	maxVal <- log(max(c(1/exp(min(hmat)), exp(max(hmat))))); print(paste("maxVal", maxVal))
	if(is.null(hmcols) | is.null(bk)) {
		print("Computing Heatmap Colors")
		bk <- c(seq(-1*maxVal, maxVal, 0.01)) #seq(min(hmat), max(hmat), 0.01)
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

getOrder <- function(hd, b1p, b1m, b2p, b2m) {
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

ord_b7_GDNF1h <- getOrder(hdu, "MCF-7_B7_GDNF_1_hr_plus.bw", "MCF-7_B7_GDNF_1_hr_minus.bw", "MCF-7_B7_GDNF_0_hr_plus.bw", "MCF-7_B7_GDNF_0_hr_minus.bw")
ord_b7_GDNF1h_dn <- getOrder(hdd, "MCF-7_B7_GDNF_1_hr_plus.bw", "MCF-7_B7_GDNF_1_hr_minus.bw", "MCF-7_B7_GDNF_0_hr_plus.bw", "MCF-7_B7_GDNF_0_hr_minus.bw")
ord_b7_GDNF1h_nrg <- getOrder(nrg, "MCF-7_B7_GDNF_1_hr_plus.bw", "MCF-7_B7_GDNF_1_hr_minus.bw", "MCF-7_B7_GDNF_0_hr_plus.bw", "MCF-7_B7_GDNF_0_hr_minus.bw")

## Heatmaps

#B7
plus_E2 <- calcHeatmap(hdu, "MCF-7_B7_GDNF_1_hr_plus.bw", "MCF-7_B7_GDNF_0_hr_plus.bw", "MCF-7_B7_GDNF_1_hr_minus.bw", "MCF-7_B7_GDNF_0_hr_minus.bw", "B7-GDNF-genes-up", hm_order=ord_b7_GDNF1h)
plus_E2 <- calcHeatmap(hdd, "MCF-7_B7_GDNF_1_hr_plus.bw", "MCF-7_B7_GDNF_0_hr_plus.bw", "MCF-7_B7_GDNF_1_hr_minus.bw", "MCF-7_B7_GDNF_0_hr_minus.bw", "B7-GDNF-genes-dn", hm_order=ord_b7_GDNF1h_dn)
nrg_E2 <- calcHeatmap(nrg, "MCF-7_B7_GDNF_1_hr_plus.bw", "MCF-7_B7_GDNF_0_hr_plus.bw", "MCF-7_B7_GDNF_1_hr_minus.bw", "MCF-7_B7_GDNF_0_hr_minus.bw", "B7-GDNF-genes-nonreg", hm_order=ord_b7_GDNF1h_nrg)

plus_E2 <- calcHeatmap(hdu, "MCF-7_C11_GDNF_1_hr_plus.bw", "MCF-7_C11_GDNF_0_hr_plus.bw", "MCF-7_C11_GDNF_1_hr_minus.bw", "MCF-7_C11_GDNF_0_hr_minus.bw", "C11-GDNF-genes-up", hm_order=ord_b7_GDNF1h)
plus_E2 <- calcHeatmap(hdd, "MCF-7_C11_GDNF_1_hr_plus.bw", "MCF-7_C11_GDNF_0_hr_plus.bw", "MCF-7_C11_GDNF_1_hr_minus.bw", "MCF-7_C11_GDNF_0_hr_minus.bw", "C11-GDNF-genes-dn", hm_order=ord_b7_GDNF1h_dn)



