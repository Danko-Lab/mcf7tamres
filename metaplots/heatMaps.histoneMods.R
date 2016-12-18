##
## Creates heatmaps on human dREG-HD sites.

require(bigWig)

dist <- 10000

## Load dREG-HD sites.
hd <- read.table("../dreg/dREG.HD.rerun.Nov.30.2016/mcf7.bedgraph.gz.bed.gz_dREG_HD_relaxed.bed"); hd <- hd[hd[,2]-dist > 0,]; hd <- hd[hd$V1 != "chrY",]
print(NROW(hd))

writeHeatmap<- function(hMarkFile, name, hm_order= NULL, 
							cols= c("white","#00A63E"), dist= 10000, step=100,
							path="/local/storage/data/hg19/mcf7/histones/") {
	## Load mark.
	hMark <- load.bigWig(paste(path, hMarkFile, sep=""))  #"/local/storage/data/hg19/cd4/epiRoadmap_histone/H3K27ac.bw")

	## Get a matrix of counts.
	hCountMatrix <- bed.step.bpQuery.bigWig(hMark, center.bed(hd, dist, dist), step=step, abs.value=TRUE)
	hmat <- log(matrix(unlist(hCountMatrix), nrow= NROW(hd), byrow=TRUE)+1)
	if(is.null(hm_order)) {
	  hm_order <- order(rowSums(hmat[,(NCOL(hmat)/2 -10):(NCOL(hmat)/2 +10)]), decreasing=TRUE)
	}
	hmat <- hmat[hm_order,]

	## Average by rows of 10.
	navg <- 5 ## Average every navg rows
	avgMat <- t(sapply(1:floor(NROW(hmat)/navg), function(x) {colMeans(hmat[((x-1)*navg+1):min(NROW(hmat),(x*navg)),])}))
	hmat <- avgMat

	## Write out a heatmap.
	library(pheatmap)
	bk <- seq(min(hmat), max(hmat), 0.01)
	hmcols <- colorRampPalette(cols)(length(bk)-1) # red

	png(paste(name,".png",sep=""), width=250, height = 400)     # width and height are in pixels
		pheatmap( hmat, cluster_rows = FALSE, cluster_cols = FALSE, col= hmcols, breaks = bk, legend=FALSE, show_rownames=FALSE, show_colnames=FALSE )
	dev.off()
	
	return(hm_order)
}

ord <- writeHeatmap("wgEncodeSydhHistoneMcf7H3k27acUcdSig.bigWig", "H3K27ac", cols=c("white","#00A63E"))
sup <- writeHeatmap("GSM945269_hg19_wgEncodeUwHistoneMcf7H3k4me3StdRawRep1.bigWig", "H3K4me3", hm_order= ord, cols=c("white","#fe0000"))
sup <- writeHeatmap("wgEncodeSydhHistoneMcf7H3k36me3bUcdSig.bigWig", "H3K36me3", hm_order= ord, cols=c("white","#0022fe"))
sup <- writeHeatmap("GSM1024784_hg19_wgEncodeUwDnaseMcf7Est100nm1hRawRep1.bigWig", "DNase-I", hm_order= ord, cols=c("white","#000000"), path="/local/storage/data/hg19/mcf7/dnase/")

pth= "/local/storage/projects/mcf7tamres/dreg/"
sup <- writeHeatmap("mcf7.plus.bw", "PROseq.plus", hm_order= ord, cols=c("white","#fe0000"), path=pth)
sup <- writeHeatmap("mcf7.minus.bw", "PROseq.minus", hm_order= ord, cols=c("white","#0000fe"), path=pth)


