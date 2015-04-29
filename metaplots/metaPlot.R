##
## R script to generat meta plots around features of interest ... ER binding sites, and others.

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

doit <- function(bed, stp, halfWindow, ...) {
	bed <- bed[grep("Un|random", bed$V1, invert=TRUE),]
	bed <- center.bed(bed, upstreamWindow= halfWindow, downstreamWindow= halfWindow)
	bed_rev <- bed; bed_rev[bed[,6] == "+",6] <- "-"; bed_rev[bed[,6] == "-",6] <- "+"

	B7_meta_p <- metaprofile.bigWig(bed, B7_pl, B7_mn, step=stp)
	B7_meta_m <- rev(metaprofile.bigWig(bed_rev, B7_pl, B7_mn, step=stp))

        C11_meta_p <- metaprofile.bigWig(bed, C11_pl, C11_mn, step=stp)
        C11_meta_m <- rev(metaprofile.bigWig(bed_rev, C11_pl, C11_mn, step=stp))

        G11_meta_p <- metaprofile.bigWig(bed, G11_pl, G11_mn, step=stp)
        G11_meta_m <- rev(metaprofile.bigWig(bed_rev, G11_pl, G11_mn, step=stp))

        H9_meta_p <- metaprofile.bigWig(bed, H9_pl, H9_mn, step=stp)
        H9_meta_m <- rev(metaprofile.bigWig(bed_rev, H9_pl, H9_mn, step=stp))

        B7b_meta_p <- metaprofile.bigWig(bed, B7_pl, B7_mn, step=stp)
        B7b_meta_m <- rev(metaprofile.bigWig(bed_rev, B7_pl, B7_mn, step=stp))

        G11b_meta_p <- metaprofile.bigWig(bed, G11_BBCA_pl, G11_BBCA_mn, step=stp)
        G11b_meta_m <- rev(metaprofile.bigWig(bed_rev, G11_BBCA_pl, G11_BBCA_mn, step=stp))

        N = length(B7_meta_p$middle)
        x = ((1:N) - N/2)* stp #1:N*stp
        ylim=c(-1*max(c(B7_meta_m$middle, C11_meta_m$middle, G11_meta_m$middle, H9_meta_m$middle, B7b_meta_m$middle, G11b_meta_m$middle)), 
				max(B7_meta_p$middle, C11_meta_p$middle, G11_meta_p$middle, H9_meta_p$middle, B7b_meta_p$middle, G11b_meta_p$middle))
	
	par(mfrow=c(1,2))
	plot(-500, -500, ylim=ylim, xlim=c(min(x), max(x)), xlab= "Distance [bp]", ylab= "Signal.")
	lines(x, B7_meta_p$middle, col="#C14F4D")
	lines(x, -1*B7_meta_m$middle, col="#C14F4D")
	lines(x, C11_meta_p$middle, col="dark red")
	lines(x, -1*C11_meta_m$middle, col="dark red")
        lines(x, G11_meta_p$middle, col="#7E75B1")
        lines(x, -1*G11_meta_m$middle, col="#7E75B1")
        lines(x, H9_meta_p$middle, col="dark blue")
        lines(x, -1*H9_meta_m$middle, col="dark blue")

        plot(-500, -500, ylim=ylim, xlim=c(min(x), max(x)), xlab= "Distance [bp]", ylab= "Signal.")
        lines(x, B7_meta_p$middle, col="dark gray")
        lines(x, -1*B7_meta_m$middle, col="dark gray")
        lines(x, B7b_meta_p$middle, col="#98E676")
        lines(x, -1*B7b_meta_m$middle, col="#98E676")
        lines(x, G11_meta_p$middle, col="dark gray")
        lines(x, -1*G11_meta_m$middle, col="dark gray")
        lines(x, G11b_meta_p$middle, col="#9F97D4")
        lines(x, -1*G11b_meta_m$middle, col="#9F97D4")

# 	plot.metaprofile(B7_meta_p, minus.profile=B7_meta_m, X0=halfWindow/stp, ylim=ylim)
#       plot.metaprofile(C11_meta_p, minus.profile=C11_meta_m, X0=halfWindow/stp, ylim=ylim)
#       plot.metaprofile(G11_meta_p, minus.profile=G11_meta_m, X0=halfWindow/stp, ylim=ylim)
#       plot.metaprofile(H9_meta_p, minus.profile=H9_meta_m, X0=halfWindow/stp, ylim=ylim)

#       plot.metaprofile(B7_meta_p, minus.profile=B7_meta_m, X0=halfWindow/stp, ylim=ylim)
#       plot.metaprofile(B7b_meta_p, minus.profile=B7b_meta_m, X0=halfWindow/stp, ylim=ylim)
#       plot.metaprofile(G11_meta_p, minus.profile=G11_meta_m, X0=halfWindow/stp, ylim=ylim)
#       plot.metaprofile(G11b_meta_p, minus.profile=G11b_meta_m, X0=halfWindow/stp, ylim=ylim)

        #plot.metaprofile(H_meta_p, minus.profile=H_meta_m, X0=halfWindow/stp, ylim=ylim)
	#lines(x, C_meta_p$middle, col="#17b92b")
	#lines(x, -1*C_meta_m$middle, col="#17b92b")
        #lines(x, M_meta_p$middle, col="#5b6c0c")
        #lines(x, -1*M_meta_m$middle, col="#5b6c0c")
}

## Read in BED files and print...
erbs <- read.table("ERaBS.bed")
doit(erbs, stp=20, halfWindow=2000)

