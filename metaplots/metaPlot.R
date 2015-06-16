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
C11_pl <- load.bigWig(paste(path, "5089_5598_18733_H73KFGXX_MCF7_C11_TGACCA_R1_plus.bw",sep=""))
C11_mn <- load.bigWig(paste(path, "5089_5598_18734_H73KFGXX_MCF7_C11_ACAGTG_R1_minus.bw",sep=""))

B7_BBCA_pl <- load.bigWig(paste(path, "5089_5598_18730_H73KFGXX_MCF7_B7_BBCA_ATCACG_R1_plus.bw",sep=""))
B7_BBCA_mn <- load.bigWig(paste(path, "5089_5598_18730_H73KFGXX_MCF7_B7_BBCA_ATCACG_R1_minus.bw",sep=""))
G11_BBCA_pl <- load.bigWig(paste(path, "5089_5598_18731_H73KFGXX_MCF7_G11_BBCA_CGATGT_R1_plus.bw", sep=""))
G11_BBCA_mn <- load.bigWig(paste(path, "5089_5598_18731_H73KFGXX_MCF7_G11_BBCA_CGATGT_R1_minus.bw", sep=""))

doit <- function(bed, stp, halfWindow, ...) {
	bed <- bed[grep("Un|random", bed$V1, invert=TRUE),]
	bed <- center.bed(bed, upstreamWindow= halfWindow, downstreamWindow= halfWindow)
	bed_rev <- bed; bed_rev[bed[,6] == "+",6] <- "-"; bed_rev[bed[,6] == "-",6] <- "+"

	B7_meta_p <- metaprofile.bigWig(bed, B7_pl, B7_mn, step=stp)$middle/ stp / (B7_pl$basesCovered * B7_pl$mean + abs(B7_mn$basesCovered * B7_mn$mean)) * 1e6
	B7_meta_m <- metaprofile.bigWig(bed_rev, B7_pl, B7_mn, step=stp)$middle/ stp/ (B7_pl$basesCovered * B7_pl$mean + abs(B7_mn$basesCovered * B7_mn$mean)) * 1e6

        C11_meta_p <- metaprofile.bigWig(bed, C11_pl, C11_mn, step=stp)$middle/ stp/ (C11_pl$basesCovered * C11_pl$mean + abs(C11_mn$basesCovered * C11_mn$mean)) * 1e6
        C11_meta_m <- metaprofile.bigWig(bed_rev, C11_pl, C11_mn, step=stp)$middle/ stp/ (C11_pl$basesCovered * C11_pl$mean + abs(C11_mn$basesCovered * C11_mn$mean)) * 1e6

        G11_meta_p <- metaprofile.bigWig(bed, G11_pl, G11_mn, step=stp)$middle/ stp/ (G11_pl$basesCovered * G11_pl$mean + abs(G11_mn$basesCovered * G11_mn$mean)) * 1e6
        G11_meta_m <- metaprofile.bigWig(bed_rev, G11_pl, G11_mn, step=stp)$middle/ stp/ (G11_pl$basesCovered * G11_pl$mean + abs(G11_mn$basesCovered * G11_mn$mean)) * 1e6

        H9_meta_p <- metaprofile.bigWig(bed, H9_pl, H9_mn, step=stp)$middle/ stp/ (H9_pl$basesCovered * H9_pl$mean + abs(H9_mn$basesCovered * H9_mn$mean)) * 1e6
        H9_meta_m <- metaprofile.bigWig(bed_rev, H9_pl, H9_mn, step=stp)$middle/ stp/ (H9_pl$basesCovered * H9_pl$mean + abs(H9_mn$basesCovered * H9_mn$mean)) * 1e6

        B7b_meta_p <- metaprofile.bigWig(bed, B7_pl, B7_mn, step=stp)$middle/ stp/ (B7_BBCA_pl$basesCovered * B7_BBCA_pl$mean + abs(B7_BBCA_mn$basesCovered * B7_BBCA_mn$mean)) * 1e6
        B7b_meta_m <- metaprofile.bigWig(bed_rev, B7_pl, B7_mn, step=stp)$middle/ stp/ (B7_BBCA_pl$basesCovered * B7_BBCA_pl$mean + abs(B7_BBCA_mn$basesCovered * B7_BBCA_mn$mean)) * 1e6

        G11b_meta_p <- metaprofile.bigWig(bed, G11_BBCA_pl, G11_BBCA_mn, step=stp)$middle/ stp/ (G11_BBCA_pl$basesCovered * G11_BBCA_pl$mean + abs(G11_BBCA_mn$basesCovered * G11_BBCA_mn$mean)) * 1e6
        G11b_meta_m <- metaprofile.bigWig(bed_rev, G11_BBCA_pl, G11_BBCA_mn, step=stp)$middle/ stp/ (G11_BBCA_pl$basesCovered * G11_BBCA_pl$mean + abs(G11_BBCA_mn$basesCovered * G11_BBCA_mn$mean)) * 1e6

        N = length(B7_meta_p)
        x = ((1:N) - N/2)* stp #1:N*stp
        ylim=c(-1*max(c(B7_meta_m, C11_meta_m, G11_meta_m, H9_meta_m, B7b_meta_m, G11b_meta_m)), 
				max(B7_meta_p, C11_meta_p, G11_meta_p, H9_meta_p, B7b_meta_p, G11b_meta_p))
	
	par(mfrow=c(1,2))
	plot(-500, -500, ylim=ylim, xlim=c(min(x), max(x)), xlab= "Distance [bp]", ylab= "Signal.", main="TamS v. TamR")
	lines(x, B7_meta_p, col="#C14F4D")
	lines(x, rev(-1*B7_meta_m), col="#C14F4D")
	lines(x, C11_meta_p, col="dark red")
	lines(x, rev(-1*C11_meta_m), col="dark red")
        lines(x, G11_meta_p, col="#7E75B1")
        lines(x, rev(-1*G11_meta_m), col="#7E75B1")
        lines(x, H9_meta_p, col="dark blue")
        lines(x, rev(-1*H9_meta_m), col="dark blue")

        plot(-500, -500, ylim=ylim, xlim=c(min(x), max(x)), xlab= "Distance [bp]", ylab= "Signal.", main="BBCA/ DMSO")
        lines(x, B7_meta_p, col="dark gray")
        lines(x, rev(-1*B7_meta_m), col="dark gray")
        lines(x, B7b_meta_p, col="black")
        lines(x, rev(-1*B7b_meta_m), col="black")
        lines(x, G11_meta_p, col="pink")
        lines(x, rev(-1*G11_meta_m), col="pink")
        lines(x, G11b_meta_p, col="red")
        lines(x, rev(-1*G11b_meta_m), col="red")

        par(mfrow=c(1,1))
 	plot.metaprofile(B7_meta_p, minus.profile=B7_meta_m, X0=halfWindow/stp, ylim=ylim)
        plot.metaprofile(C11_meta_p, minus.profile=C11_meta_m, X0=halfWindow/stp, ylim=ylim)
        plot.metaprofile(G11_meta_p, minus.profile=G11_meta_m, X0=halfWindow/stp, ylim=ylim)
        plot.metaprofile(H9_meta_p, minus.profile=H9_meta_m, X0=halfWindow/stp, ylim=ylim)

        plot.metaprofile(B7_meta_p, minus.profile=B7_meta_m, X0=halfWindow/stp, ylim=ylim)
        plot.metaprofile(B7b_meta_p, minus.profile=B7b_meta_m, X0=halfWindow/stp, ylim=ylim)
        plot.metaprofile(G11_meta_p, minus.profile=G11_meta_m, X0=halfWindow/stp, ylim=ylim)
        plot.metaprofile(G11b_meta_p, minus.profile=G11b_meta_m, X0=halfWindow/stp, ylim=ylim)
}

## Read in BED files and print...
erbs <- read.table("ERaBS.bed")
pdf("ERBS.meta.pdf")
doit(erbs, stp=20, halfWindow=2000)
dev.off()

