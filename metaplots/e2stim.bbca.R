##
## Looking at the effects of BBCA on E2-induced eRNAs.

require(bigWig)

path="/local/storage/projects/mcf7tamres/data/"
B7_c_pl  <- load.bigWig(paste(path, "B7_CTRL_plus.bw", sep=""))
B7_c_mn  <- load.bigWig(paste(path, "B7_CTRL_minus.bw", sep=""))
B7_e2_pl <- load.bigWig(paste(path, "B7_E2_plus.bw", sep=""))
B7_e2_mn <- load.bigWig(paste(path, "B7_E2_minus.bw", sep=""))
B7_e2b_pl<- load.bigWig(paste(path, "B7_E2_BBCA_plus.bw", sep=""))
B7_e2b_mn<- load.bigWig(paste(path, "B7_E2_BBCA_minus.bw", sep=""))

G11_c_pl  <- load.bigWig(paste(path, "G11_CTRL_plus.bw", sep=""))
G11_c_mn  <- load.bigWig(paste(path, "G11_CTRL_minus.bw", sep=""))
G11_e2_pl <- load.bigWig(paste(path, "G11_E2_plus.bw", sep=""))
G11_e2_mn <- load.bigWig(paste(path, "G11_E2_minus.bw", sep=""))
G11_e2b_pl<- load.bigWig(paste(path, "G11_E2_BBCA_plus.bw", sep=""))
G11_e2b_mn<- load.bigWig(paste(path, "G11_E2_BBCA_minus.bw", sep=""))

doit <- function(bed, stp, halfWindow, c_pl, c_mn, e_pl, e_mn, e_b_pl, e_b_mn, ...) {
        bed <- bed[grep("Un|random", bed$V1, invert=TRUE),]
        bed <- center.bed(bed, upstreamWindow= halfWindow, downstreamWindow= halfWindow)
	
        c_p <- metaprofile.bigWig(bed, c_pl, c_mn, step=stp)$middle/ stp / (c_pl$basesCovered * c_pl$mean + abs(c_mn$basesCovered * c_mn$mean)) * 1e6
        c_m <- metaprofile.bigWig(bed, c_mn, c_pl, step=stp)$middle/ stp/ (c_pl$basesCovered * c_pl$mean + abs(c_mn$basesCovered * c_mn$mean)) * 1e6
        e_p <- metaprofile.bigWig(bed, e_pl, e_mn, step=stp)$middle/ stp / (e_pl$basesCovered * e_pl$mean + abs(e_mn$basesCovered * e_mn$mean)) * 1e6
        e_m <- metaprofile.bigWig(bed, e_mn, e_pl, step=stp)$middle/ stp/ (e_pl$basesCovered * e_pl$mean + abs(e_mn$basesCovered * e_mn$mean)) * 1e6
        e_b_p <- metaprofile.bigWig(bed, e_b_pl, e_b_mn, step=stp)$middle/ stp / (e_b_pl$basesCovered * e_b_pl$mean + abs(e_b_mn$basesCovered * e_b_mn$mean)) * 1e6
        e_b_m <- metaprofile.bigWig(bed, e_b_mn, e_b_pl, step=stp)$middle/ stp/ (e_b_pl$basesCovered * e_b_pl$mean + abs(e_b_mn$basesCovered * e_b_mn$mean)) * 1e6

        N = length(c_p)
        x = ((1:N) - N/2)* stp #1:N*stp
        ylim=c(-1*max(c(c_m, e_m, e_b_m)), max(c_p, e_p, e_b_p))

        plot(-500, -500, ylim=ylim, xlim=c(min(x), max(x)), xlab= "Distance [bp]", ylab= "Signal [RPKM].", ...)
        lines(x, c_p, col="#C14F4D")
        lines(x, (-1*c_m), col="#C14F4D")
        lines(x, e_p, col="dark red")
        lines(x, (-1*e_m), col="dark red")
        lines(x, e_b_p, col="#7E75B1")
        lines(x, (-1*e_b_m), col="#7E75B1")

}

erbs <- read.table("ERaBS.bed")
pdf("E2_trt_BBCA.pdf")
doit(erbs, stp=20, halfWindow=2000, B7_c_pl, B7_c_mn, B7_e2_pl, B7_e2_mn, B7_e2b_pl, B7_e2b_mn, main="B7 E2 stimulation")
doit(erbs, stp=20, halfWindow=2000, G11_c_pl, G11_c_mn, G11_e2_pl, G11_e2_mn, G11_e2b_pl, G11_e2b_mn, main="G11 E2 stimulation")
dev.off()
