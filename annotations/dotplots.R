## Plot out RPKM normalized data from TamR/ TamS lines.

refGene <- read.table("tuSelecter/final_tus.txt", header=TRUE)
refGene <- rbind(refGene, read.table("tuSelecter/final_tus.ESR1_GREB1.txt", header=TRUE))

refGene <- refGene[grep("random|Un|hap", refGene$TXCHROM, invert=TRUE),]
refGene <- refGene[(refGene$TXEND-refGene$TXSTART)>4000,]

bodies <- refGene
bodies$TXSTART[bodies$TXSTRAND == "+"] <-bodies$TXSTART[bodies$TXSTRAND == "+"]+1000
bodies$TXEND[bodies$TXSTRAND == "-"] <- bodies$TXEND[bodies$TXSTRAND == "-"]-1000

tss <- refGene
tss$TXEND[tss$TXSTRAND == "+"] <-tss$TXSTART[tss$TXSTRAND == "+"]+1000
tss$TXSTART[tss$TXSTRAND == "-"] <- tss$TXEND[tss$TXSTRAND == "-"]-1000

## Read in bigWigs.
require(bigWig)

getCounts <- function(bwplus_name, bwminus_name, path="../data/", gb= bodies, rpkm=FALSE) {
	bwplus <- load.bigWig(paste(path, bwplus_name, sep=""))
	bwminus<- load.bigWig(paste(path, bwminus_name,sep=""))

        counts <- bed6.region.bpQuery.bigWig(bwplus, bwminus, gb, abs.value = TRUE)
        if(rpkm==TRUE) {
                counts <- counts * (1000/(gb[,3]-gb[,2])) * (1e6/(abs(bwplus$mean)*bwplus$basesCovered+abs(bwminus$mean)*bwminus$basesCovered))
        }
        return(counts)
}

B7 <- getCounts("sB7_pl.bw", "sB7_mn.bw", rpkm=TRUE)
G11<- getCounts("rG11_pl.bw", "rG11_mn.bw", rpkm=TRUE)
H9 <- getCounts("rH9_pl.bw", "rH9_mn.bw", rpkm=TRUE)
C11<- getCounts("sC11_pl.bw", "sC11_mn.bw", rpkm=TRUE)

rpkm_norm <- cbind(B7, C11, H9, G11)

pdf("ret.circplot.pdf")
source("../lib/circplot.R")
cd.circplot(log(rpkm_norm[bodies$GENENAME == "RET",],2), c("TamS", "TamS", "TamR", "TamR"), lims=c(3,-3))
dev.off()


## Get the first 60 kb of gene body.
bodies60 <- refGene
bodies60$TXSTART[bodies60$TXSTRAND == "+"] <-bodies60$TXSTART[bodies60$TXSTRAND == "+"]+1000
bodies60$TXEND[bodies60$TXSTRAND == "-"] <- bodies60$TXEND[bodies60$TXSTRAND == "-"]-1000
size <- bodies60$TXEND-bodies60$TXSTART
maxsize <- 60000
bodies60$TXEND[bodies60$TXSTRAND == "+" & size>maxsize] <- refGene$TXSTART[bodies60$TXSTRAND == "+" & size>maxsize]+maxsize
bodies60$TXSTART[bodies60$TXSTRAND == "-" & size>maxsize] <- refGene$TXEND[bodies60$TXSTRAND == "-" & size>maxsize]-maxsize

gdnf_rpkm_norm <- data.frame(
	B7_0h= getCounts("MCF-7_B7_GDNF_0_hr_plus.bw", "MCF-7_B7_GDNF_0_hr_minus.bw", gb=bodies60, rpkm=TRUE),
	C11_0h=getCounts("MCF-7_C11_GDNF_0_hr_plus.bw", "MCF-7_C11_GDNF_0_hr_minus.bw", gb=bodies60, rpkm=TRUE),
	G11_0h=getCounts("MCF-7_G11_GDNF_0_hr_plus.bw", "MCF-7_G11_GDNF_0_hr_minus.bw", gb=bodies60, rpkm=TRUE),
	H9_0h= getCounts("MCF-7_H9_GDNF_0_hr_plus.bw", "MCF-7_H9_GDNF_0_hr_minus.bw", gb=bodies60, rpkm=TRUE),

        B7_1h= getCounts("MCF-7_B7_GDNF_1_hr_plus.bw", "MCF-7_B7_GDNF_1_hr_minus.bw", gb=bodies60, rpkm=TRUE),
        C11_1h=getCounts("MCF-7_C11_GDNF_1_hr_plus.bw", "MCF-7_C11_GDNF_1_hr_minus.bw", gb=bodies60, rpkm=TRUE),
        G11_1h=getCounts("MCF-7_G11_GDNF_1_hr_plus.bw", "MCF-7_G11_GDNF_1_hr_minus.bw", gb=bodies60, rpkm=TRUE),
        H9_1h= getCounts("MCF-7_H9_GDNF_1_hr_plus.bw", "MCF-7_H9_GDNF_1_hr_minus.bw", gb=bodies60, rpkm=TRUE),

        B7_24h= getCounts("MCF-7_B7_GDNF_24_hr_plus.bw", "MCF-7_B7_GDNF_24_hr_minus.bw", gb=bodies60, rpkm=TRUE),
        C11_24h=getCounts("MCF-7_C11_GDNF_24_hr_plus.bw", "MCF-7_C11_GDNF_24_hr_minus.bw", gb=bodies60, rpkm=TRUE),
        G11_24h=getCounts("MCF-7_G11_GDNF_24_hr_plus.bw", "MCF-7_G11_GDNF_24_hr_minus.bw", gb=bodies60, rpkm=TRUE),
        H9_24h= getCounts("MCF-7_H9_GDNF_24_hr_plus.bw", "MCF-7_H9_GDNF_24_hr_minus.bw", gb=bodies60, rpkm=TRUE))

gdnf_tamsr <- rep(c("TamS", "TamS", "TamR", "TamR"),3)
gdnf_trt   <- c(rep("0h",4), rep("1h",4), rep("24h",4))
gdnf_clone <- rep(c("B7", "C11", "G11", "H9"), 3)

pdf("gdnf_stim.circplot.pdf")#, height=5, width=5)

cd.circplot(t(log(gdnf_rpkm_norm[bodies$GENENAME == "ESR1.iso2", gdnf_tamsr=="TamS"],2)), gdnf_trt[gdnf_tamsr=="TamS"])
cd.circplot(t(log(gdnf_rpkm_norm[bodies$GENENAME == "ESR1.iso2", gdnf_tamsr=="TamR"],2)), gdnf_trt[gdnf_tamsr=="TamR"])
#cd.circplot(t(log(gdnf_rpkm_norm[bodies$GENENAME == "ESR1", gdnf_tamsr=="TamS"],2)), gdnf_trt[gdnf_tamsr=="TamS"])
#cd.circplot(t(log(gdnf_rpkm_norm[bodies$GENENAME == "ESR1", gdnf_tamsr=="TamR"],2)), gdnf_trt[gdnf_tamsr=="TamR"])
#cd.circplot(t(log(gdnf_rpkm_norm[bodies$GENENAME == "ESR1", gdnf_clone=="C11"],2)), gdnf_trt[gdnf_clone=="C11"])

cd.circplot(t(log(gdnf_rpkm_norm[bodies$GENENAME == "GREB1", gdnf_tamsr=="TamS"],2)), gdnf_trt[gdnf_tamsr=="TamS"])
cd.circplot(t(log(gdnf_rpkm_norm[bodies$GENENAME == "PGR", gdnf_tamsr=="TamS"],2)), gdnf_trt[gdnf_tamsr=="TamS"])
cd.circplot(t(log(gdnf_rpkm_norm[bodies$GENENAME == "ELOVL2", gdnf_tamsr=="TamS"],2)), gdnf_trt[gdnf_tamsr=="TamS"])
cd.circplot(t(log(gdnf_rpkm_norm[bodies$GENENAME == "NOS1AP", gdnf_tamsr=="TamS"],2)), gdnf_trt[gdnf_tamsr=="TamS"])

dev.off()

#

