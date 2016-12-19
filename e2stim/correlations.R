## get correlations between MCF-7 samples

## Read in refseq genes.
refGene <- read.table("../annotations/refGene.bed.gz")
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

GC <- function(prefix, path="../data/", rpkm=TRUE) {
 pl <- load.bigWig(paste(path, prefix, "_plus.bw", sep=""))
 mn <- load.bigWig(paste(path, prefix, "_minus.bw", sep=""))

 ## Count reads in each ...
 cnts <- bed6.region.bpQuery.bigWig(pl, mn, bodies, abs.value = TRUE) #/(pl$)
        if(rpkm==TRUE) {
                cnts <- cnts * (1000/(bodies[,3]-bodies[,2])) * (1e6/(abs(pl$mean)*pl$basesCovered+abs(mn$mean)*mn$basesCovered))
        }
 cnts
}

rpkm <- cbind(B7=GC("B7_CTRL"),
                          B7.e2=GC("B7_E2"),
                          B7.tam=GC("B7_TAM"),
                          B7.Be2=GC("B7_E2_BBCA"),
                          B7.Btam=GC("B7_TAM_BBCA"),

                          G11=GC("G11_CTRL"),
                          G11.e2=GC("G11_E2"),
                          G11.tam=GC("G11_TAM"),
                          G11.Be2=GC("G11_E2_BBCA"),
                          G11.Btam=GC("G11_TAM_BBCA"))

## Now get the correlations between tam and untreated.
cor.test(rpkm[,"B7"], rpkm[,"B7.tam"], method="spearman")
cor.test(rpkm[,"G11"], rpkm[,"G11.tam"], method="spearman")

#cor.test(log(rpkm[,"B7"]), log(rpkm[,"B7.tam"]), method="pearson")
#cor.test(log(rpkm[,"G11"]), log(rpkm[,"G11.tam"]), method="pearson")

## Draw density scatterplots.
source("~/NHP/lib/densScatterplot.R")

pdf("S2B.UNT-TAM.densityScatterplots.pdf")

densScatterplot(rpkm[,"B7"], rpkm[,"B7.tam"], uselog=TRUE, xlab="PRO-seq, Untreated [RPKM]", ylab="PRO-seq, Tam [rpkm]", main="B7")
densScatterplot(rpkm[,"G11"], rpkm[,"G11.tam"], uselog=TRUE, xlab="PRO-seq, Untreated [RPKM]", ylab="PRO-seq, Tam [rpkm]", main="G11")

dev.off()

densScatterplot(rpkm[,"B7"], rpkm[,"B7.e2"], uselog=TRUE, xlab="PRO-seq, Untreated [RPKM]", ylab="PRO-seq, E2 [rpkm]", main="B7")
densScatterplot(rpkm[,"G11"], rpkm[,"G11.e2"], uselog=TRUE, xlab="PRO-seq, Untreated [RPKM]", ylab="PRO-seq, E2 [rpkm]", main="G11")



