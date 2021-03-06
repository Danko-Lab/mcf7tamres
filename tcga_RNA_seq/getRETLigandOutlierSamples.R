#!/usr/bin/bash

esr1 <- read.table("ESR1.tmp")
donorID <- unique(esr1$V4) ## Use sampleID

read.dataset<- function(MARK, IDs) {
  data <- read.table(paste(MARK, ".tmp", sep=""))
  return(data[match(IDs, data$V4),])
}

esr1 <- read.dataset("ESR1", donorID)
ret  <- read.dataset("RET",  donorID)
egr1 <- read.dataset("EGR1", donorID)
gdnf <- read.dataset("GDNF", donorID)
nrtn <- read.dataset("NRTN", donorID)
artn <- read.dataset("ARTN", donorID)
pspn <- read.dataset("PSPN", donorID)

gfra1<- read.dataset("GFRA1", donorID)
gfra2<- read.dataset("GFRA2", donorID)
gfra3<- read.dataset("GFRA3", donorID)
gfra4<- read.dataset("GFRA4", donorID)

retL <- data.frame(GDNF= gdnf$V9, NRTN= nrtn$V9, ARTN= artn$V9, PSPN= pspn$V9, RET= ret$V9, EGR1= egr1$V9, ESR1= esr1$V9, GFRA1= gfra1$V9, GFRA2= gfra2$V9, GFRA3= gfra3$V9, GFRA4= gfra4$V9)
rownames(retL) <- donorID
cor(retL)

## Find ER+ samples.
ERpos <- -5

pc <- 2.208627e-09

pdf("ERpos.pdf")
  hist(log(retL$ESR1,10), 20, main="ESR1", xlab="ESR1 expression"); abline(v=-5, lty="dashed", col="blue")
  sum(log(retL$ESR1,10) > ERpos)/NROW(retL)
dev.off()

retL <- retL[log(retL$ESR1,10)> ERpos,] ## REMOVE ER- patients

## Now standardize.
retL.med <- sapply(1:NCOL(retL), function(i) {median(log(retL[retL[,i]>0,i], 10), na.rm=TRUE)} )
retL.var <- sapply(1:NCOL(retL), function(i) {var(log(retL[retL[,i]>0,i], 10), na.rm=TRUE)} )
retL.IQR <- sapply(1:NCOL(retL), function(i) {IQR(log(retL[retL[,i]>0,i], 10), na.rm=TRUE)} )
retL.2550 <- 2.5*sapply(1:NCOL(retL), function(i) {quantile(log(retL[retL[,i]>0,i], 10), 0.5, na.rm=TRUE)-quantile(log(retL[retL[,i]>0,i], 10), 0.25, na.rm=TRUE)}  )
retL.050 <- sapply(1:NCOL(retL), function(i) {quantile(log(retL[retL[,i]>0,i], 10), 0.5, na.rm=TRUE)-quantile(log(retL[retL[,i]>0,i], 10), 0, na.rm=TRUE)} )

retL.stdZ <- scale(log(retL,10), center=retL.med, scale=scale=sqrt(retL.var))
TH= 1.28 ## >90% of data, assuming standard normal. ##1.645 ## >0.95% of data, assuming standard normal.

#CUTOFF_VALS <- TH*sqrt(retL.var) + retL.med ## Assume standard normal

CUTOFF_VALS <- retL.med + retL.2550 ## median+IQR

#CUTOFF_VALS <- retL.med + retL.050

pdf("Ret.Ligand.Expression.Distribution.pdf")

 par(mfrow=c(2,2))
 hist(log(retL$GDNF, 10), 20, main="GDNF", xlab="GDNF Expression"); abline(v=CUTOFF_VALS[1], lty="dashed", col="blue")
 hist(log(retL$NRTN, 10), 20, main="NRTN", xlab="NRTN Expression"); abline(v=CUTOFF_VALS[2], lty="dashed", col="blue")
 hist(log(retL$ARTN, 10), 20, main="ARTN", xlab="ARTN Expression"); abline(v=CUTOFF_VALS[3], lty="dashed", col="blue")
 hist(log(retL$PSPN, 10), 20, main="PSPN", xlab="PSPN Expression"); abline(v=CUTOFF_VALS[4], lty="dashed", col="blue")

# par(mfrow=c(2,2))
# hist(retL.stdZ[,"GDNF"], 20, main="GDNF", xlab="GDNF Expression"); abline(v=TH, lty="dashed", col="blue")
# hist(retL.stdZ[,"NRTN"], 20, main="NRTN", xlab="NRTN Expression"); abline(v=TH, lty="dashed", col="blue")
# hist(retL.stdZ[,"ARTN"], 20, main="ARTN", xlab="ARTN Expression"); abline(v=TH, lty="dashed", col="blue")
# hist(retL.stdZ[,"PSPN"], 20, main="PSPN", xlab="PSPN Expression"); abline(v=TH, lty="dashed", col="blue")

dev.off()

sum(log(retL$GDNF,10) > CUTOFF_VALS[1] & log(retL$NRTN,10) < CUTOFF_VALS[2] &  log(retL$ARTN,10) < CUTOFF_VALS[3] & log(retL$PSPN,10) < CUTOFF_VALS[4])
sum(log(retL$NRTN,10) > CUTOFF_VALS[2] & log(retL$GDNF,10) < CUTOFF_VALS[1] &  log(retL$ARTN,10) < CUTOFF_VALS[3] & log(retL$PSPN,10) < CUTOFF_VALS[4])
sum(log(retL$ARTN,10) > CUTOFF_VALS[3] & log(retL$NRTN,10) < CUTOFF_VALS[2] &  log(retL$GDNF,10) < CUTOFF_VALS[1] & log(retL$PSPN,10) < CUTOFF_VALS[4])
sum(log(retL$PSPN,10) > CUTOFF_VALS[4] & log(retL$NRTN,10) < CUTOFF_VALS[2] &  log(retL$ARTN,10) < CUTOFF_VALS[3] & log(retL$GDNF,10) < CUTOFF_VALS[1])

any <- sum(log(retL$GDNF,10) > CUTOFF_VALS[1] | log(retL$NRTN,10) > CUTOFF_VALS[2] | log(retL$ARTN,10) > CUTOFF_VALS[3] | log(retL$PSPN,10) > CUTOFF_VALS[4]); any

multiple <- data.frame(GDNF= (log(retL$GDNF,10) > CUTOFF_VALS[1]), NRTN= (log(retL$NRTN,10) > CUTOFF_VALS[2]), ARTN= (log(retL$ARTN,10) > CUTOFF_VALS[3]), PSPN= (log(retL$PSPN,10) > CUTOFF_VALS[4]))
sum(rowSums(multiple) > 1)

sum(log(retL$GDNF,10) > CUTOFF_VALS[1] | log(retL$NRTN,10) > CUTOFF_VALS[2] | log(retL$ARTN,10) > CUTOFF_VALS[3] | log(retL$PSPN,10) > CUTOFF_VALS[4])/ NROW(retL)

## Make pie chart.


## Make scatter plots.
cor(retL)

getMod <- function(val1, val2) {
  glm(log(val1[val1>0 & val2>0],10)~log(val2[val1>0 & val2>0],10))
} 

pdf("ESR1.RET.coR.DensityScatterplot.pdf", height=5, width=5)
 source("../lib/densScatterplot.R")
 densScatterplot(retL$ESR1, retL$RET, uselog=TRUE, xlab="ESR1", ylab="RET");  mod1 <- getMod(retL$RET, retL$ESR1); abline(mod1)
 densScatterplot(retL$ESR1, retL$GDNF, uselog=TRUE, xlab="ESR1", ylab="GDNF"); #mod1 <- getMod(retL$GDNF, retL$ESR1); abline(mod1)
 densScatterplot(retL$ESR1, retL$NRTN, uselog=TRUE, xlab="ESR1", ylab="NRTN"); mod1 <- getMod(retL$NRTN, retL$ESR1); abline(mod1)
 densScatterplot(retL$ESR1, retL$ARTN, uselog=TRUE, xlab="ESR1", ylab="ARTN"); #mod1 <- getMod(retL$ARTN, retL$ESR1); abline(mod1)
 densScatterplot(retL$ESR1, retL$PSPN, uselog=TRUE, xlab="ESR1", ylab="PSPN"); #mod1 <- getMod(retL$PSPN, retL$ESR1); abline(mod1)

 densScatterplot(retL$ESR1, retL$GFRA1, uselog=TRUE, xlab="ESR1", ylab="GFRA1");  mod1 <- getMod(retL$GFRA1, retL$ESR1); abline(mod1)
 densScatterplot(retL$ESR1, retL$GFRA2, uselog=TRUE, xlab="ESR1", ylab="GFRA2")
 densScatterplot(retL$ESR1, retL$GFRA3, uselog=TRUE, xlab="ESR1", ylab="GFRA3")
 densScatterplot(retL$ESR1, retL$GFRA4, uselog=TRUE, xlab="ESR1", ylab="GFRA4")

dev.off()

pdf("ESR1.EGR1.pdf", height=5, width=5)
 densScatterplot(retL$ESR1, retL$EGR1, uselog=TRUE, xlab="ESR1", ylab="EGR1")
# densScatterplot(retL$ESR1, retL$EGR1, uselog=TRUE, xlab="ESR1", ylab="EGR1");mod1 <- glm(log(egr1$V9[indx6][ERpos],10)~log(esr1$V9[ERpos],10));abline(mod1)
dev.off()

## Now make a violin plot.
require(vioplot)
pdf("S4B.RETLexpr.TCGA.pdf", height=4, width=8)
vioplot(log(retL$RET+pc, 10), 
	log(retL$GFRA1+pc, 10), 
        log(retL$GDNF+pc, 10),         

	log(retL$GFRA2+pc, 10), 
        log(retL$NRTN+pc, 10),         

	log(retL$GFRA3+pc, 10), 
        log(retL$ARTN+pc, 10),

	log(retL$GFRA4+pc, 10), 
        log(retL$PSPN+pc, 10),

	names=c("RET", "GFRA1", "GDNF", "GFRA2", "NRTN", "GFRA3", "ARTN", "GFRA4", "PSPN"), 
	col=c("dark green"))
dev.off()

