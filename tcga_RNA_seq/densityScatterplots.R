## Script to create scatterplots.

esr1 <- read.table("ESR1.tmp")
ret <- read.table("RET.tmp")
gdnf <- read.table("GDNF.tmp")
nrtn <- read.table("NRTN.tmp")
artn <- read.table("ARTN.tmp")
msmb <- read.table("PSPN.tmp")
egr1 <- read.table("EGR1.tmp")

indx <- match(esr1$V1, ret$V1)
indx2<-match(esr1$V1, gdnf$V1)
indx3<-match(esr1$V1, nrtn$V1)
indx4<-match(esr1$V1, artn$V1)
indx5<-match(esr1$V1, msmb$V1)
indx6<-match(esr1$V1, egr1$V1)

sum(esr1$V1 == ret$V1)
sum(esr1$V1 == ret$V1[indx])
sum(esr1$V1 == gdnf$V1[indx2])
sum(esr1$V1 == egr1$V1[indx6])

ERpos <- esr1$V9 > 5e-5

sum(esr1$V1[ERpos] == egr1$V1[indx6][ERpos])

cor.test(esr1$V9, ret$V9[indx], method="spearman") #[Rho= 0.4500471; p= 1.195204e-60]
cor.test(esr1$V9, gdnf$V9[indx2], method="spearman") #[Rho= -0.0358615; p= 0.2154]
cor.test(esr1$V9, nrtn$V9[indx3], method="spearman") #
cor.test(esr1$V9, artn$V9[indx4], method="spearman") #
cor.test(esr1$V9, msmb$V9[indx5], method="spearman") #
cor.test(esr1$V9, egr1$V9[indx6], method="spearman") #

cor.test(esr1$V9[ERpos], ret$V9[indx][ERpos], method="spearman") # [Rho= 0.1232322; p= 0.0005438]
cor.test(esr1$V9[ERpos] , gdnf$V9[indx2][ERpos] , method="spearman") #
cor.test(esr1$V9[ERpos] , nrtn$V9[indx3][ERpos], method="spearman") #
cor.test(esr1$V9[ERpos] , artn$V9[indx4][ERpos], method="spearman") #
cor.test(esr1$V9[ERpos] , msmb$V9[indx5][ERpos], method="spearman") #
cor.test(esr1$V9[ERpos] , egr1$V9[indx6][ERpos], method="spearman") # [Rho= -0.2191036; p= 5.625e-10]

pdf("ESR1.RET.DensityScatterplot.pdf", height=5, width=5)
 source("../lib/densScatterplot.R")
 densScatterplot(esr1$V9, ret$V9[indx], uselog=TRUE, xlab="ESR1", ylab="RET"); mod1 <- glm(log(ret$V9[indx][ERpos],10)~log(esr1$V9[ERpos],10));abline(mod1)
 densScatterplot(esr1$V9, gdnf$V9[indx2], uselog=TRUE, xlab="ESR1", ylab="GDNF")
 densScatterplot(esr1$V9, nrtn$V9[indx3], uselog=TRUE, xlab="ESR1", ylab="NRTN")
 densScatterplot(esr1$V9, artn$V9[indx4], uselog=TRUE, xlab="ESR1", ylab="ARTN")
 densScatterplot(esr1$V9, msmb$V9[indx5], uselog=TRUE, xlab="ESR1", ylab="PSPN")
dev.off()

pdf("ESR1.EGR1.pdf", height=5, width=5)
 densScatterplot(esr1$V9, egr1$V9[indx6], uselog=TRUE, xlab="ESR1", ylab="EGR1")
 densScatterplot(esr1$V9[ERpos], egr1$V9[indx6][ERpos], uselog=TRUE, xlab="ESR1", ylab="EGR1");mod1 <- glm(log(egr1$V9[indx6][ERpos],10)~log(esr1$V9[ERpos],10));abline(mod1)
dev.off()

 cor.test(gdnf$V9[indx2], egr1$V9[indx6], method="spearman")
 cor.test(nrtn$V9[indx3], egr1$V9[indx6], method="spearman")
 cor.test(artn$V9[indx4], egr1$V9[indx6], method="spearman")
 cor.test(msmb$V9[indx5], egr1$V9[indx6], method="spearman")

 densScatterplot(artn$V9[indx4], egr1$V9[indx6], uselog=TRUE, xlab="ARTN", ylab="EGR1")

