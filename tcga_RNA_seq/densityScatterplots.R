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

ERpos <- esr1$V9 > 1e-5

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

## Variation on densScatterplot that turns ER- genes gray.
densScatterplot.GrayERneg <- function(x1, x2, uselog=FALSE, n=256, ...) {
  df <- data.frame(x1, x2)

  if(uselog) {
    x <- densCols(log(x1,10),log(x2,10), colramp=colorRampPalette(c("black", "white")))
  } else {
    x <- densCols(x1,x2, colramp=colorRampPalette(c("black", "white")))
  }
  df$dens <- col2rgb(x)[1,] + 1L

  ## Map densities to colors
#  cols <-  colorRampPalette(c("#000099", "#00FEFF", "#45FE4F", "#FCFF00", "#FF9400", "#FF3100"))(256)
  cols <- colorRampPalette(c("dark gray", "#000099", "#45FE4F", "#FCFF00", "#FF9400", "#FF3100"))(n)
  df$col <- cols[df$dens]
  df$col[!ERpos] <- "light gray"

  ## Plot it, reordering rows so that densest points are plotted on top
  if(uselog) {
    plot(x2~x1, data=df[order(df$dens),], pch=20, col=col, cex=2, log="xy", ...)
  } else {
    plot(x2~x1, data=df[order(df$dens),], pch=20, col=col, cex=2, ...)
  }
}

##require(groHMM) # In groHMM, but not exported.
tlsDeming <- function(x,y,d=1) {
 sxx <- 1/(NROW(x)-1) * sum((x-mean(x))^2)
 sxy <- 1/(NROW(x)-1) * sum((x-mean(x))*(y-mean(y)))
 syy <- 1/(NROW(y)-1) * sum((y-mean(y))^2)
 X <- (syy-d*sxx+sqrt((syy-d*sxx)^2+4*d*(sxy^2)))/(2*sxy)
 int <- mean(y) - X * mean(x)
 return(c(int,X))
}

pdf("ESR1.EGR1.pdf", height=5, width=5)
 densScatterplot.GrayERneg(esr1$V9, egr1$V9[indx6], uselog=TRUE, xlab="ESR1", ylab="EGR1"); abline(v=1e-5); 
# densScatterplot.GrayERneg(log(esr1$V9, 10), log(egr1$V9[indx6], 10), uselog=FALSE, xlab="ESR1", ylab="EGR1"); abline(v=-5); 
 mod1 <- tlsDeming(y= log(egr1$V9[indx6][ERpos],10), x= log(esr1$V9[ERpos],10));abline(mod1)

 densScatterplot(esr1$V9[ERpos], egr1$V9[indx6][ERpos], uselog=TRUE, xlab="ESR1", ylab="EGR1");mod1 <- glm(log(egr1$V9[indx6][ERpos],10)~log(esr1$V9[ERpos],10));abline(mod1)
dev.off()

 cor.test(gdnf$V9[indx2], egr1$V9[indx6], method="spearman")
 cor.test(nrtn$V9[indx3], egr1$V9[indx6], method="spearman")
 cor.test(artn$V9[indx4], egr1$V9[indx6], method="spearman")
 cor.test(msmb$V9[indx5], egr1$V9[indx6], method="spearman")

 densScatterplot(artn$V9[indx4], egr1$V9[indx6], uselog=TRUE, xlab="ARTN", ylab="EGR1")

