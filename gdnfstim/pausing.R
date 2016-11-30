## Quantile normalize and get data!

## Libraries
library(preprocessCore)

## Load data, remove 0s, quantile normalize.
load("Counts.RData")
gbd <- gene_body_density
tsd <- tss_density
keep <- rowSums(gbd)>0 & rowSums(tsd)>0
gbd  <- normalize.quantiles(gbd[keep,])
tsd  <- normalize.quantiles(tsd[keep,])

## Compute pausing indices.
PI_B7_0h <- tsd[,1]/ gbd[,1]
PI_B7_1h <- tsd[,2]/ gbd[,2]
PI_C11_0h<- tsd[,4]/ gbd[,4]
PI_C11_1h<- tsd[,5]/ gbd[,5]

cngPI_B7 <- log(PI_B7_1h)-log(PI_B7_0h)
cngPI_C11<- log(PI_C11_1h)-log(PI_C11_0h)
cngPI_TamS <- rowMeans(cbind(cngPI_B7, cngPI_C11))

## Get pvalues.  Get up/ down regulated genes.
pv <- read.table("GDNF_treatment.tsv", header=TRUE)
sum(pv$GENEID == as.character(bodies$GENEID))/NROW(bodies) ## Sanity check!  MUST be 1.
pv <- pv[keep,]

upreg <- !is.na(pv$FDR_GDNF_1h_S) & pv$FDR_GDNF_1h_S < 0.01 & pv$FC_GDNF_1h_S > 0
dnreg <- !is.na(pv$FDR_GDNF_1h_S) & pv$FDR_GDNF_1h_S < 0.01 & pv$FC_GDNF_1h_S < 0
noreg <- !is.na(pv$FDR_GDNF_1h_S) & pv$FDR_GDNF_1h_S > 0.5  & abs(pv$FC_GDNF_1h_S) < 0.1

## Compare changes in PI between groups.
boxplot(cngPI_B7[noreg], cngPI_C11[noreg], cngPI_B7[upreg], cngPI_C11[upreg], cngPI_B7[dnreg], cngPI_C11[dnreg])
boxplot((tsd[noreg,2]-tsd[noreg,1]), (tsd[upreg,2]-tsd[upreg,1]), (tsd[dnreg,2]-tsd[dnreg,1]), outline=FALSE)
boxplot((gbd[noreg,2]-gbd[noreg,1]), (gbd[upreg,2]-gbd[upreg,1]), (gbd[dnreg,2]-gbd[dnreg,1]), outline=FALSE)

pdf(file="PausingIndex.pdf", width=4, height=6)
  labels <- paste(c("No change", "Up-regulated", "Down-regulated"))
  boxplot(cngPI_TamS[noreg], cngPI_TamS[upreg], cngPI_TamS[dnreg], ylab="Change in pausing index", xlab="", xaxt="n", outline=FALSE)
  abline(h=0)
#  axis(1, labels=FALSE)
  text(x= seq_along(labels), y=par("usr")[3]-0.1, srt=35, adj=1, labels= labels, xpd= TRUE)
dev.off()

wilcox.test(cngPI_TamS[noreg], cngPI_TamS[upreg])
wilcox.test(cngPI_TamS[noreg], cngPI_TamS[upreg])$p.value

wilcox.test(cngPI_TamS[noreg], cngPI_TamS[dnreg])
wilcox.test(cngPI_TamS[noreg], cngPI_TamS[dnreg])$p.value

