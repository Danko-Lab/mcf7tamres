## Q: Do Letrozole non-responders have a higher amount of RET or RET lignads??
## A: RET lignads!
#
# Data from: GDS3116; GSE20181; 
# Paper: http://www.nature.com.proxy.library.cornell.edu/tpj/journal/v12/n1/full/tpj201067a.html 
#
#GDNF: 221359_at
#NRTN: 210683_at
#ARTN: 210237_at
#PSPN: 221373_x_at
#RET: 211421_s_at

rl <- read.table("RET_ligands.tsv", sep="\t", header=TRUE)
r  <- read.table("RET.tsv", sep="\t", header=TRUE)

## Gain power by averaging patients.
df <- data.frame(Patient= unique(rl$Patient))
df$resp <- rl$Response[match(df$Patient, rl$Patient)]
df$gdnf <- sapply(df$Patient, function(x) {mean(rl$GDNF[rl$Patient == x])})
df$nrtn <- sapply(df$Patient, function(x) {mean(rl$NRTN[rl$Patient == x])})
df$artn <- sapply(df$Patient, function(x) {mean(rl$ARTN[rl$Patient == x])})
df$pspn <- sapply(df$Patient, function(x) {mean(rl$PSPN[rl$Patient == x])})
df$ret  <- sapply(df$Patient, function(x) {mean(r$RET[rl$Patient == x])})
#df$egr1 <- sapply(df$Patient, function(x) {mean(egr$EGR1[egr$Patient == x])})

IQRls <- function(x) {(quantile(x, 0.5)-quantile(x, 0))}

df$gdnf <- scale(df$gdnf, center=median(df$gdnf), scale=IQRls(df$gdnf))
df$nrtn <- scale(df$nrtn, center=median(df$nrtn), scale=IQRls(df$nrtn))
df$artn <- scale(df$artn, center=median(df$artn), scale=IQRls(df$artn))
df$pspn <- scale(df$pspn, center=median(df$pspn), scale=IQRls(df$pspn))

datacols <- 3:6 #c(3:6) #3:6

wilcox.test(rowSums(df[df$resp==" responder",datacols]), rowSums(df[df$resp==" nonresponder",datacols]), alternative= "less")
wilcox.test(df[df$resp==" responder","ret"], df[df$resp==" nonresponder","ret"])

pdf("RET-Ligand-Score.pdf", height=6, width=4)

	boxplot(rowSums(df[df$resp==" responder",datacols]), rowSums(df[df$resp==" nonresponder",datacols]), 
		names=c("Responder", "Non-Responder"), outline=FALSE, ylim=c(-3, 12.5), ylab="RET Ligands Score", main="RET ligands")
	stripchart(list(rowSums(df[df$resp==" responder",datacols]), rowSums(df[df$resp==" nonresponder",datacols])), 
		vertical=TRUE, method="jitter", add=TRUE, pch=20, col="blue", cex=2)

	boxplot(df$ret[df$resp==" responder"], df$ret[df$resp==" nonresponder"], 
		names=c("Responder", "Non-Responder"), outline=FALSE, ylim=c(0, 610), ylab="RET expression level", main="RET")
	stripchart(list(df[df$resp==" responder","ret"], df[df$resp==" nonresponder","ret"]), 
		vertical=TRUE, method="jitter", add=TRUE, pch=20, col="blue", cex=2)

dev.off()

## Plot separately by time point.
pdf("RET-Ligand-RET.pdf", height=6, width=12)

df <- data.frame(Patient= rl$Patient, time= rl$Time, resp= rl$Response, gdnf= rl$GDNF, nrtn=rl$NRTN, artn=rl$ARTN, pspn=rl$PSPN)

df$gdnf <- scale(df$gdnf, center=median(df$gdnf), scale=IQRls(df$gdnf))
df$nrtn <- scale(df$nrtn, center=median(df$nrtn), scale=IQRls(df$nrtn))
df$artn <- scale(df$artn, center=median(df$artn), scale=IQRls(df$artn))
df$pspn <- scale(df$pspn, center=median(df$pspn), scale=IQRls(df$pspn))

datacols <- 4:7

wilcox.test(rowSums(df[df$time=="t90" & df$resp==" responder",datacols]), rowSums(df[df$time=="t90" & df$resp==" nonresponder",datacols]))
wilcox.test(rowSums(df[df$time=="t10" & df$resp==" responder",datacols]), rowSums(df[df$time=="t10" & df$resp==" nonresponder",datacols]))
wilcox.test(rowSums(df[df$time=="t0" & df$resp==" responder",datacols]), rowSums(df[df$time=="t0" & df$resp==" nonresponder",datacols]))

data <- list(rowSums(df[df$time=="t0" & df$resp==" responder",datacols]), rowSums(df[df$time=="t0" & df$resp==" nonresponder",datacols]),
	rowSums(df[df$time=="t10" & df$resp==" responder",datacols]), rowSums(df[df$time=="t10" & df$resp==" nonresponder",datacols]),
	rowSums(df[df$time=="t90" & df$resp==" responder",datacols]), rowSums(df[df$time=="t90" & df$resp==" nonresponder",datacols]))

boxplot(data, names=c("0 days, respond", "0 days, non-respond", "10 days, respond", "10 days, non-respond", "90 days, respond", "90 days, non-respond" ))

## Then do RET.
wilcox.test(r$RET[r$Time=="t0" & r$Response==" responder"], r$RET[r$Time=="t0" & r$Response==" nonresponder"])
wilcox.test(r$RET[r$Time=="t10" & r$Response==" responder"], r$RET[r$Time=="t10" & r$Response==" nonresponder"])
wilcox.test(r$RET[r$Time=="t90" & r$Response==" responder"], r$RET[r$Time=="t90" & r$Response==" nonresponder"])

data <- list(r$RET[r$Time=="t0" & r$Response==" responder"], r$RET[r$Time=="t0" & r$Response==" nonresponder"],
	r$RET[r$Time=="t10" & r$Response==" responder"], r$RET[r$Time=="t10" & r$Response==" nonresponder"],
	r$RET[r$Time=="t90" & r$Response==" responder"], r$RET[r$Time=="t90" & r$Response==" nonresponder"])

boxplot(data, names=c("0 days, respond", "0 days, non-respond", "10 days, respond", "10 days, non-respond", "90 days, respond", "90 days, non-respond" ))

dev.off()

## Now look at EGR1 expression over time.
egr<- read.table("EGR1.tsv", sep="\t", header=TRUE)

wilcox.test(egr$EGR1[egr$Time == "t0"], egr$EGR1[egr$Time == "t90"])

pdf("EGR1-time.pdf", height=6, width=4)

boxplot(egr$EGR1[egr$Time == "t0"], egr$EGR1[egr$Time == "t90"],
	names=c("pre-treatment", "90 days"), outline=FALSE, ylab="EGR1 expression level", main="EGR1 Expression")
stripchart(list(egr$EGR1[egr$Time == "t0"], egr$EGR1[egr$Time == "t90"]), 
		vertical=TRUE, method="jitter", add=TRUE, pch=20, col="blue", cex=2)

dev.off()


wilcox.test(egr$EGR1[egr$Time=="t0" & egr$Response==" responder"], egr$EGR1[egr$Time=="t0" & egr$Response==" nonresponder"])
wilcox.test(egr$EGR1[egr$Time=="t10" & egr$Response==" responder"], egr$EGR1[egr$Time=="t10" & egr$Response==" nonresponder"])
wilcox.test(egr$EGR1[egr$Time=="t90" & egr$Response==" responder"], egr$EGR1[egr$Time=="t90" & egr$Response==" nonresponder"])

data <- list(egr$EGR1[egr$Time=="t0" & egr$Response==" responder"], egr$EGR1[egr$Time=="t0" & egr$Response==" nonresponder"],
	egr$EGR1[egr$Time=="t10" & egr$Response==" responder"], egr$EGR1[egr$Time=="t10" & egr$Response==" nonresponder"],
	egr$EGR1[egr$Time=="t90" & egr$Response==" responder"], egr$EGR1[egr$Time=="t90" & egr$Response==" nonresponder"])

boxplot(data, names=c("0 days, respond", "0 days, non-respond", "10 days, respond", "10 days, non-respond", "90 days, respond", "90 days, non-respond" ))
