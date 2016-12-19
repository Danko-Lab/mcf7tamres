## Compare fold-changes following GDNF treatment between TamS and TamR lines.
gdnf <- read.table("Signif.Changes.TamRes.tsv", header=TRUE)

indx_1h <- (!is.na(gdnf[,"FDR_GDNF_1h_S"]) & gdnf[,"FDR_GDNF_1h_S"]< 0.01) | (!is.na(gdnf[,"FDR_GDNF_1h_R"]) & gdnf[,"FDR_GDNF_1h_R"]< 0.01)
indx_24h <- (!is.na(gdnf[,"FDR_GDNF_24h_S"]) & gdnf[,"FDR_GDNF_24h_S"]< 0.01) | (!is.na(gdnf[,"FDR_GDNF_24h_R"]) & gdnf[,"FDR_GDNF_24h_R"]< 0.01)

pdf("S5C.GDNF-Fold-Changes.TamS-TamR.pdf")

cor.test(gdnf[indx_1h, "FC_GDNF_1h_S"], gdnf[indx_1h, "FC_GDNF_1h_R"])
plot(gdnf[indx_1h, "FC_GDNF_1h_S"], gdnf[indx_1h, "FC_GDNF_1h_R"], pch=19, xlab="Log-2 fold-change 1h GDNF, TamS", ylab="Log-2 fold-change 1h GDNF, TamR"); abline(0,1); abline(h=0); abline(v=0)

cor.test(gdnf[indx_24h, "FC_GDNF_24h_S"], gdnf[indx_24h, "FC_GDNF_24h_R"])
plot(gdnf[indx_24h, "FC_GDNF_24h_S"], gdnf[indx_24h, "FC_GDNF_24h_R"], pch=19, xlab="Log-2 fold-change 24h GDNF, TamS", ylab="Log-2 fold-change 24h GDNF, TamR"); abline(0,1); abline(h=0); abline(v=0)

dev.off()

