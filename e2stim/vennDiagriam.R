require(VennDiagram)

#################################
## E2+BBCA
bbca <- read.table("BBCA.Responding.tsv", header=TRUE); bbca <- bbca[bbca$FDR_bbca < 0.01,]
NROW(unique(bbca$V7)) ## 614

e2 <- read.table("GSE27463_RefSeq.reg.tsv.gz", header=TRUE); e2 <- e2[e2$E2.40m.qVal < 0.01,]
NROW(unique(e2$GeneSymbol)) ## 2445

both<- list(bbca=unique(bbca$V7), e2=unique(e2$GeneSymbol))

vd <- calculate.overlap(both)
NROW(vd$a3) ##129

venn.diagram(both, "E2.BBCA.Venn.png", imagetype="png")

#################################
## BBCA+GDNF 1 hour
gdnf <- read.table("../gdnfstim/GDNF_treatment.tsv", header=TRUE); 
gdnf1 <- gdnf[gdnf$FDR_GDNF_1h_S < 0.01,]
NROW(unique(gdnf1$V7)) ## 2170

both<- list(bbca=unique(bbca$V7), gdnf_1h=unique(gdnf1$V7))
vd <- calculate.overlap(both)
NROW(vd$a3) ##121
venn.diagram(both, "GDNF_1h.BBCA.Venn.png", imagetype="png")

#################################
## Aaand GDNF24 hours.
gdnf24<- gdnf[gdnf$FDR_GDNF_24h_S < 0.01,]
NROW(unique(gdnf24$V7)) ## 821

both<- list(bbca=unique(bbca$V7), gdnf_24h=unique(gdnf24$V7))
vd <- calculate.overlap(both)
NROW(vd$a3) ##164
venn.diagram(both, "GDNF_24h.BBCA.Venn.png", imagetype="png")

#################################
## GDNF24 and E2
both<- list(e2=unique(e2$GeneSymbol), gdnf_24h=unique(gdnf24$V7))
vd <- calculate.overlap(both)
NROW(vd$a3) ##219
venn.diagram(both, "GDNF_24h.E2.Venn.png", imagetype="png")

#################################
## GDNF1 and GDNF24
both<- list(gdnf_1h=unique(gdnf1$V7), gdnf_24h=unique(gdnf24$V7))
vd <- calculate.overlap(both)
NROW(vd$a3) ##227
venn.diagram(both, "GDNF_24h.GDNF_1h.Venn.png", imagetype="png")

#################################
## ALL Conditions
venn.diagram(list(gdnf_1h=unique(gdnf1$V7), gdnf_24h=unique(gdnf24$V7),bbca=unique(bbca$V7), e2=unique(e2$GeneSymbol)), "All.Venn.png", imagetype="png")

venn.diagram(list(gdnf_1h=unique(gdnf1$V7), gdnf_24h=unique(gdnf24$V7),e2=unique(e2$GeneSymbol)), "GDNF.1.24.E2.Venn.png", imagetype="png")




