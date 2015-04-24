##
## Gets motifs enriched in TREs that change in TAM resistence lines.

tres <- read.table("TRE.Signif.Changes.TamRes.tsv", header=TRUE)
tres <- tres[tres$V5>0.85,]

tre_up_tam <- tres[tres$FDR_TAM < 0.01 & tres$FC_TAM < 0,]
tre_dn_tam <- tres[tres$FDR_TAM < 0.01 & tres$FC_TAM > 0,]

tre_unc <- tres[tres$FDR_TAM > 0.2 & abs(tres$FC_TAM)< 0.25, ]

source("center_dREG_site.R")
tre_up <- center_dREG_site(tre_up_tam, "mcf7.plus.bw", "mcf7.minus.bw", readThresh=10)
tre_dn <- center_dREG_site(tre_dn_tam, "mcf7.plus.bw", "mcf7.minus.bw", readThresh=10)

tre_unc <- center_dREG_site(tre_unc, "mcf7.plus.bw", "mcf7.minus.bw", readThresh=10)


require(rtfbsdb)

db.human <- CisBP.extdata("Homo_sapiens")
tfs <- CisBP.find(db.human)

tfs <- tfbs.getDistanceMatrix(tfs, ncores=25)
tfs <- tfbs.getExpression(tfs, "mcf7.plus.bw", "mcf7.minus.bw")   ## Get expressed TFs ... 

tfs <- tfbs.clusterMotifs(tfs, method="agnes", pdf.heatmap="motifs/MCF7.heatmap.pdf", group.k= )

comparative_scanDb_rtfbs()
