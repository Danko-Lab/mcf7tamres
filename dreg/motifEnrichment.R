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
hg19 <- "/local/storage/data/hg19/hg19.2bit"

db.human <- CisBP.extdata("Homo_sapiens")
tfs <- tfbs.createFromCisBP(db.human, file.bigwig.plus= "mcf7.plus.bw", file.bigwig.minus= "mcf7.minus.bw", file.gencode.gtf="/local/storage/data/hg19/all/gencode/gencode.v19.annotation.gtf.gz")
tfs <- tfbs.getDistanceMatrix(tfs, ncores=25)
clu <- tfbs.clusterMotifs(tfs, method="agnes", pdf.heatmap="motifs/MCF7.heatmap.pdf")

for(i in c(5, 6, 7, 8, 9, 10)) {
 motif_up <- tfbs.compareTFsite(tfs, hg19, tre_up, tre_unc, background.correction=TRUE, threshold=i)
 motif_dn <- tfbs.compareTFsite(tfs, hg19, tre_dn, tre_unc, background.correction=TRUE, threshold=i)

 print(paste("Log score: ",i,sep=""))

 print("TRE UP in Tam Res.")
 print(head(motif_up[motif_up$es.ratio>1,][order(motif_up$assoc.pvalue[motif_up$es.ratio>1]),], 10))

 print("TRE DOWN in Tam Res.")
 print(head(motif_dn[motif_dn$es.ratio>1,][order(motif_dn$assoc.pvalue[motif_dn$es.ratio>1]),], 10))
}

save.image("MCF7db.RData")

## Draw logos for several genes.

pdf("logos.pdf")
 tfbs.drawLogo(tfs, index=which(tfs@TFID=="M4376_1.01"))
 tfbs.drawLogo(tfs, index=which(tfs@TFID=="M4440_1.01"))
 tfbs.drawLogo(tfs, index=which(tfs@TFID=="M4400_1.01"))
dev.off()


