##
## Gets motifs enriched in TREs that change in TAM resistence lines.

tres <- read.table("TRE.Signif.Changes.TamRes.tsv", header=TRUE)
#tres <- tres[tres$V5>0.85,]

tre_up_tam <- tres[tres$FDR_TAM < 0.01 & tres$FC_TAM < 0,]
tre_dn_tam <- tres[tres$FDR_TAM < 0.01 & tres$FC_TAM > 0,]

tre_unc <- tres[tres$FDR_TAM > 0.2 & abs(tres$FC_TAM)< 0.25, ]



#source("center_dREG_site.R")
tre_up <- tre_up_tam #center_dREG_site(tre_up_tam, "mcf7.plus.bw", "mcf7.minus.bw", readThresh=10)
tre_dn <- tre_dn_tam #center_dREG_site(tre_dn_tam, "mcf7.plus.bw", "mcf7.minus.bw", readThresh=10)
tre_unc <- tre_unc   #center_dREG_site(tre_unc, "mcf7.plus.bw", "mcf7.minus.bw", readThresh=10)

require(rtfbsdb)
hg19 <- "/local/storage/data/hg19/hg19.2bit"

db.human <- CisBP.extdata("Homo_sapiens")
tfs <- tfbs.createFromCisBP(db.human, file.bigwig.plus= "mcf7.plus.bw", file.bigwig.minus= "mcf7.minus.bw", file.gencode.gtf="/local/storage/data/hg19/all/gencode/gencode.v19.annotation.gtf.gz")
tfs <- tfbs.getDistanceMatrix(tfs, ncores=25)
clu <- tfbs.clusterMotifs(tfs, method="agnes", pdf.heatmap="motifs/MCF7.heatmap.pdf")
use_indx <- tfbs.selectByGeneExp(tfs, clu) 

for(i in c(7, 8, 9)) {
 motif_up <- tfbs.compareTFsite(tfs, file.twoBit= hg19, positive.bed= tre_up, negative.bed= tre_unc, use.cluster=clu, threshold=i, gc.correction=TRUE, pv.adj="bonferroni")
 motif_dn <- tfbs.compareTFsite(tfs, file.twoBit= hg19, positive.bed= tre_dn, negative.bed= tre_unc, use.cluster=clu, threshold=i, gc.correction=TRUE, pv.adj="bonferroni")

 print(paste("Log score: ",i,sep=""))

 print("TRE UP in Tam Res.")
 print(head(motif_up$result[motif_up$result$es.ratio>1,][order(motif_up$result$pv.adj[motif_up$result$es.ratio>1]),], 10))

 print("TRE DOWN in Tam Res.")
 print(head(motif_dn$result[motif_dn$result$es.ratio>1,][order(motif_dn$result$pv.adj[motif_dn$result$es.ratio>1]),], 10))

}

save.image("MCF7db.RData")

## Draw logos for several genes.

pdf("logos.pdf")
 tfbs.drawLogo(tfs, index=which(tfs@TFID=="M4376_1.01"))
 tfbs.drawLogo(tfs, index=which(tfs@TFID=="M4440_1.01"))
 tfbs.drawLogo(tfs, index=which(tfs@TFID=="M4400_1.01"))
dev.off()

## TFAP2C
motif_up$result[motif_up$result$tf.name=="TFAP2C",]
motif_dn$result[motif_dn$result$tf.name=="TFAP2C",]

tfbs.drawLogo(tfs, index=which(tfs@TFID=="M2190_1.01"))

