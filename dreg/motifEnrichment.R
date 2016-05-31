##
## Gets motifs enriched in TREs that change in TAM resistence lines.

tres <- read.table("TRE.Signif.Changes.TamRes.tsv", header=TRUE)
#tres <- tres[tres$V5>0.85,]

sort_bed <- function(bed) {
  require(gtools)
  ord <- mixedorder(paste(bed[,1],bed[,2],bed[,3], sep="-"))
  return(bed[ord,])
}

tre_up_tam <- tres[tres$FDR_TAM < 0.01 & tres$FC_TAM < 0,]
tre_dn_tam <- tres[tres$FDR_TAM < 0.01 & tres$FC_TAM > 0,]

tre_unc <- tres[tres$FDR_TAM > 0.2 & abs(tres$FC_TAM)< 0.25, ]

#source("center_dREG_site.R")
tre_up <- sort_bed(tre_up_tam) #center_dREG_site(tre_up_tam, "mcf7.plus.bw", "mcf7.minus.bw", readThresh=10)
tre_dn <- sort_bed(tre_dn_tam) #center_dREG_site(tre_dn_tam, "mcf7.plus.bw", "mcf7.minus.bw", readThresh=10)
tre_unc <- sort_bed(tre_unc)   #center_dREG_site(tre_unc, "mcf7.plus.bw", "mcf7.minus.bw", readThresh=10)

require(rtfbsdb)

createDB <- function() {
 hg19 <- "/local/storage/data/hg19/hg19.2bit"
 
 db.human <- CisBP.extdata("Homo_sapiens")
 tfs <- tfbs.createFromCisBP(db.human);
 
 # Clustering... and expression
 tfs <- tfbs.clusterMotifs(tfs, method="apcluster", ncores=10)
 tfs <- tfbs.getExpression(tfs, hg19, "/local/storage/data/hg19/all/gencode/gencode.v19.annotation.gtf.gz", "mcf7.plus.bw", "mcf7.minus.bw", ncores=10);
 tfs <- tfbs.selectByGeneExp(tfs)
 
 save.image(file="APCluster.rdata")
}

load("APCluster.rdata")

options=list(abline = NULL,title  = "",xlab   = "Order",ylab   = "-log10(p-value)",y.max  = NULL,top.motif.labels = 5,bottom.motif.labels = 5, color.scheme = 2, width = 4, height = 7, zoom.tick = 1, zoom.label = 1,zoom.motif.logo = 1.5, zoom.legend.label=1,zoom.motif.label = 1 );

for(i in c(7, 7.5, 8, 9)) {
 motif_cng<- tfbs.enrichmentTest(tfs, file.twoBit= hg19, positive.bed= sort_bed(rbind(tre_up, tre_dn)), negative.bed= tre_unc, threshold=i, use.cluster=TRUE, gc.correction=TRUE, ncores=20, pv.adj="bonferroni")
 motif_up <- tfbs.enrichmentTest(tfs, file.twoBit= hg19, positive.bed= tre_up, negative.bed= tre_unc, threshold=i, use.cluster=TRUE, gc.correction=TRUE, ncores=20, pv.adj="bonferroni")
 motif_dn <- tfbs.enrichmentTest(tfs, file.twoBit= hg19, positive.bed= tre_dn, negative.bed= tre_unc, threshold=i, use.cluster=TRUE, gc.correction=TRUE, ncores=20, pv.adj="bonferroni")

 print(paste("Log score: ",i,sep=""))

 print("TREs enriched in any changed sute following Tam Res.")
 print(head(motif_cng$result[motif_cng$result$fe.ratio>1,][order(motif_cng$result$pv.adj[motif_cng$result$fe.ratio>1]),], 10))
 tfbs.reportEnrichment(tfs, motif_cng, file.pdf=paste("Motif.cng.",i,".pdf", sep=""), sig.only=TRUE, report.title="TEST FULL", enrichment.type="enriched", pv.threshold= 0.1);
 tfbs.plotEnrichment(tfs, motif_cng, file.pdf=paste("Motif.cng.QQ.",i,".pdf", sep=""), enrichment.type="enriched", options= options)

 print("TRE UP in Tam Res.")
 print(head(motif_up$result[motif_up$result$fe.ratio>1,][order(motif_up$result$pv.adj[motif_up$result$fe.ratio>1]),], 10))
 tfbs.reportEnrichment(tfs, motif_up, file.pdf=paste("Motif.up.",i,".pdf", sep=""), sig.only=TRUE, report.title="TEST FULL", enrichment.type="enriched", pv.threshold= 0.1);
 tfbs.plotEnrichment(tfs, motif_up, file.pdf=paste("Motif.up.QQ.",i,".pdf", sep=""), enrichment.type="enriched", options= options)

 print("TRE DOWN in Tam Res.")
 print(head(motif_dn$result[motif_dn$result$fe.ratio>1,][order(motif_dn$result$pv.adj[motif_dn$result$fe.ratio>1]),], 10))
 tfbs.reportEnrichment(tfs, motif_dn, file.pdf=paste("Motif.dn.",i,".pdf", sep=""), sig.only=TRUE, report.title="TEST FULL", enrichment.type="enriched", pv.threshold= 0.1);
 tfbs.plotEnrichment(tfs, motif_dn, file.pdf=paste("Motif.dn.QQ.",i,".pdf", sep=""), enrichment.type="enriched", options= options)
}

save.image("MCF7db.RData")

q("no")

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

