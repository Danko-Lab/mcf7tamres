##
## Gets motifs enriched in TREs that change in TAM resistence lines.

load("APCluster.rdata")

tres <- read.table("GDNF_treatment.dREG-HD.tsv", header=TRUE)
#tres <- tres[tres$V5>0.85,]

sort_bed <- function(bed) {
  require(gtools)
  ord <- mixedorder(paste(bed[,1],bed[,2],bed[,3], sep="-"))
  return(bed[ord,])
}

tre_cng_1h <- tres[!is.na(tres$FDR_GDNF_1h) & tres$FDR_GDNF_1h < 0.01,]
tre_cng_24h <- tres[!is.na(tres$FDR_GDNF_1h) & tres$FDR_GDNF_24h < 0.01,] 

tre_unc <- tres[!is.na(tres$FDR_GDNF_1h) & !is.na(tres$FDR_GDNF_24h) & tres$FDR_GDNF_1h > 0.2 & tres$FDR_GDNF_24h > 0.2 & abs(tres$FC_GDNF_1h)< 0.5 & abs(tres$FC_GDNF_24h)< 0.5, ]

#source("center_dREG_site.R")
tre_cng_1h <- sort_bed(tre_cng_1h) #center_dREG_site(tre_up_tam, "mcf7.plus.bw", "mcf7.minus.bw", readThresh=10)
tre_cng_24h <- sort_bed(tre_cng_24h) #center_dREG_site(tre_dn_tam, "mcf7.plus.bw", "mcf7.minus.bw", readThresh=10)
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

## Now do the motif search

options=list(abline = NULL,title  = "",xlab   = "Order",ylab   = "-log10(p-value)",y.max  = NULL,top.motif.labels = 5,bottom.motif.labels = 5, color.scheme = 2, width = 4, height = 7, zoom.tick = 1, zoom.label = 1,zoom.motif.logo = 1.5, zoom.legend.label=1,zoom.motif.label = 1 );

#for(i in c(7, 7.5, 8, 9)) {
 i <- 7.5

 motif_1h  <- tfbs.enrichmentTest(tfs, file.twoBit= hg19, positive.bed= tre_cng_1h, negative.bed= tre_unc, threshold=i, use.cluster=TRUE, gc.correction=TRUE, ncores=4, pv.adj="bonferroni", gc.min.sample= 500)
 motif_24h <- tfbs.enrichmentTest(tfs, file.twoBit= hg19, positive.bed= tre_cng_24h, negative.bed= tre_unc, threshold=i, use.cluster=TRUE, gc.correction=TRUE, ncores=4, pv.adj="bonferroni", gc.min.sample= 500)

 print(paste("Log score: ",i,sep=""))

 print("Motifs enriched in TREs which change at 1h.")
 tfbs.reportEnrichment(tfs, motif_1h, file.pdf=paste("Motif.GDNF-1h.",i,".pdf", sep=""), sig.only=TRUE, report.title="TEST FULL", enrichment.type="enriched", pv.threshold= 0.1);
 tfbs.plotEnrichment(tfs, motif_1h, file.pdf=paste("Motif.GDNF-1h.QQ.",i,".pdf", sep=""), enrichment.type="enriched", options= options)

 print("TRE DOWN in Tam Res.")
 tfbs.reportEnrichment(tfs, motif_24h, file.pdf=paste("Motif.GDNF-24h.",i,".pdf", sep=""), sig.only=TRUE, report.title="TEST FULL", enrichment.type="enriched", pv.threshold= 0.1);
 tfbs.plotEnrichment(tfs, motif_24h, file.pdf=paste("Motif.GDNF-24h.QQ.",i,".pdf", sep=""), enrichment.type="enriched", options= options)

#}


