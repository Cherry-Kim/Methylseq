STEP1_rGREAT <- function(){
	library(rGREAT)
	bed = read.table("metilene_Tumor_Normal.filter_qval.0.05.rGREAT.input.bed", header=T)
	dim(bed)
	job = submitGreatJob(bed, includeCuratedRegDoms = FALSE,  species='hg38', rule = "oneClosest")
	res = plotRegionGeneAssociationGraphs(job)
	write.table(res, "rGREAT_gene.out", quote=FALSE, row.names = FALSE)
}

STEP1_ChIPseeker <- function(){
        library(ChIPseeker)
        library(org.Hs.eg.db)
        library(TxDb.Hsapiens.UCSC.hg38.knownGene)
        txdb = TxDb.Hsapiens.UCSC.hg38.knownGene

	bed = read.table("metilene_Tumor_Normal.filter_qval.0.05.rGREAT.input.bed", header=T)
	#head(bed)
  	# chr   start     end
	#1 chr1  205339  205340

        peaks.gr = makeGRangesFromDataFrame(bed, keep.extra.columns=TRUE)
        peakAnno  = annotatePeak(peaks.gr,tssRegion=c(-3000,3000), TxDb=txdb, annoDb="org.Hs.eg.db")
        plotAnnoPie(peakAnno,cex =1.1) #legend size
        peak.anno = as.data.frame(peakAnno)
        write.table(peak.anno, file='metilene_Tumor_Normal.filter_qval.0.05.anno.txt',row.names=FALSE,quote=FALSE,sep='\t')
}
#STEP1_rGREAT()
STEP1_ChIPseeker()
	#python 1-2.anno.py -> 
	data <- read.delim('metilene_Tumor_Normal.filter_qval.0.05.head.out',header=T)
	data2 <- read.delim('metilene_Tumor_Normal.filter_qval.0.05.anno2.txt',header=T)
	a <- cbind(data, data2)
	write.table(a, file='metilene_Tumor_Normal.filter_qval.0.05.anno3.txt',row.names=FALSE,quote=FALSE,sep='\t')
