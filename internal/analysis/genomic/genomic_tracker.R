stuff <- read.table("../differential/results_lfc/combined_clones.txt", header=TRUE,stringsAsFactors=FALSE)
chosen <- stuff$ENSEMBL[stuff$adj.P.Val==0]

################### Integrating ChIP-seq data (O'Geen) #######################

require(GenomicFeatures)
# txdb <- makeTxDbFromUCSC(genome='hg19',tablename='ensGene') 
# saveDb(txdb, "hg19.ensGene")
txdb <- loadDb("hg19.ensGene")

out <- genes(txdb)
out <- out[names(out) %in% chosen]
strand(out) <- "*"
length(out)

# Pulling out cas9, checking whether or not they overlap the genes.
cas9.data <- read.table("GSM1496583_293T_noED_VEGFA_rep1_peaks.txt.gz")
peak.sitesI <- GRanges(cas9.data[,1], IRanges(cas9.data[,2], cas9.data[,3]))
cas9.data <- read.table("GSM1496584_293T_noED_VEGFA_rep2_peaks.txt.gz")
peak.sitesII <- GRanges(cas9.data[,1], IRanges(cas9.data[,2], cas9.data[,3]))

peak.sites <- c(peak.sitesI, peak.sitesII)
findOverlaps(peak.sites, resize(out, fix="center", width(out)+2e4)) # No overlaps within 10 kbp on either side.

################### Integrating ChIP-seq data (Thakore)  #######################

# Pulling out Thakore peaks.
all.peaks <- list()
for (x in list.files("Thakore_peaks", pattern="_peaks.xls", full=TRUE)) {
    current <- read.table(x, header=TRUE)
    cur.peaks <- GRanges(current[,1], IRanges(current[,2], current[,3]))
    cur.peaks <- cur.peaks[current[,7] >= 10]
    all.peaks[[length(all.peaks)+1]] <- cur.peaks
}
suppressWarnings(all.peaks <- do.call(c, all.peaks))

# Lifting back to hg19.
require(rtracklayer)
hg.chain <- import.chain("hg38ToHg19.over.chain")
suppressWarnings(peak.sites <- liftOver(all.peaks, hg.chain))

findOverlaps(peak.sites, resize(out, fix="center", width(out)+2e4)) # No overlaps within 10 kbp on either side.

################# Integrating RNA-seq data #########################

thakore.genes <- read.csv("thakore.csv", stringsAsFactors=FALSE, header=FALSE)
intersect(chosen, thakore.genes[,2])

################# Checking if they have common function ##################

require(quantsmooth)
CHR <- as.character(sub("^chr", "", seqnames(out)))
MapInfo <- (start(out) + end(out))/2
pdf("chromosomal.pdf", width=5, height=10)
chrompos<-prepareGenomePlot(data.frame(CHR,MapInfo),paintCytobands = TRUE, units="hg19", organism="hsa", sex=TRUE)
points(chrompos[,2],chrompos[,1]+0.1,pch=25,col="red",bg="red")
dev.off()
