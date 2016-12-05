#################################################################################
# Intersects genes for use later.

for (mode in c("results_de", "results_lfc")) {
    if (mode=="results_de") {
        extra <- "de"
    } else {
        extra <- "lfc"
    }
    
    Ambion <- read.table(file.path(mode, "siRNA_289_vs_Ambion.txt"), header=TRUE, stringsAsFactors=FALSE)
    Dharmacon <- read.table(file.path(mode, "siRNA_289_vs_Dharmacon.txt"), header=TRUE, stringsAsFactors=FALSE)
    m <- match(Ambion$ENSEMBL, Dharmacon$ENSEMBL)
    Dharmacon <- Dharmacon[m,]
    pval <- pmax(Ambion$P.Value, Dharmacon$P.Value)
    o <- order(pval)
    write.table(file=file.path(mode, "combined_siRNA_289.txt"), row.names=FALSE, sep="\t", quote=FALSE,
                data.frame(Ambion[,c("ENSEMBL", "SYMBOL")], Ambion.logFC=Ambion$logFC, Dharmacon.logFC=Dharmacon$logFC,
                           P.Value=pval, adj.P.Val=p.adjust(pval, method="BH"))[o,])
 
    LNAA <- read.table(file.path(mode, "LNA_289.2_vs_controlA.txt"), header=TRUE, stringsAsFactors=FALSE)
    LNAB <- read.table(file.path(mode, "LNA_289.2_vs_controlB.txt"), header=TRUE, stringsAsFactors=FALSE)
    m <- match(LNAA$ENSEMBL, LNAB$ENSEMBL)
    LNAB <- LNAB[m,]
    pval <- pmax(LNAA$P.Value, LNAB$P.Value)
    o <- order(pval)
    write.table(file=file.path(mode, "combined_LNA_289.txt"), row.names=FALSE, sep="\t", quote=FALSE,
                data.frame(LNAA[,c("ENSEMBL", "SYMBOL")], LNAA.logFC=LNAA$logFC, LNAB.logFC=LNAB$logFC,
                           P.Value=pval, adj.P.Val=p.adjust(pval, method="BH"))[o,])
   
    guide1 <- read.table(file.path(mode, "CRISPRi_289.1_vs_negguide2.txt"), header=TRUE, stringsAsFactors=FALSE)
    guide7 <- read.table(file.path(mode, "CRISPRi_289.9_vs_negguide2.txt"), header=TRUE, stringsAsFactors=FALSE)
    m <- match(guide1$ENSEMBL, guide7$ENSEMBL)
    guide7 <- guide7[m,]
    pval <- pmax(guide1$P.Value, guide7$P.Value)
    o <- order(pval)
    write.table(file=file.path(mode, "combined_CRISPRi_289.txt"), row.names=FALSE, sep="\t", quote=FALSE,
                data.frame(guide1[,c("ENSEMBL", "SYMBOL")], Guide1.logFC=guide1$logFC, Guide7.logFC=guide7$logFC,
                           P.Value=pval, adj.P.Val=p.adjust(pval, method="BH"))[o,])

    h19v1 <- read.table(file.path(mode, "CRISPRi_H19.2_vs_negguide1.txt"), header=TRUE, stringsAsFactors=FALSE)
    h19v2 <- read.table(file.path(mode, "CRISPRi_H19.2_vs_negguide2.txt"), header=TRUE, stringsAsFactors=FALSE)
    m <- match(h19v1$ENSEMBL, h19v2$ENSEMBL)
    h19v2 <- h19v2[m,]
    pval <- pmax(h19v1$P.Value, h19v2$P.Value)
    o <- order(pval)
    write.table(file=file.path(mode, "combined_CRISPRi_H19.txt"), row.names=FALSE, sep="\t", quote=FALSE,
                data.frame(h19v1[,c("ENSEMBL", "SYMBOL")], Versus1.logFC=h19v1$logFC, Versus2.logFC=h19v2$logFC,
                           P.Value=pval, adj.P.Val=p.adjust(pval, method="BH"))[o,])
}




