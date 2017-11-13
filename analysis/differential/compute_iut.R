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
    m <- match(rownames(Ambion), rownames(Dharmacon))
    Dharmacon <- Dharmacon[m,]
    pval <- pmax(Ambion$P.Value, Dharmacon$P.Value)
    o <- order(pval)
    write.table(file=file.path(mode, "combined_siRNA_289.txt"), col.names=NA, sep="\t", quote=FALSE,
                data.frame(Symbol=Ambion$Symbol, Ambion.logFC=Ambion$logFC, Dharmacon.logFC=Dharmacon$logFC,
                           P.Value=pval, adj.P.Val=p.adjust(pval, method="BH"), row.names=rownames(Ambion))[o,])
 
    LNAA <- read.table(file.path(mode, "LNA_289.2_vs_controlA.txt"), header=TRUE, stringsAsFactors=FALSE)
    LNAB <- read.table(file.path(mode, "LNA_289.2_vs_controlB.txt"), header=TRUE, stringsAsFactors=FALSE)
    m <- match(rownames(LNAA), rownames(LNAB))
    LNAB <- LNAB[m,]
    pval <- pmax(LNAA$P.Value, LNAB$P.Value)
    o <- order(pval)
    write.table(file=file.path(mode, "combined_LNA_289.txt"), col.names=NA, sep="\t", quote=FALSE,
                data.frame(Symbol=LNAA$Symbol, LNAA.logFC=LNAA$logFC, LNAB.logFC=LNAB$logFC,
                           P.Value=pval, adj.P.Val=p.adjust(pval, method="BH"), row.names=rownames(LNAA))[o,])
   
    guide1 <- read.table(file.path(mode, "CRISPRi_289.1_vs_negguide2.txt"), header=TRUE, stringsAsFactors=FALSE)
    guide7 <- read.table(file.path(mode, "CRISPRi_289.9_vs_negguide2.txt"), header=TRUE, stringsAsFactors=FALSE)
    m <- match(rownames(guide1), rownames(guide7))
    guide7 <- guide7[m,]
    pval <- pmax(guide1$P.Value, guide7$P.Value)
    o <- order(pval)
    write.table(file=file.path(mode, "combined_CRISPRi_289.txt"), col.names=NA, sep="\t", quote=FALSE,
                data.frame(Symbol=guide1$Symbol, Guide1.logFC=guide1$logFC, Guide7.logFC=guide7$logFC,
                           P.Value=pval, adj.P.Val=p.adjust(pval, method="BH"), row.names=rownames(guide1))[o,])

    guide1_het <- read.table(file.path(mode, "CRISPRi_het_289.9_vs_Nc.1.txt"), header=TRUE, stringsAsFactors=FALSE)
    guide2_het <- read.table(file.path(mode, "CRISPRi_het_289.9_vs_Nc.2.txt"), header=TRUE, stringsAsFactors=FALSE)
    m <- match(rownames(guide1_het), rownames(guide2_het))
    guide2_het <- guide2_het[m,]
    pval <- pmax(guide1_het$P.Value, guide2_het$P.Value)
    o <- order(pval)
    write.table(file=file.path(mode, "combined_CRISPRi_het_289.txt"), col.names=NA, sep="\t", quote=FALSE,
                data.frame(Symbol=guide1_het$Symbol, Versus1.logFC=guide1_het$logFC, Versus2.logFC=guide2_het$logFC,
                           P.Value=pval, adj.P.Val=p.adjust(pval, method="BH"), row.names=rownames(guide1_het))[o,])

    h19v1 <- read.table(file.path(mode, "CRISPRi_H19.2_vs_negguide1.txt"), header=TRUE, stringsAsFactors=FALSE)
    h19v2 <- read.table(file.path(mode, "CRISPRi_H19.2_vs_negguide2.txt"), header=TRUE, stringsAsFactors=FALSE)
    m <- match(rownames(h19v1), rownames(h19v2))
    h19v2 <- h19v2[m,]
    pval <- pmax(h19v1$P.Value, h19v2$P.Value)
    o <- order(pval)
    write.table(file=file.path(mode, "combined_CRISPRi_H19.txt"), col.names=NA, sep="\t", quote=FALSE,
                data.frame(Symbol=h19v1$Symbol, Versus1.logFC=h19v1$logFC, Versus2.logFC=h19v2$logFC,
                           P.Value=pval, adj.P.Val=p.adjust(pval, method="BH"), row.names=rownames(h19v1))[o,])

    h19v1_het <- read.table(file.path(mode, "CRISPRi_het_H19_vs_Nc.1.txt"), header=TRUE, stringsAsFactors=FALSE)
    h19v2_het <- read.table(file.path(mode, "CRISPRi_het_H19_vs_Nc.2.txt"), header=TRUE, stringsAsFactors=FALSE)
    m <- match(rownames(h19v1_het), rownames(h19v2_het))
    h19v2_het <- h19v2_het[m,]
    pval <- pmax(h19v1_het$P.Value, h19v2_het$P.Value)
    o <- order(pval)
    write.table(file=file.path(mode, "combined_CRISPRi_het_H19.txt"), col.names=NA, sep="\t", quote=FALSE,
                data.frame(Symbol=h19v1_het$Symbol, Versus1.logFC=h19v1_het$logFC, Versus2.logFC=h19v2_het$logFC,
                           P.Value=pval, adj.P.Val=p.adjust(pval, method="BH"), row.names=rownames(h19v1_het))[o,])
}




