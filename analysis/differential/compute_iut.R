#################################################################################
# Intersects genes for use later.

getIUTP <- function(tab1, tab2) {
    lfc1 <- tab1$logFC
    lfc2 <- tab2$logFC
    p1 <- tab1$P.Value
    p2 <- tab2$P.Value

    # Converting into one-sided p-values.
    up.p1 <- ifelse(lfc1 > 0, p1/2, 1-p1/2)
    up.p2 <- ifelse(lfc2 > 0, p2/2, 1-p2/2)
    down.p1 <- ifelse(lfc1 < 0, p1/2, 1-p1/2) # technically '1 - up.p1', but avoid numerical imprecision.
    down.p2 <- ifelse(lfc2 < 0, p2/2, 1-p2/2)
    
    # Combining and then converting back to two-sided p-values via '*2'.
    up.p <- pmax(up.p1, up.p2)
    down.p <- pmax(down.p1, down.p2)
    pval <- pmin(up.p, down.p, 0.5)*2

    return(data.frame(
        AveExpr=(tab1$AveExpr + tab2$AveExpr)/2,
        logFC=(tab1$logFC + tab2$logFC)/2,
        P.Value=pval))
}

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
    combined <- getIUTP(Ambion, Dharmacon)
    o <- order(combined$P.Value)
    write.table(file=file.path(mode, "combined_siRNA_289.txt"), col.names=NA, sep="\t", quote=FALSE,
                data.frame(Symbol=Ambion$Symbol, Ambion.logFC=Ambion$logFC, Dharmacon.logFC=Dharmacon$logFC,
                           combined, adj.P.Val=p.adjust(combined$P.Value, method="BH"), 
                           row.names=rownames(Ambion))[o,])
 
    LNAA <- read.table(file.path(mode, "LNA_289.2_vs_controlA.txt"), header=TRUE, stringsAsFactors=FALSE)
    LNAB <- read.table(file.path(mode, "LNA_289.2_vs_controlB.txt"), header=TRUE, stringsAsFactors=FALSE)
    m <- match(rownames(LNAA), rownames(LNAB))
    LNAB <- LNAB[m,]
    combined <- getIUTP(LNAA, LNAB)
    o <- order(combined$P.Value)
    write.table(file=file.path(mode, "combined_LNA_289.txt"), col.names=NA, sep="\t", quote=FALSE,
                data.frame(Symbol=LNAA$Symbol, LNAA.logFC=LNAA$logFC, LNAB.logFC=LNAB$logFC,
                           combined, adj.P.Val=p.adjust(combined$P.Value, method="BH"), 
                           row.names=rownames(LNAA))[o,])
   
    guide1 <- read.table(file.path(mode, "CRISPRi_289.1_vs_negguide2.txt"), header=TRUE, stringsAsFactors=FALSE)
    guide7 <- read.table(file.path(mode, "CRISPRi_289.9_vs_negguide2.txt"), header=TRUE, stringsAsFactors=FALSE)
    m <- match(rownames(guide1), rownames(guide7))
    guide7 <- guide7[m,]
    combined <- getIUTP(guide1, guide7)
    o <- order(combined$P.Value)
    write.table(file=file.path(mode, "combined_CRISPRi_289.txt"), col.names=NA, sep="\t", quote=FALSE,
                data.frame(Symbol=guide1$Symbol, Guide1.logFC=guide1$logFC, Guide7.logFC=guide7$logFC,
                           combined, adj.P.Val=p.adjust(combined$P.Value, method="BH"), 
                           row.names=rownames(guide1))[o,])

    guide1_het <- read.table(file.path(mode, "CRISPRi_het_289.9_vs_Nc.1.txt"), header=TRUE, stringsAsFactors=FALSE)
    guide2_het <- read.table(file.path(mode, "CRISPRi_het_289.9_vs_Nc.2.txt"), header=TRUE, stringsAsFactors=FALSE)
    m <- match(rownames(guide1_het), rownames(guide2_het))
    guide2_het <- guide2_het[m,]
    combined <- getIUTP(guide1_het, guide2_het)
    o <- order(combined$P.Value)
    write.table(file=file.path(mode, "combined_CRISPRi_het_289.txt"), col.names=NA, sep="\t", quote=FALSE,
                data.frame(Symbol=guide1_het$Symbol, Versus1.logFC=guide1_het$logFC, Versus2.logFC=guide2_het$logFC,
                           combined, adj.P.Val=p.adjust(combined$P.Value, method="BH"), 
                           row.names=rownames(guide1_het))[o,])

    h19v1 <- read.table(file.path(mode, "CRISPRi_H19.2_vs_negguide1.txt"), header=TRUE, stringsAsFactors=FALSE)
    h19v2 <- read.table(file.path(mode, "CRISPRi_H19.2_vs_negguide2.txt"), header=TRUE, stringsAsFactors=FALSE)
    m <- match(rownames(h19v1), rownames(h19v2))
    h19v2 <- h19v2[m,]
    combined <- getIUTP(h19v1, h19v2)
    o <- order(combined$P.Value)
    write.table(file=file.path(mode, "combined_CRISPRi_H19.txt"), col.names=NA, sep="\t", quote=FALSE,
                data.frame(Symbol=h19v1$Symbol, Versus1.logFC=h19v1$logFC, Versus2.logFC=h19v2$logFC,
                           combined, adj.P.Val=p.adjust(combined$P.Value, method="BH"), 
                           row.names=rownames(h19v1))[o,])

    h19v1_het <- read.table(file.path(mode, "CRISPRi_het_H19_vs_Nc.1.txt"), header=TRUE, stringsAsFactors=FALSE)
    h19v2_het <- read.table(file.path(mode, "CRISPRi_het_H19_vs_Nc.2.txt"), header=TRUE, stringsAsFactors=FALSE)
    m <- match(rownames(h19v1_het), rownames(h19v2_het))
    h19v2_het <- h19v2_het[m,]
    combined <- getIUTP(h19v1_het, h19v2_het)
    o <- order(combined$P.Value)
    write.table(file=file.path(mode, "combined_CRISPRi_het_H19.txt"), col.names=NA, sep="\t", quote=FALSE,
                data.frame(Symbol=h19v1_het$Symbol, Versus1.logFC=h19v1_het$logFC, Versus2.logFC=h19v2_het$logFC,
                           combined, adj.P.Val=p.adjust(combined$P.Value, method="BH"), 
                           row.names=rownames(h19v1_het))[o,])
}




