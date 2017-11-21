# Makes volcano plots, highlighting the targeted lncRNA in each case. 
# We use the average log-fold change and the maximum p-value in cases where we have two comparisons.

makeVolcano <- function(pvalues, logFCs, chosen, threshold, ylim=c(0, 30), xlim=c(-6, 6), chosen.col="red", chosen.pch=2) {
    lp <- -log10(pvalues)
    par(mar=c(5.1, 5.1, 4.1, 2.1))
    plot(logFCs, lp, pch=16, col="grey50", cex.axis=1.2, cex.lab=1.4, xlab=expression(Log[2]~"fold change"), 
         ylab=expression("-"*Log[10]~"P-value"), xlim=xlim, ylim=ylim, cex=1.5)
    points(logFCs[chosen], lp[chosen], pch=16, col=chosen.col, cex=chosen.pch)
    abline(v=-0.5, col="red", lwd=2, lty=2)
    abline(v=0.5, col="red", lwd=2, lty=2)
    abline(h=-log10(threshold), col="red", lwd=2, lty=2)
}

for (mode in c("results_de", "results_lfc")) {
    if (mode=="results_de") {
        nextra <- "de"
    } else {
        nextra <- "lfc"
    }

    siRNA <- read.table(file.path(mode, "combined_siRNA_289.txt"), header=TRUE, stringsAsFactors=FALSE)
    LNA <- read.table(file.path(mode, "combined_LNA_289.txt"), header=TRUE, stringsAsFactors=FALSE)
    CRISPRi <- read.table(file.path(mode, "combined_CRISPRi_289.txt"), header=TRUE, stringsAsFactors=FALSE)
    CRISPRi.het <- read.table(file.path(mode, "combined_CRISPRi_het_289.txt"), header=TRUE, stringsAsFactors=FALSE)
 
    combined.p <- c(siRNA$P.Value, LNA$P.Value, CRISPRi$P.Value, CRISPRi.het$P.Value)
    is.sig <- p.adjust(combined.p, method="BH") <= 0.05
    threshold <- max(combined.p[is.sig])
    saveRDS(file=sprintf("threshold_%s_289.rds", nextra), threshold)

    H19 <- read.table(file.path(mode, "combined_CRISPRi_H19.txt"), header=TRUE, stringsAsFactors=FALSE)
    H19.het <- read.table(file.path(mode, "combined_CRISPRi_het_H19.txt"), header=TRUE, stringsAsFactors=FALSE)
    h19.p <- c(H19$P.Value, H19.het$P.Value)
    h19.sig <- p.adjust(h19.p, method="BH") <= 0.05
    h19.threshold <- max(h19.p[h19.sig])
    saveRDS(file=sprintf("threshold_%s_H19.rds", nextra), h19.threshold)

    for (zoom in c(TRUE, FALSE)) {
        if (zoom) {
            extra <- paste0(nextra, "_zoom")
            ylim <- c(0, 5)
        } else {
            extra <- nextra
            ylim <- c(0, 30)
        }
        
        pdf(sprintf("pics/siRNA_volcano_289_%s.pdf", extra))
        makeVolcano(siRNA$P.Value, (siRNA$Ambion.logFC + siRNA$Dharmacon.logFC)/2, 
                    which(rownames(siRNA)=="ENSG00000234771"), threshold=threshold, ylim=ylim)
        dev.off()

        pdf(sprintf("pics/LNA_volcano_289_%s.pdf", extra))
        makeVolcano(LNA$P.Value, (LNA$LNAA.logFC + LNA$LNAB.logFC)/2, 
                    which(rownames(LNA)=="ENSG00000234771"), threshold=threshold, ylim=ylim)
        dev.off()

        pdf(sprintf("pics/CRISPRi_volcano_289_%s.pdf", extra))
        makeVolcano(CRISPRi$P.Value, (CRISPRi$Guide1.logFC + CRISPRi$Guide7.logFC)/2,
                    which(rownames(CRISPRi)=="ENSG00000234771"), threshold=threshold, ylim=ylim)
        dev.off()

        pdf(sprintf("pics/CRISPRi_het_volcano_289_%s.pdf", extra))
        makeVolcano(CRISPRi.het$P.Value, (CRISPRi.het$Versus1.logFC + CRISPRi.het$Versus2.logFC)/2,
                    which(rownames(CRISPRi.het)=="ENSG00000234771"), threshold=threshold, ylim=ylim)
        dev.off()

        pdf(sprintf("pics/CRISPRi_volcano_H19_%s.pdf", extra))
        makeVolcano(H19$P.Value, (H19$Versus1.logFC + H19$Versus2.logFC)/2,
                    which(rownames(H19)=="ENSG00000130600"), threshold=h19.threshold, ylim=ylim)
        dev.off()

        pdf(sprintf("pics/CRISPRi_het_volcano_H19_%s.pdf", extra))
        makeVolcano(H19.het$P.Value, (H19.het$Versus1.logFC + H19.het$Versus2.logFC)/2,
                    which(rownames(H19.het)=="ENSG00000130600"), threshold=h19.threshold, ylim=ylim)
        dev.off()
    }

    # Also comparing controls.
    LNAcon <- read.table(file.path(mode, "LNA_controlB_vs_trans.txt"), header=TRUE, stringsAsFactors=FALSE)
    CRISPRicon <- read.table(file.path(mode, "CRISPRi_clone2_vs_cells_I.txt"), header=TRUE, stringsAsFactors=FALSE)
    CRISPRi.hetcon <- read.table(file.path(mode, "CRISPRi_het_BFP_vs_cells.txt"), header=TRUE, stringsAsFactors=FALSE)
    flowthresh <- readRDS(sprintf("threshold_%s_control.rds", nextra))

    lna.sig <- LNAcon$P.Value <= flowthresh
    crispri.sig <- CRISPRicon$P.Value <= flowthresh
    crispri.het.sig <- CRISPRi.hetcon$P.Value <= flowthresh

    pdf(sprintf("pics/LNA_volcano_controlB_vs_trans_%s.pdf", nextra))
    makeVolcano(LNAcon$P.Value, LNAcon$logFC, lna.sig, threshold=flowthresh, xlim=c(-10, 10), chosen.pch=1.5, chosen.col="dodgerblue")
    abline(v=-3, col="black", lwd=3, lty=3)
    abline(v=3, col="black", lwd=3, lty=3)
    dev.off()

    pdf(sprintf("pics/CRISPRi_volcano_clone2_vs_cells_I.pdf", nextra))
    makeVolcano(CRISPRicon$P.Value, CRISPRicon$logFC, crispri.sig, threshold=flowthresh, xlim=c(-10, 10), chosen.pch=1.5, chosen.col="orange")
    abline(v=-3, col="black", lwd=3, lty=3)
    abline(v=3, col="black", lwd=3, lty=3)
    dev.off()

    pdf(sprintf("pics/CRISPRi_het_volcano_BFP_vs_cells.pdf", nextra))
    makeVolcano(CRISPRi.hetcon$P.Value, CRISPRi.hetcon$logFC, crispri.het.sig, threshold=flowthresh, xlim=c(-10, 10), chosen.pch=1.5, chosen.col="violet")
    abline(v=-3, col="black", lwd=3, lty=3)
    abline(v=3, col="black", lwd=3, lty=3)
    dev.off()


}

