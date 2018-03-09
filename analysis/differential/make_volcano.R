# Makes volcano plots, highlighting the targeted lncRNA in each case. 
# We use the average log-fold change and the maximum p-value in cases where we have two comparisons.

makeVolcano <- function(pvalues, logFCs, chosen, threshold, ylim=c(0, 30), xlim=c(-6, 6), sig.col="black", chosen.col="red", chosen.pch=2) {
    lp <- -log10(pvalues)
    sig <- pvalues <= threshold
    par(mar=c(5.1, 5.1, 4.1, 2.1))
    plot(logFCs[!sig], lp[!sig], pch=16, col="grey50", cex.axis=1.2, cex.lab=1.4, xlab=expression(Log[2]~"fold change"), 
         ylab=expression("-"*Log[10]~"P-value"), xlim=xlim, ylim=ylim, cex=1.5)
    points(logFCs[sig], lp[sig], pch=16, cex=1.5, col=sig.col)
    points(logFCs[chosen], lp[chosen], pch=16, col=chosen.col, cex=chosen.pch)
    abline(v=-0.5, col="red", lwd=2, lty=2)
    abline(v=0.5, col="red", lwd=2, lty=2)
    abline(h=-log10(threshold), col="red", lwd=2, lty=2)
}

for (mode in c("results_de", "results_lfc")) {
    if (mode=="results_de") {
        extra <- "de"
    } else {
        extra <- "lfc"
    }

    siRNA <- read.table(file.path(mode, "combined_siRNA_289.txt"), header=TRUE, stringsAsFactors=FALSE)
    LNA <- read.table(file.path(mode, "combined_LNA_289.txt"), header=TRUE, stringsAsFactors=FALSE)
    CRISPRi <- read.table(file.path(mode, "combined_CRISPRi_289.txt"), header=TRUE, stringsAsFactors=FALSE)
    CRISPRi.het <- read.table(file.path(mode, "combined_CRISPRi_het_289.txt"), header=TRUE, stringsAsFactors=FALSE)
 
    combined.p <- c(siRNA$P.Value, LNA$P.Value, CRISPRi$P.Value, CRISPRi.het$P.Value)
    is.sig <- p.adjust(combined.p, method="BH") <= 0.05
    threshold <- max(combined.p[is.sig])
    saveRDS(file=sprintf("threshold_%s_289.rds", extra), threshold)

    H19 <- read.table(file.path(mode, "combined_CRISPRi_H19.txt"), header=TRUE, stringsAsFactors=FALSE)
    H19.het <- read.table(file.path(mode, "combined_CRISPRi_het_H19.txt"), header=TRUE, stringsAsFactors=FALSE)
    h19.p <- c(H19$P.Value, H19.het$P.Value)
    h19.sig <- p.adjust(h19.p, method="BH") <= 0.05
    h19.threshold <- max(h19.p[h19.sig])
    saveRDS(file=sprintf("threshold_%s_H19.rds", extra), h19.threshold)

    ylim <- c(0, 30)
    pdf(sprintf("pics/siRNA_volcano_289_%s.pdf", extra))
    makeVolcano(siRNA$P.Value, siRNA$logFC, which(rownames(siRNA)=="ENSG00000234771"), threshold=threshold, ylim=ylim, sig.col="purple")
    dev.off()

    pdf(sprintf("pics/LNA_volcano_289_%s.pdf", extra))
    makeVolcano(LNA$P.Value, LNA$logFC, which(rownames(LNA)=="ENSG00000234771"), threshold=threshold, ylim=ylim, sig.col="dodgerblue")
    dev.off()

    pdf(sprintf("pics/CRISPRi_volcano_289_%s.pdf", extra))
    makeVolcano(CRISPRi$P.Value, CRISPRi$logFC, which(rownames(CRISPRi)=="ENSG00000234771"), threshold=threshold, ylim=ylim, sig.col="orange")
    dev.off()

    pdf(sprintf("pics/CRISPRi_het_volcano_289_%s.pdf", extra))
    makeVolcano(CRISPRi.het$P.Value, CRISPRi.het$logFC, which(rownames(CRISPRi.het)=="ENSG00000234771"), threshold=threshold, ylim=ylim, sig.col="brown")
    dev.off()

    pdf(sprintf("pics/CRISPRi_volcano_H19_%s.pdf", extra))
    makeVolcano(H19$P.Value, H19$logFC, which(rownames(H19)=="ENSG00000130600"), threshold=h19.threshold, ylim=ylim, sig.col="orange")
    dev.off()

    pdf(sprintf("pics/CRISPRi_het_volcano_H19_%s.pdf", extra))
    makeVolcano(H19.het$P.Value, H19.het$logFC, which(rownames(H19.het)=="ENSG00000130600"), threshold=h19.threshold, ylim=ylim, sig.col="orange")
    dev.off()

    # Also comparing controls.
    flowthresh <- readRDS(sprintf("threshold_%s_control.rds", extra))

    for (rnaif in c("siRNA_Ambion_vs_Dharmacon.txt", "siRNA_Ambion_vs_trans.txt", "siRNA_Dharmacon_vs_trans.txt")) {
        siRNAcon <- read.table(file.path(mode, rnaif), header=TRUE, stringsAsFactors=FALSE)
        rnai.sig <- siRNAcon$P.Value <= flowthresh

        pdf(sprintf("pics/siRNA_volcano_%s_%s.pdf", sub(".txt$", "", sub("^siRNA_", "", rnaif)), extra))
        makeVolcano(siRNAcon$P.Value, siRNAcon$logFC, NULL, threshold=flowthresh, sig.col="dodgerblue")
        dev.off()
    }

    for (lnaf in c("LNA_controlB_vs_trans.txt", "LNA_controlA_vs_controlB.txt", "LNA_controlA_vs_trans.txt", "LNA_trans_vs_cells.txt")) {
        LNAcon <- read.table(file.path(mode, lnaf), header=TRUE, stringsAsFactors=FALSE)
        lna.sig <- LNAcon$P.Value <= flowthresh

        pdf(sprintf("pics/LNA_volcano_%s_%s.pdf", sub(".txt$", "", sub("^LNA_", "", lnaf)), extra))
        makeVolcano(LNAcon$P.Value, LNAcon$logFC, NULL, threshold=flowthresh, sig.col="dodgerblue")
        dev.off()
    }

    for (crisprif in c("CRISPRi_het_Nc.1_vs_BFP.txt", "CRISPRi_het_Nc.1_vs_Nc.2.txt", "CRISPRi_het_Nc.2_vs_BFP.txt", "CRISPRi_het_BFP_vs_cells.txt")) {
        CRISPRicon <- read.table(file.path(mode, crisprif), header=TRUE, stringsAsFactors=FALSE)
        crispri.sig <- CRISPRicon$P.Value <= flowthresh

        pdf(sprintf("pics/CRISPRi_het_volcano_%s_%s.pdf", sub(".txt$", "", sub("^CRISPRi_het_", "", crisprif)), extra))
        makeVolcano(CRISPRicon$P.Value, CRISPRicon$logFC, NULL, threshold=flowthresh, sig.col="orange")
        dev.off()
    }

    for (crisprif in c("CRISPRi_clone2_vs_cells_I.txt", "CRISPRi_negguide1_vs_clone2.txt", "CRISPRi_negguide1_vs_negguide2.txt", "CRISPRi_negguide2_vs_clone2_II.txt")) { 
        CRISPRicon <- read.table(file.path(mode, crisprif), header=TRUE, stringsAsFactors=FALSE)
        crispri.sig <- CRISPRicon$P.Value <= flowthresh

        pdf(sprintf("pics/CRISPRi_volcano_%s_%s.pdf", sub(".txt$", "", sub("^CRISPRi_", "", crisprif)), extra))
        makeVolcano(CRISPRicon$P.Value, CRISPRicon$logFC, NULL, threshold=flowthresh, sig.col="orange")
        dev.off()
    }
}

