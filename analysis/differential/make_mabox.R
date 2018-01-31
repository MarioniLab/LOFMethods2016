# This creates a series of smoothed scatter plots showing the distribution of log-fold changes 
# and mean abundances for upregulated, downregulated and other genes in each comparison.

FUN <- function(fname, threshold, main) {
    tab <- read.table(file.path("results_lfc", fname), header=TRUE, stringsAsFactors=FALSE, sep="\t")
    x <- tab$AveExpr
    y <- tab$logFC
    xlim <- range(x)
    ylim <- range(y)
    sig <- tab$P.Value < threshold
    up <- y > 0

    layout(rbind(1, c(2,3,4)), widths=c(3,2,2), heights=c(1,10))
    old <- par()$mar
    par(mar=c(0,0,0,0))
    plot(0,0,axes=FALSE, xlab="", ylab="", xaxt="n", yaxt="n", type="n")
    text(0,0, main, font=2, cex=1.6)

    par(mar=c(5.1, 4.1, 0, 2.1))
    smoothScatter(x[!sig], y[!sig], xlim=xlim, ylim=ylim, colramp = colorRampPalette(c("white", "black")),
                  nrpoints=0, xlab="Average log-expression", ylab="Log-fold change",
                  cex.axis=1.2, cex.lab=1.6)
    points(x[sig & up], y[sig & up], col=rgb(1,0,0,0.5), pch=16)
    points(x[sig & !up], y[sig & !up], col=rgb(0,0,1,0.5), pch=16)

    boxplot(list(`non-DE`=y[!sig],
                 Up=y[sig & up],
                 Down=y[sig & !up]), 
            border=c("black", "red", "blue"), ylab="Log-fold change", 
            cex.axis=1.2, cex.names=1.4, cex.lab=1.6)

    boxplot(list(`non-DE`=x[!sig],
                 Up=x[sig & up],
                 Down=x[sig & !up]), 
            border=c("black", "red", "blue"), ylab="Average log-expression", 
            cex.axis=1.2, cex.names=1.4, cex.lab=1.6)

    return(NULL)
}

#################
# Running through all the control comparisons used in the flow charts with LFC.

con_threshold <- readRDS("threshold_lfc_control.rds")

# CRISPRi (clonal):

pdf("pics/CRISPRi_MA_clone2_vs_cells_I.pdf", width=10, height=5)
FUN("CRISPRi_clone2_vs_cells_I.txt", con_threshold, main="CRISPRi (clonal), clone 2 vs cells")
dev.off()

pdf("pics/CRISPRi_MA_negguide1_vs_clone2.pdf", width=10, height=5)
FUN("CRISPRi_negguide1_vs_clone2.txt", con_threshold, main="CRISPRi (clonal), negative guide 1 vs cells")
dev.off()

pdf("pics/CRISPRi_MA_negguide2_vs_clone2_II.pdf", width=10, height=5)
FUN("CRISPRi_negguide2_vs_clone2_II.txt", con_threshold, main="CRISPRi (clonal), negative guide 2 vs cells")
dev.off()

pdf("pics/CRISPRi_MA_negguide1_vs_negguide2.pdf", width=10, height=5)
FUN("CRISPRi_negguide1_vs_negguide2.txt", con_threshold, main="CRISPRi (clonal), negative guide 1 vs negative guide 2")
dev.off()

# LNA:

pdf("pics/LNA_MA_trans_vs_cells.pdf", width=10, height=5)
FUN("LNA_trans_vs_cells.txt", con_threshold, main="LNA/RNAi, transfection control vs cells")
dev.off()

pdf("pics/LNA_MA_controlA_vs_trans.pdf", width=10, height=5)
FUN("LNA_controlA_vs_trans.txt", con_threshold, main="LNA, negative control A vs transfection control")
dev.off()

pdf("pics/LNA_MA_controlB_vs_trans.pdf", width=10, height=5)
FUN("LNA_controlB_vs_trans.txt", con_threshold, main="LNA, negative control B vs transfection control")
dev.off()

pdf("pics/LNA_MA_controlA_vs_controlB.pdf", width=10, height=5)
FUN("LNA_controlA_vs_controlB.txt", con_threshold, main="LNA, negative control A vs negative control B")
dev.off()

# RNAi:

pdf("pics/siRNA_MA_Ambion_vs_trans.pdf", width=10, height=5)
FUN("siRNA_Ambion_vs_trans.txt", con_threshold, main="LNA, Ambion control vs transfection control")
dev.off()

pdf("pics/siRNA_MA_Dharmacon_vs_trans.pdf", width=10, height=5)
FUN("siRNA_Dharmacon_vs_trans.txt", con_threshold, main="siRNA, Dharmacon control vs transfection control")
dev.off()

pdf("pics/siRNA_MA_Ambion_vs_Dharmacon.pdf", width=10, height=5)
FUN("siRNA_Ambion_vs_Dharmacon.txt", con_threshold, main="siRNA, Ambion control vs Dharmacon control")
dev.off()

# CRISPRi (heterogeneous):

pdf("pics/CRISPRi_het_MA_BFP_vs_cells.pdf", width=10, height=5)
FUN("CRISPRi_het_BFP_vs_cells.txt", con_threshold, main="CRISPRi (heterogeneous), BFP-expressing vs cells")
dev.off()

pdf("pics/CRISPRi_het_MA_Nc.1_vs_BFP.pdf", width=10, height=5)
FUN("CRISPRi_het_Nc.1_vs_BFP.txt", con_threshold, main="CRISPRi (heterogeneous), negative guide 1 vs untreated")
dev.off()

pdf("pics/CRISPRi_het_MA_Nc.2_vs_BFP.pdf", width=10, height=5)
FUN("CRISPRi_het_Nc.2_vs_BFP.txt", con_threshold, main="CRISPRi (heterogeneous), negative guide 2 vs untreated")
dev.off()

pdf("pics/CRISPRi_het_MA_Nc.1_vs_Nc.2.pdf", width=10, height=5)
FUN("CRISPRi_het_Nc.1_vs_Nc.2.txt", con_threshold, main="CRISPRi (heterogeneous), negative guide 1 vs negative guide 2")
dev.off()

#################
# Running through all the 289 comparisons with LFC.

threshold_289 <- readRDS("threshold_lfc_289.rds")

pdf("pics/CRISPRi_het_MA_289.pdf", width=10, height=5)
FUN("combined_CRISPRi_het_289.txt", threshold_289, main="CRISPRi (heterogeneous), effect of 289 knockdown")
dev.off()

pdf("pics/CRISPRi_MA_289.pdf", width=10, height=5)
FUN("combined_CRISPRi_289.txt", threshold_289, main="CRISPRi (clonal), effect of 289 knockdown")
dev.off()

pdf("pics/siRNA_MA_289.pdf", width=10, height=5)
FUN("combined_siRNA_289.txt", threshold_289, main="RNAi, effect of 289 knockdown")
dev.off()

pdf("pics/LNA_MA_289.pdf", width=10, height=5)
FUN("combined_LNA_289.txt", threshold_289, main="LNA, effect of 289 knockdown")
dev.off()
