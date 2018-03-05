# Setting up the threshold.
threshold <- readRDS("../../analysis/differential/threshold_lfc_control.rds")
dir.create("pics")

PROCESSOR <- function(pos, neg, target, largest=5, ..., col="red") {
    positive <- read.table(pos, header=TRUE)
    negative <- read.table(neg, header=TRUE)

    shared <- intersect(rownames(positive), rownames(negative))
    positive <- positive[shared,]
    negative <- negative[shared,]
    
    chosen <- which(positive$Symbol == target)
    is.P <- positive$P.Value <= threshold
    is.N <- negative$P.Value <= threshold
    
    X <- negative$logFC
    Y <- positive$logFC
    plot(X, Y, col="grey80", pch=16, xlim=c(-largest, largest), ylim=c(-largest, largest), ...)

    PnN <- is.P & !is.N
    points(X[PnN], Y[PnN], col=col, pch=16)
    NnP <- !is.P & is.N
    points(X[NnP], Y[NnP], col=col, bg="white", pch=21)
    NP <- is.P & is.N
    points(X[NP], Y[NP], col=col, pch=4)

    arrows(X[chosen] - 0.5, Y[chosen] + 0.5, X[chosen], Y[chosen], lwd=2, length=0.1)

    legend("topright", pch=c(16, 21, 4), legend=c(sprintf("DE in Y only (%i)", sum(PnN)),
                                                  sprintf("DE in X only (%i)", sum(NnP)),
                                                  sprintf("DE in both (%i)", sum(NP))), col=col)
    return(positive[PnN,])
}

# RNAi (TOG)
pdf("pics/TOG_RNAi.pdf")
tog.rnai <- PROCESSOR("results_lfc/siRNA_TOG.txt", "../../analysis/differential/results_lfc/siRNA_Ambion_vs_Dharmacon.txt",
                      xlab="Ambion versus Dharmacon", ylab="ch-TOG siRNA versus Dharmacon", 
                      target="CKAP5", main="RNAi (ch-TOG)", col="purple")
dev.off()

# CRISPRi (TOG)
pdf("pics/TOG_CRISPRi.pdf")
tog.crispri <- PROCESSOR("results_lfc/CRISPRi_TOG.txt", "../../analysis/differential/results_lfc/CRISPRi_het_Nc.1_vs_Nc.2.txt",
                         xlab="Negative guide 1 versus 2", ylab="ch-TOG guide versus negative guide 2", 
                         target="CKAP5", main="CRISPRi (ch-TOG)", col="brown")
dev.off()

# LNA (MALAT1)
pdf("pics/MALAT1_LNA.pdf")
out <- PROCESSOR("results_lfc/LNA_MALAT1.txt",
                 "../../analysis/differential/results_lfc/LNA_controlA_vs_controlB.txt",
                 xlab="Negative control A versus B", 
                 ylab="MALAT1 LNA versus control A", 
                 target="MALAT1", main="LNA (MALAT1)", col="dodgerblue")
dev.off()

# LNA (CRISPRi)
pdf("pics/MALAT1_CRISPRi.pdf")
out <- PROCESSOR("results_lfc/CRISPRi_MALAT1.txt",
                 "../../analysis/differential/results_lfc/CRISPRi_het_Nc.1_vs_Nc.2.txt",
                 xlab="Negative control A versus B", 
                 ylab="MALAT1 LNA versus control A", 
                 target="MALAT1", main="CRISPRi (MALAT1)", col="brown")
dev.off()

# Drugs.
pdf("pics/Monastrol.pdf")
out <- PROCESSOR("results_lfc/siRNA_TOG.txt",  "results_lfc/Monastrol.txt",
                 xlab="Monastrol vs unstreated", 
                 ylab="ch-TOG RNAi versus control Dharmacon", 
                 target="CKAP5", main="RNAi (ch-TOG)", col="brown")

out <- PROCESSOR("results_lfc/CRISPRi_TOG.txt",  "results_lfc/Monastrol.txt",
                 xlab="Monastrol vs unstreated", 
                 ylab="ch-TOG CRISPRi versus control guide", 
                 target="CKAP5", main="CRISPRi (ch-TOG)", col="brown")
dev.off()

pdf("pics/NOCO.pdf")
out <- PROCESSOR("results_lfc/siRNA_TOG.txt",  "results_lfc/NOCO.txt",
                 xlab="NOCO vs unstreated", 
                 ylab="ch-TOG RNAi versus control Dharmacon", 
                 target="CKAP5", main="RNAi (ch-TOG)", col="brown")

out <- PROCESSOR("results_lfc/CRISPRi_TOG.txt",  "results_lfc/NOCO.txt",
                 xlab="NOCO vs unstreated", 
                 ylab="ch-TOG CRISPRi versus control guide", 
                 target="CKAP5", main="CRISPRi (ch-TOG)", col="brown")
dev.off()
