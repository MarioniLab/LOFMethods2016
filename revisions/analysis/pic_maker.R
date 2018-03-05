all.results <- list(TOG.RNAi=read.table("results_lfc/siRNA_TOG.txt", header=TRUE),
                    TOG.CRISPRi=read.table("results_lfc/CRISPRi_TOG.txt", header=TRUE),
                    MALAT1.LNA=read.table("results_lfc/LNA_MALAT1.txt", header=TRUE),
                    MALAT1.CRISPRi=read.table("results_lfc/CRISPRi_MALAT1.txt", header=TRUE))

###############################

all.p <- unlist(lapply(all.results, FUN="[[", i="P.Value"))
adj.p <- p.adjust(all.p, method="BH")
threshold <- max(all.p[adj.p <= 0.05])

# Setting up a function to compare between the different technologies.
dir.create("pics")
PROCESSOR <- function(positive, negative, target, pos.lab="Y", neg.lab="X", largest=5, ..., col="red", arrow.col="black") {
    shared <- intersect(rownames(positive), rownames(negative))
    positive <- positive[shared,]
    negative <- negative[shared,]
    
    chosen <- which(positive$Symbol == target)
    is.P <- positive$P.Value <= threshold
    is.N <- negative$P.Value <= threshold
    
    X <- negative$logFC
    Y <- positive$logFC
    none <- !is.P & !is.N
    plot(X[none], Y[none], col="grey80", pch=16, 
         xlim=c(-largest, largest), 
         ylim=c(-largest, largest), ...)

    PnN <- is.P & !is.N
    points(X[PnN], Y[PnN], col=col, pch=16)
    NnP <- !is.P & is.N
    points(X[NnP], Y[NnP], col=col, bg="white", pch=21)
    NP <- is.P & is.N
    points(X[NP], Y[NP], col=col, pch=4)

    arrows(X[chosen] - 0.5, Y[chosen], X[chosen] - 0.1, Y[chosen], lwd=2, length=0.1, col=arrow.col)

    legend("topright", pch=c(16, 21, 4), legend=c(sprintf("DE in %s only (%i)", pos.lab, sum(PnN)),
                                                  sprintf("DE in %s only (%i)", neg.lab, sum(NnP)),
                                                  sprintf("DE in both (%i)", sum(NP))), col=col)
    return(NULL)
}

# Comparisons between methods, within the same gene.
pdf("pics/TOG_methods.pdf")
PROCESSOR(all.results$TOG.CRISPRi, all.results$TOG.RNAi,
          xlab="siRNA versus Dharmacon", ylab="CRISPRi vs negative guide", 
          target="CKAP5", main="ch-TOG", pos.lab="CRISPRi", neg.lab="RNAi", col="black", arrow.col="red")
dev.off()

pdf("pics/MALAT1_methods.pdf")
PROCESSOR(all.results$MALAT1.CRISPRi, all.results$MALAT1.LNA,
          xlab="LNA versus negative control A", ylab="CRISPRi vs negative guide", 
          target="MALAT1", main="MALAT1", pos.lab="CRISPRi", neg.lab="LNA", col="black", arrow.col="red")
dev.off()

# Drugs.
pdf("pics/Monastrol.pdf")
monastrol <- read.table("results_lfc/Monastrol.txt", header=TRUE)
PROCESSOR(all.results$TOG.RNAi, monastrol,
          xlab="Monastrol vs unstreated", 
          ylab="ch-TOG RNAi versus control Dharmacon", 
          target="CKAP5", main="RNAi (ch-TOG)", pos.lab="CRISPRi", neg.lab="treated", col="brown")

PROCESSOR(all.results$TOG.CRISPRi, monastrol,
          xlab="Monastrol vs unstreated", 
          ylab="ch-TOG CRISPRi versus control guide", 
          target="CKAP5", main="CRISPRi (ch-TOG)", pos.lab="CRISPRi", neg.lab="treated", col="brown")
dev.off()

pdf("pics/NOCO.pdf")
monastrol <- read.table("results_lfc/NOCO.txt", header=TRUE)
PROCESSOR(all.results$TOG.RNAi, monastrol,
          xlab="NOCO vs unstreated", 
          ylab="ch-TOG RNAi versus control Dharmacon", 
          target="CKAP5", main="RNAi (ch-TOG)", pos.lab="CRISPRi", neg.lab="treated", col="brown")

PROCESSOR(all.results$TOG.CRISPRi, monastrol,
          xlab="NOCO vs unstreated", 
          ylab="ch-TOG CRISPRi versus control guide", 
          target="CKAP5", main="CRISPRi (ch-TOG)", pos.lab="CRISPRi", neg.lab="treated", col="brown")
dev.off()

###############################

# Comparisons between positive and negative comparisons.

ref.results <- list(RNAi=read.table("../../analysis/differential/results_lfc/siRNA_Ambion_vs_Dharmacon.txt", header=TRUE),
                    LNA=read.table("../../analysis/differential/results_lfc/LNA_controlA_vs_controlB.txt", header=TRUE),
                    CRISPRi=read.table("../../analysis/differential/results_lfc/CRISPRi_het_Nc.1_vs_Nc.2.txt", header=TRUE))

GENERATE_CURVE <- function(pos, neg) {
    keep <- intersect(rownames(pos), rownames(neg))
    findInterval(sort(neg$P.Value[rownames(neg) %in% keep]), sort(pos$P.Value[rownames(pos) %in% keep]))
}

pdf("roc.pdf")
rnai <- GENERATE_CURVE(all.results$TOG.RNAi, ref.results$RNAi)
plot(seq_along(rnai), rnai, type="l", lwd=2, col="purple", xlim=c(0, 100), ylim=c(0, 1000),
     ylab="Number of DE genes upon knockdown", xlab="Number of DE genes between negative controls")
lna <- GENERATE_CURVE(all.results$MALAT1.LNA, ref.results$LNA)
lines(seq_along(lna), lna, lwd=2, col="dodgerblue")
crispri <- GENERATE_CURVE(all.results$TOG.CRISPRi, ref.results$CRISPRi)
lines(seq_along(crispri), crispri, lwd=2, col="brown")
crispri <- GENERATE_CURVE(all.results$MALAT1.CRISPRi, ref.results$CRISPRi)
lines(seq_along(crispri), crispri, lwd=2, col="brown", lty=2)

legend("topleft", col=c("purple", "dodgerblue", "brown", "brown"),
       lwd=2, lty=c(1,1,1,2), legend=c("RNAi (ch-TOG)", "LNA (MALAT1)", "CRISPRi (ch-TOG)", "CRISPRi (MALAT1)"))
dev.off()
