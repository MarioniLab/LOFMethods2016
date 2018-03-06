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

# Generating MA plots.

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

pdf("pics/ma.pdf", width=10, height=5)
FUN("siRNA_TOG.txt", threshold, main="RNAi (ch-TOG)")
FUN("CRISPRi_TOG.txt", threshold, main="CRISPRi (ch-TOG)")
FUN("LNA_MALAT1.txt", threshold, main="LNA (MALAT1)")
FUN("CRISPRi_MALAT1.txt", threshold, main="CRISPRi (MALAT1)")
dev.off()
