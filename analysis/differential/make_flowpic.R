#############################################################################
# Get number of DE genes.

getNumDE <- function(dir, fnames) {
    collected.p <- collected.i <- list()
    for (f in seq_along(fnames)) {
        res <- read.table(file.path(dir, fnames[[f]]), header=TRUE, sep="\t", stringsAsFactors=FALSE)
        collected.p[[f]] <- res$P.Value
        collected.i[[f]] <- rep(f, nrow(res))
    }
    collected.p <- unlist(collected.p)
    is.sig <- p.adjust(collected.p, method="BH") <= 0.05
    collected.i <- unlist(collected.i)
    out <- tabulate(collected.i[is.sig], nbins=length(fnames))
    names(out) <- names(fnames)
    return(list(numbers=out, threshold=max(collected.p[is.sig])))
}

textAtMid <- function(x0, y0, x1, y1, ..., adj=0) {
    m <- (y1-y0)/(x1-x0)
    b <- y0 - x0*m
    mid <- (x0 + x1)/2
    text(mid, mid*m + b + adj, ...)
}

for (mode in c("results_de", "results_lfc")) {
    if (mode=="results_de") {
        extra <- "de"
        scaling <- 1000
    } else {
        extra <- "lfc"
        scaling <- 100
    }
    
    de.out <- getNumDE(mode, list(
                               Clone2vCells="CRISPRi_clone2_vs_cells_I.txt",
                               Guide1vClone2="CRISPRi_negguide1_vs_clone2.txt",
                               Guide2vClone2="CRISPRi_negguide2_vs_clone2_II.txt",
                               Guide1vGuide2="CRISPRi_negguide1_vs_negguide2.txt",
#                               G289.1vGuide2="CRISPRi_289.1_vs_negguide2.txt",
#                               G289.9vGuide2="CRISPRi_289.9_vs_negguide2.txt",
                               TransvCells="LNA_trans_vs_cells.txt",
                               ConAvTrans="LNA_controlA_vs_trans.txt",
                               ConBvTrans="LNA_controlB_vs_trans.txt",
                               ConAvConB="LNA_controlA_vs_controlB.txt",
#                               LNA289.2vConA="LNA_289.2_vs_controlA.txt",
#                               LNA289.2vConB="LNA_289.2_vs_controlB.txt",
                               AmbionvTrans="siRNA_Ambion_vs_trans.txt",
                               DharmaconvTrans="siRNA_Dharmacon_vs_trans.txt",
                               AmbionvDharmacon="siRNA_Ambion_vs_Dharmacon.txt",
#                               si289vAmbion="siRNA_289_vs_Ambion.txt",
#                               si289vDharmacon="siRNA_289_vs_Dharmacon.txt",
                               HetvCells="CRISPRi_het_BFP_vs_cells.txt",
                               Guide1vHet="CRISPRi_het_Nc.1_vs_BFP.txt",
                               Guide2vHet="CRISPRi_het_Nc.2_vs_BFP.txt",
                               Het.Guide1vGuide2="CRISPRi_het_Nc.1_vs_Nc.2.txt"
                               ))
    all.nDE <- de.out$numbers
    saveRDS(de.out$threshold, paste0("threshold_", extra, "_control.rds"))

    ############################################################################
    # Generates a flow diagram for CRISPRi (clonal).

    pdf(sprintf("pics/CRISPRi_flow_%s.pdf", extra), width=20, height=8)
    par(mar=c(0,0,0,0))
    plot(0, 0, type="n", axes=FALSE, ylab="", xlab="", xlim=c(0, 100), ylim=c(-15, 30))

    last.pos <- 0
    mid.height <- 10
    box.sides <- 10
    shift <- 20

    rect(last.pos, mid.height, last.pos + box.sides, mid.height + box.sides, col="orange")
    text(last.pos + box.sides/2, mid.height + box.sides/2, "Cells", cex=1.4)

    last.pos <- last.pos + shift
    rect(last.pos, mid.height, last.pos + box.sides, mid.height + box.sides, col="orange")
    text(last.pos + box.sides/2, mid.height + box.sides/2, "dCas9-KRAB\nclone 2", cex=1.4)

    last.pos <- last.pos + shift
    upper.height <- mid.height + 10
    rect(last.pos, upper.height, last.pos + box.sides, upper.height + box.sides, col="orange")
    text(last.pos + box.sides/2, upper.height + box.sides/2, "Negative\nguide 1", cex=1.4)
    lower.height <- mid.height - 10
    rect(last.pos, lower.height, last.pos + box.sides, lower.height + box.sides, col="orange")
    text(last.pos + box.sides/2, lower.height + box.sides/2, "Negative\nguide 2", cex=1.4)

#    last.pos <- last.pos + shift
#    lower.heightA <- lower.height + 7
#    rect(last.pos, lower.heightA, last.pos + box.sides, lower.heightA + box.sides, col="orange")
#    text(last.pos + box.sides/2, lower.heightA + box.sides/2, "289 guide 1", cex=1.4)
#    lower.heightB <- lower.height - 7
#    rect(last.pos, lower.heightB, last.pos + box.sides, lower.heightB + box.sides, col="orange")
#    text(last.pos + box.sides/2, lower.heightB + box.sides/2, "289 guide 9", cex=1.4)

    # Adding connecting lines.
    last.pos <- 0
    nDE <- all.nDE[["Clone2vCells"]]
    lwidth <- max(0.1, nDE/scaling)
    rect(last.pos + box.sides, mid.height + box.sides/2 - lwidth/2, last.pos + shift, mid.height + box.sides/2 + lwidth/2, col="black", border=NA) # DE between clone and cells.
    text(last.pos + (shift+box.sides)/2, mid.height + box.sides/2 - lwidth/2, pos=1, nDE, cex=1.4)

    last.pos <- last.pos + shift
    nDE <- all.nDE[["Guide1vClone2"]]
    lwidth <- max(0.1, nDE/scaling)
    polygon(rep(c(last.pos + box.sides, last.pos + shift), each=2),
            c(mid.height + box.sides, mid.height + box.sides - lwidth, upper.height + box.sides/2 - lwidth/2, upper.height + box.sides/2 + lwidth/2),
            col="black", border=NA) 
    textAtMid(last.pos + box.sides, mid.height + box.sides, last.pos + shift, upper.height + box.sides/2 + lwidth/2, nDE, cex=1.4, pos=3, adj=0.5)

    nDE <- all.nDE[["Guide2vClone2"]]
    lwidth <- max(0.1, nDE/scaling)
    polygon(rep(c(last.pos + box.sides, last.pos + shift), each=2),
            c(mid.height + lwidth, mid.height, lower.height + box.sides/2 - lwidth/2, lower.height + box.sides/2 + lwidth/2),
            col="black", border=NA) 
    textAtMid(last.pos + box.sides, mid.height, last.pos + shift, lower.height + box.sides/2 - lwidth/2, nDE, cex=1.4, pos=1)

    last.pos <- last.pos + shift
    nDE <- all.nDE[["Guide1vGuide2"]]
    lwidth <- max(0.1, nDE/scaling)
    rect(last.pos + box.sides/2 - lwidth/2, lower.height + box.sides, last.pos + box.sides/2 + lwidth/2, upper.height, col="black")
    text(last.pos + box.sides/2 + lwidth/2, (lower.height + box.sides + upper.height)/2, pos=4, paste0(nDE, "*"), cex=1.4)

#    nDE <- all.nDE[["G289.1vGuide2"]]
#    lwidth <- max(0.1, nDE/scaling)
#    polygon(rep(c(last.pos + box.sides, last.pos + shift), each=2),
#            c(lower.height + box.sides, lower.height + box.sides - lwidth, lower.heightA + box.sides/2 - lwidth/2, lower.heightA + box.sides/2 + lwidth/2),
#            col="black", border=NA) # DE between guide 1 and clone.
#    textAtMid(last.pos + box.sides, lower.height + box.sides, last.pos + shift, lower.heightA + box.sides/2 + lwidth/2, nDE, cex=1.4, pos=3)
#
#    nDE <- all.nDE[["G289.9vGuide2"]]
#    lwidth <- max(0.1, nDE/scaling)
#    polygon(rep(c(last.pos + box.sides, last.pos + shift), each=2),
#            c(lower.height + lwidth, lower.height, lower.heightB + box.sides/2 - lwidth/2, lower.heightB + box.sides/2 + lwidth/2),
#            col="black", border=NA) # DE between guide 1 and clone.
#    textAtMid(last.pos + box.sides, lower.height, last.pos + shift, lower.heightB + box.sides/2 - lwidth/2, nDE, cex=1.4, pos=1)

    # Adding explanation.
    last.pos <- 0
    extras <- 1
    line.pos <- -8

    lines(rep(c(last.pos + box.sides-extras, last.pos+shift+extras), each=2),
          c(line.pos, line.pos - extras, line.pos - extras, line.pos), lwd=1.5)
    text(last.pos+(box.sides+shift)/2, line.pos-extras-1, pos=1, "dCas9-KRAB +\ncloning effects", cex=1.4)

    last.pos <- last.pos + shift
    lines(rep(c(last.pos + box.sides-extras, last.pos+shift+extras), each=2),
          c(line.pos, line.pos - extras, line.pos - extras, line.pos), lwd=1.5)
    text(last.pos+(box.sides+shift)/2, line.pos-extras-1, pos=1, "Transduction +\noff-target effects", cex=1.4)

#    last.pos <- last.pos + shift
#    lines(rep(c(last.pos + box.sides-extras, last.pos+shift+extras), each=2),
#          c(line.pos, line.pos - extras, line.pos - extras, line.pos), lwd=1.5)
#    text(last.pos+(box.sides+shift)/2, line.pos-extras-1, pos=1, "Knockdown +\noff-target effects", cex=1.4)

    dev.off()

    ############################################################################
    # Generates a flow diagram for CRISPRi (heterogeneous).

    pdf(sprintf("pics/CRISPRi_het_flow_%s.pdf", extra), width=20, height=8)
    par(mar=c(0,0,0,0))
    plot(0, 0, type="n", axes=FALSE, ylab="", xlab="", xlim=c(0, 100), ylim=c(-15, 30))

    last.pos <- 0
    mid.height <- 10
    box.sides <- 10
    shift <- 20

    rect(last.pos, mid.height, last.pos + box.sides, mid.height + box.sides, col="violet")
    text(last.pos + box.sides/2, mid.height + box.sides/2, "Cells", cex=1.4)

    last.pos <- last.pos + shift
    rect(last.pos, mid.height, last.pos + box.sides, mid.height + box.sides, col="violet")
    text(last.pos + box.sides/2, mid.height + box.sides/2, "dCas9-KRAB\nBFP-positive", cex=1.4)

    last.pos <- last.pos + shift
    upper.height <- mid.height + 10
    rect(last.pos, upper.height, last.pos + box.sides, upper.height + box.sides, col="violet")
    text(last.pos + box.sides/2, upper.height + box.sides/2, "Negative\nguide 1", cex=1.4)
    lower.height <- mid.height - 10
    rect(last.pos, lower.height, last.pos + box.sides, lower.height + box.sides, col="violet")
    text(last.pos + box.sides/2, lower.height + box.sides/2, "Negative\nguide 2", cex=1.4)

#    last.pos <- last.pos + shift
#    lower.heightA <- lower.height + 7
#    rect(last.pos, lower.heightA, last.pos + box.sides, lower.heightA + box.sides, col="violet")
#    text(last.pos + box.sides/2, lower.heightA + box.sides/2, "289 guide 1", cex=1.4)
#    lower.heightB <- lower.height - 7
#    rect(last.pos, lower.heightB, last.pos + box.sides, lower.heightB + box.sides, col="violet")
#    text(last.pos + box.sides/2, lower.heightB + box.sides/2, "289 guide 9", cex=1.4)

    # Adding connecting lines.
    last.pos <- 0
    nDE <- all.nDE[["HetvCells"]]
    lwidth <- max(0.1, nDE/scaling)
    rect(last.pos + box.sides, mid.height + box.sides/2 - lwidth/2, last.pos + shift, mid.height + box.sides/2 + lwidth/2, col="black", border=NA) # DE between clone and cells.
    text(last.pos + (shift+box.sides)/2, mid.height + box.sides/2 - lwidth/2, pos=1, nDE, cex=1.4)

    last.pos <- last.pos + shift
    nDE <- all.nDE[["Guide1vHet"]]
    lwidth <- max(0.1, nDE/scaling)
    polygon(rep(c(last.pos + box.sides, last.pos + shift), each=2),
            c(mid.height + box.sides, mid.height + box.sides - lwidth, upper.height + box.sides/2 - lwidth/2, upper.height + box.sides/2 + lwidth/2),
            col="black", border=NA) 
    textAtMid(last.pos + box.sides, mid.height + box.sides, last.pos + shift, upper.height + box.sides/2 + lwidth/2, nDE, cex=1.4, pos=3, adj=0.5)

    nDE <- all.nDE[["Guide2vHet"]]
    lwidth <- max(0.1, nDE/scaling)
    polygon(rep(c(last.pos + box.sides, last.pos + shift), each=2),
            c(mid.height + lwidth, mid.height, lower.height + box.sides/2 - lwidth/2, lower.height + box.sides/2 + lwidth/2),
            col="black", border=NA) 
    textAtMid(last.pos + box.sides, mid.height, last.pos + shift, lower.height + box.sides/2 - lwidth/2, nDE, cex=1.4, pos=1)

    last.pos <- last.pos + shift
    nDE <- all.nDE[["Guide1vGuide2"]]
    lwidth <- max(0.1, nDE/scaling)
    rect(last.pos + box.sides/2 - lwidth/2, lower.height + box.sides, last.pos + box.sides/2 + lwidth/2, upper.height, col="black")
    text(last.pos + box.sides/2 + lwidth/2, (lower.height + box.sides + upper.height)/2, pos=4, paste0(nDE, "*"), cex=1.4)

#    nDE <- all.nDE[["G289.1vGuide2"]]
#    lwidth <- max(0.1, nDE/scaling)
#    polygon(rep(c(last.pos + box.sides, last.pos + shift), each=2),
#            c(lower.height + box.sides, lower.height + box.sides - lwidth, lower.heightA + box.sides/2 - lwidth/2, lower.heightA + box.sides/2 + lwidth/2),
#            col="black", border=NA) # DE between guide 1 and clone.
#    textAtMid(last.pos + box.sides, lower.height + box.sides, last.pos + shift, lower.heightA + box.sides/2 + lwidth/2, nDE, cex=1.4, pos=3)
#
#    nDE <- all.nDE[["G289.9vGuide2"]]
#    lwidth <- max(0.1, nDE/scaling)
#    polygon(rep(c(last.pos + box.sides, last.pos + shift), each=2),
#            c(lower.height + lwidth, lower.height, lower.heightB + box.sides/2 - lwidth/2, lower.heightB + box.sides/2 + lwidth/2),
#            col="black", border=NA) # DE between guide 1 and clone.
#    textAtMid(last.pos + box.sides, lower.height, last.pos + shift, lower.heightB + box.sides/2 - lwidth/2, nDE, cex=1.4, pos=1)

    # Adding explanation.
    last.pos <- 0
    extras <- 1
    line.pos <- -8

    lines(rep(c(last.pos + box.sides-extras, last.pos+shift+extras), each=2),
          c(line.pos, line.pos - extras, line.pos - extras, line.pos), lwd=1.5)
    text(last.pos+(box.sides+shift)/2, line.pos-extras-1, pos=1, "dCas9-KRAB\neffects", cex=1.4)

    last.pos <- last.pos + shift
    lines(rep(c(last.pos + box.sides-extras, last.pos+shift+extras), each=2),
          c(line.pos, line.pos - extras, line.pos - extras, line.pos), lwd=1.5)
    text(last.pos+(box.sides+shift)/2, line.pos-extras-1, pos=1, "Transduction +\noff-target effects", cex=1.4)

#    last.pos <- last.pos + shift
#    lines(rep(c(last.pos + box.sides-extras, last.pos+shift+extras), each=2),
#          c(line.pos, line.pos - extras, line.pos - extras, line.pos), lwd=1.5)
#    text(last.pos+(box.sides+shift)/2, line.pos-extras-1, pos=1, "Knockdown +\noff-target effects", cex=1.4)

    dev.off()

    #############################################################################
    # Generates a flow diagram for LNA

    pdf(sprintf("pics/LNA_flow_%s.pdf", extra), width=20, height=8)
    par(mar=c(0,0,0,0))
    plot(0, 0, type="n", axes=FALSE, ylab="", xlab="", xlim=c(0, 100), ylim=c(-15, 30))

    last.pos <- 0
    mid.height <- 10
    box.sides <- 10
    shift <- 20

    rect(last.pos, mid.height, last.pos + box.sides, mid.height + box.sides, col="lightblue")
    text(last.pos + box.sides/2, mid.height + box.sides/2, "Cells", cex=1.4)

    last.pos <- last.pos + shift
    rect(last.pos, mid.height, last.pos + box.sides, mid.height + box.sides, col="lightblue")
    text(last.pos + box.sides/2, mid.height + box.sides/2, "Transfection\ncontrol", cex=1.4)

    last.pos <- last.pos + shift
    upper.height <- mid.height + 10
    rect(last.pos, upper.height, last.pos + box.sides, upper.height + box.sides, col="lightblue")
    text(last.pos + box.sides/2, upper.height + box.sides/2, "Negative\ncontrol A", cex=1.4)
    lower.height <- mid.height - 10
    rect(last.pos, lower.height, last.pos + box.sides, lower.height + box.sides, col="lightblue")
    text(last.pos + box.sides/2, lower.height + box.sides/2, "Negative\ncontrol B", cex=1.4)

#    last.pos <- last.pos + shift
#    rect(last.pos, mid.height, last.pos + box.sides, mid.height + box.sides, col="lightblue")
#    text(last.pos + box.sides/2, mid.height + box.sides/2, "289 LNA 2", cex=1.4)

    # Adding connecting lines.
#    all.nDE <- getNumDE(mode, list(LNA289.2vConB="LNA_289.2_vs_controlB.txt"))

    last.pos <- 0
    nDE <- all.nDE[["TransvCells"]]
    lwidth <- max(0.1, nDE/scaling)
    rect(last.pos + box.sides, mid.height + box.sides/2 - lwidth/2, last.pos + shift, mid.height + box.sides/2 + lwidth/2, col="black", border=NA) 
    text(last.pos + (shift+box.sides)/2, mid.height + box.sides/2 - lwidth/2, pos=1, nDE, cex=1.4)

    last.pos <- last.pos + shift
    nDE <- all.nDE[["ConAvTrans"]]
    lwidth <- max(0.1, nDE/scaling)
    polygon(rep(c(last.pos + box.sides, last.pos + shift), each=2),
            c(mid.height + box.sides, mid.height + box.sides - lwidth, upper.height + box.sides/2 - lwidth/2, upper.height + box.sides/2 + lwidth/2),
            col="black", border=NA) 
    textAtMid(last.pos + box.sides, mid.height + box.sides, last.pos + shift, upper.height + box.sides/2 + lwidth/2, nDE, cex=1.4, pos=3)

    nDE <- all.nDE[["ConBvTrans"]]
    lwidth <- max(0.1, nDE/scaling)
    polygon(rep(c(last.pos + box.sides, last.pos + shift), each=2),
            c(mid.height + lwidth, mid.height, lower.height + box.sides/2 - lwidth/2, lower.height + box.sides/2 + lwidth/2),
            col="black", border=NA) 
    textAtMid(last.pos + box.sides, mid.height, last.pos + shift, lower.height + box.sides/2 - lwidth/2, nDE, cex=1.4, pos=1, adj=-0.5)

    last.pos <- last.pos + shift
    nDE <- all.nDE[["ConAvConB"]]
    lwidth <- max(0.1, nDE/scaling)
    rect(last.pos + box.sides/2 - lwidth/2, lower.height + box.sides, last.pos + box.sides/2 + lwidth/2, upper.height, col="black")
    text(last.pos + box.sides/2 + lwidth/2, (lower.height + box.sides + upper.height)/2, pos=4, paste0(nDE, "*"), cex=1.4)

#    nDE <- all.nDE[["LNA289.2vConA"]]
#    lwidth <- max(0.1, nDE/scaling)
#    polygon(rep(c(last.pos + box.sides, last.pos + shift), each=2),
#            c(upper.height + box.sides/2 + lwidth/2, upper.height + box.sides/2 - lwidth/2, mid.height + box.sides - lwidth, mid.height+ box.sides),
#            col="black")
#    textAtMid(last.pos + box.sides, upper.height + box.sides/2 + lwidth/2, last.pos + shift, mid.height + box.sides, nDE, cex=1.4, pos=3, adj=0.8)
#
#    nDE <- all.nDE[["LNA289.2vConB"]]
#    lwidth <- max(0.1, nDE/scaling)
#    polygon(rep(c(last.pos + box.sides, last.pos + shift), each=2),
#            c(lower.height + box.sides/2 + lwidth/2, lower.height + box.sides/2 - lwidth/2, mid.height, mid.height+lwidth),
#            col="black")
#    textAtMid(last.pos + box.sides, lower.height + box.sides/2 - lwidth/2, last.pos + shift, mid.height, nDE, cex=1.4, pos=1, adj=-0.8)

    # Adding explanation.
    last.pos <- 0
    extras <- 1
    line.pos <- -2

    lines(rep(c(last.pos + box.sides-extras, last.pos+shift+extras), each=2),
          c(line.pos, line.pos - extras, line.pos - extras, line.pos), lwd=1.5)
    text(last.pos+(box.sides+shift)/2, line.pos-extras-1, pos=1, "Transfection effect", cex=1.4)

    last.pos <- last.pos + shift
    lines(rep(c(last.pos + box.sides-extras, last.pos+shift+extras), each=2),
          c(line.pos, line.pos - extras, line.pos - extras, line.pos), lwd=1.5)
    text(last.pos+(box.sides+shift)/2, line.pos-extras-1, pos=1, "Off-target effect", cex=1.4)

#    last.pos <- last.pos + shift
#    lines(rep(c(last.pos + box.sides-extras, last.pos+shift+extras), each=2),
#          c(line.pos, line.pos - extras, line.pos - extras, line.pos), lwd=1.5)
#    text(last.pos+(box.sides+shift)/2, line.pos-extras-1, pos=1, "Knockdown +\noff-target effects", cex=1.4)

    dev.off()

    #############################################################################
    # Generates a flow diagram for siRNA

    pdf(sprintf("pics/siRNA_flow_%s.pdf", extra), width=20, height=8)
    par(mar=c(0,0,0,0))
    plot(0, 0, type="n", axes=FALSE, ylab="", xlab="", xlim=c(0, 100), ylim=c(-15, 30))

    last.pos <- 0
    mid.height <- 10
    box.sides <- 10
    shift <- 20

    rect(last.pos, mid.height, last.pos + box.sides, mid.height + box.sides, col="salmon")
    text(last.pos + box.sides/2, mid.height + box.sides/2, "Cells", cex=1.4)

    last.pos <- last.pos + shift
    rect(last.pos, mid.height, last.pos + box.sides, mid.height + box.sides, col="salmon")
    text(last.pos + box.sides/2, mid.height + box.sides/2, "Transfection\ncontrol", cex=1.4)

    last.pos <- last.pos + shift
    upper.height <- mid.height + 10
    rect(last.pos, upper.height, last.pos + box.sides, upper.height + box.sides, col="salmon")
    text(last.pos + box.sides/2, upper.height + box.sides/2, "Ambion\ncontrol", cex=1.4)
    lower.height <- mid.height - 10
    rect(last.pos, lower.height, last.pos + box.sides, lower.height + box.sides, col="salmon")
    text(last.pos + box.sides/2, lower.height + box.sides/2, "Dharmacon\ncontrol", cex=1.4)

#    last.pos <- last.pos + shift
#    rect(last.pos, mid.height, last.pos + box.sides, mid.height + box.sides, col="salmon")
#    text(last.pos + box.sides/2, mid.height + box.sides/2, "289 siRNA", cex=1.4)

    # Adding connecting lines.
#    all.nDE <- getNumDE(mode, list(AmbionvCells="siRNA_Ambion_vs_cells.txt",
#                                   DharmaconvCells="siRNA_Dharmacon_vs_cells.txt",
#                                   AmbionvDharmacon="siRNA_Ambion_vs_Dharmacon.txt",
#                                   si289vAmbion="siRNA_289_vs_Ambion.txt",
#                                   si289vDharmacon="siRNA_289_vs_Dharmacon.txt"))

    last.pos <- 0
    nDE <- all.nDE[["TransvCells"]]
    lwidth <- max(0.1, nDE/scaling)
    rect(last.pos + box.sides, mid.height + box.sides/2 - lwidth/2, last.pos + shift, mid.height + box.sides/2 + lwidth/2, col="black", border=NA) 
    text(last.pos + (shift+box.sides)/2, mid.height + box.sides/2 - lwidth/2, pos=1, nDE, cex=1.4)

    last.pos <- last.pos + shift
    nDE <- all.nDE[["AmbionvTrans"]]
    lwidth <- max(0.1, nDE/scaling)
    polygon(rep(c(last.pos + box.sides, last.pos + shift), each=2),
            c(mid.height + box.sides, mid.height + box.sides - lwidth, upper.height + box.sides/2 - lwidth/2, upper.height + box.sides/2 + lwidth/2),
            col="black", border=NA) 
    textAtMid(last.pos + box.sides, mid.height + box.sides, last.pos + shift, upper.height + box.sides/2 + lwidth/2, nDE, cex=1.4, pos=3, adj=0.5)

    nDE <- all.nDE[["DharmaconvTrans"]]
    lwidth <- max(0.1, nDE/scaling)
    polygon(rep(c(last.pos + box.sides, last.pos + shift), each=2),
            c(mid.height + lwidth, mid.height, lower.height + box.sides/2 - lwidth/2, lower.height + box.sides/2 + lwidth/2),
            col="black", border=NA) 
    textAtMid(last.pos + box.sides, mid.height, last.pos + shift, lower.height + box.sides/2 - lwidth/2, nDE, cex=1.4, pos=1, adj=-0.5)

    last.pos <- last.pos + shift
    nDE <- all.nDE[["AmbionvDharmacon"]]
    lwidth <- max(0.1, nDE/scaling)
    rect(last.pos + box.sides/2 - lwidth/2, lower.height + box.sides, last.pos + box.sides/2 + lwidth/2, upper.height, col="black")
    text(last.pos + box.sides/2 + lwidth/2, (lower.height + box.sides + upper.height)/2, pos=4, paste0(nDE, "*"), cex=1.4)

#    nDE <- all.nDE[["si289vAmbion"]]
#    lwidth <- max(0.1, nDE/scaling)
#    polygon(rep(c(last.pos + shift, last.pos + box.sides), each=2),
#            c(mid.height + box.sides, mid.height + box.sides - lwidth, upper.height + box.sides/2 - lwidth/2, upper.height + box.sides/2 + lwidth/2),
#            col="black", border=NA) 
#    textAtMid(last.pos + box.sides, upper.height + box.sides/2 + lwidth/2, last.pos + shift, mid.height + box.sides, nDE, cex=1.4, pos=3, adj=0.8)
#
#    nDE <- all.nDE[["si289vDharmacon"]]
#    lwidth <- max(0.1, nDE/scaling)
#    polygon(rep(c(last.pos + shift, last.pos + box.sides), each=2),
#            c(mid.height + lwidth, mid.height, lower.height + box.sides/2 - lwidth/2, lower.height + box.sides/2 + lwidth/2),
#            col="black", border=NA) 
#    textAtMid(last.pos + box.sides, lower.height + box.sides/2 - lwidth/2, last.pos + shift, mid.height, nDE, cex=1.4, pos=1, adj=-0.8)

    # Adding explanation.
    last.pos <- 0
    extras <- 1
    line.pos <- -2

    lines(rep(c(last.pos + box.sides-extras, last.pos+shift+extras), each=2),
          c(line.pos, line.pos - extras, line.pos - extras, line.pos), lwd=1.5)
    text(last.pos+(box.sides+shift)/2, line.pos-extras-1, pos=1, "Transfection +\noff-target effects", cex=1.4)

#    last.pos <- last.pos + shift
#    lines(rep(c(last.pos + box.sides-extras, last.pos+shift+extras), each=2),
#          c(line.pos, line.pos - extras, line.pos - extras, line.pos), lwd=1.5)
#    text(last.pos+(box.sides+shift)/2, line.pos-extras-1, pos=1, "Knockdown + \noff-target effects", cex=1.4)

    dev.off()
}
