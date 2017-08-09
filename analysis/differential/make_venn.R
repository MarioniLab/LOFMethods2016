for (mode in c("results_de", "results_lfc")) {
    if (mode=="results_de") {
        extra <- "de"
    } else {
        extra <- "lfc"
    }

    ####################################################
    # Make a Venn diagram of clones vs cells.

    clone2 <- read.table(file.path(mode, "CRISPRi_clone2_vs_cells_II.txt"), header=TRUE, stringsAsFactors=FALSE)
    clone1 <- read.table(file.path(mode, "CRISPRi_clone1_vs_cells.txt"), header=TRUE, stringsAsFactors=FALSE)
    clone4 <- read.table(file.path(mode, "CRISPRi_clone4_vs_cells.txt"), header=TRUE, stringsAsFactors=FALSE)

    common.order <- sort(clone2$ENSEMBL)
    m2 <- match(common.order, clone2$ENSEMBL)
    m1 <- match(common.order, clone1$ENSEMBL)
    m4 <- match(common.order, clone4$ENSEMBL)

    clone2 <- clone2[m2,]
    clone1 <- clone1[m1,]
    clone4 <- clone4[m4,]

    combined.p <- c(clone1$P.Value, clone2$P.Value, clone4$P.Value)
    is.sig <- p.adjust(combined.p, method="BH") <= 0.05
    all.choices <- matrix(is.sig, nrow=nrow(clone1))
    colnames(all.choices) <- c("Clone 1", "Clone 2", "Clone 4")

    library(limma)
    pdf(sprintf("pics/venn_clones_%s.pdf", extra))
    vennDiagram(all.choices)
    dev.off()

    in.all <- rowSums(all.choices)==3L
    write.table(file=file.path(mode, "combined_clones.txt"), row.names=FALSE, sep="\t", quote=FALSE,
                data.frame(clone2[,c("ENSEMBL", "SYMBOL")], 
                           clone1.logFC=clone1$logFC, clone2.logFC=clone2$logFC, clone4.logFC=clone4$logFC, 
                           adj.P.Val=ifelse(in.all, 0, 1)))
   
    if (extra=="lfc") {
        ####################################################
        # Printing out the common genes on a heatmap.
        require(edgeR)
        y <- readRDS("object.rds")
        adjc <- cpm(y, log=TRUE, prior.count=3)

        my <- match(clone2$ENSEMBL[in.all], y$genes$ENSEMBL)
        chosen.adjc <- adjc[my,]
        my.names <- y$genes$SYMBOL[my]
        my.names[is.na(my.names)] <- y$genes$ENSEMBL[my][is.na(my.names)]
        rownames(chosen.adjc) <- my.names

        is.clone.lib <- grepl("20160713b", y$samples$group)
        chosen.adjc <- chosen.adjc[,is.clone.lib]
        full.names <- sub("\\..*", "", y$samples$group[is.clone.lib])
        full.names <- sub("Hela(Cas9)?", "", full.names)
        full.names <- sub("([^ ])([0-9])", "\\1 \\2", full.names)
        colnames(chosen.adjc) <- full.names
        
        o <- order(full.names)
        chosen.adjc <- chosen.adjc[,o]
        chosen.adjc <- chosen.adjc - rowMeans(chosen.adjc[,colnames(chosen.adjc)=="cells"]) # centre each row's Hela's to be zero.

        require(gplots)
        pdf("pics/heat_clone_common.pdf")
        heatmap.2(chosen.adjc, col=bluered, symm=TRUE, dendrogram="none", Colv=FALSE, Rowv=FALSE,
                  trace="none", margins=c(5, 8), breaks=seq(-5, 5, length.out=21),  symkey=FALSE)
        dev.off()

        # Cbinding with log-fold changes for the heterogeneous vs Cells.
        is.het.lib <- grepl("20161212", y$samples$group) & (grepl("Helahetero", y$samples$group) | grepl("HelaBFPhetero", y$samples$group))
        het.chosen.adjc <- adjc[my, is.het.lib]
        het.full.names <- sub("\\..*", "", y$samples$group[is.het.lib])
        het.full.names <- sub("Hela(Cas9)?", "", het.full.names)
        het.full.names <- sub("([^ ])([0-9])", "\\1 \\2", het.full.names)
        colnames(het.chosen.adjc) <- het.full.names
        rownames(het.chosen.adjc) <- my.names
        
        o <- order(het.full.names)
        het.chosen.adjc <- het.chosen.adjc[,o]
        het.chosen.adjc <- het.chosen.adjc - rowMeans(het.chosen.adjc[,colnames(het.chosen.adjc)=="hetero"]) 

        pdf("pics/heat_het_common.pdf")
        heatmap.2(het.chosen.adjc, col=bluered, symm=TRUE, dendrogram="none", Colv=FALSE, Rowv=FALSE,
                  trace="none", margins=c(5, 8), breaks=seq(-5, 5, length.out=21), symkey=FALSE)
        dev.off()
    }

    ####################################################
    # Make a Venn diagram of all techniques against each other.

    CRISPRi.het <- read.table(file.path(mode, "combined_CRISPRi_het_289.txt"), header=TRUE, stringsAsFactors=FALSE)
    LNA <- read.table(file.path(mode, "combined_LNA_289.txt"), header=TRUE, stringsAsFactors=FALSE)
    CRISPRi <- read.table(file.path(mode, "combined_CRISPRi_289.txt"), header=TRUE, stringsAsFactors=FALSE)

    common.order <- sort(LNA$ENSEMBL)
    ms <- match(common.order, CRISPRi.het$ENSEMBL)
    ml <- match(common.order, LNA$ENSEMBL)
    mc <- match(common.order, CRISPRi$ENSEMBL)

    CRISPRi.het <- CRISPRi.het[ms,]
    LNA <- LNA[ml,]
    CRISPRi <- CRISPRi[mc,]

    combined.p <- c(CRISPRi.het$P.Value, LNA$P.Value, CRISPRi$P.Value)
    is.sig <- p.adjust(combined.p, method="BH") <= 0.05
    all.choices <- matrix(is.sig, nrow=nrow(CRISPRi.het))
    colnames(all.choices) <- c("CRISPRi.het", "LNA", "CRISPRi")

    library(limma)
    pdf(sprintf("pics/venn_289_%s.pdf", extra))
    vennDiagram(all.choices)
    dev.off()

    ####################################################
    # Make a Venn diagram of the siRNA comparisons.

    threshold <- readRDS(paste0("threshold_", extra, ".rds"))
    Ambionv289 <- read.table(file.path(mode, "siRNA_289_vs_Ambion.txt"), header=TRUE, stringsAsFactors=FALSE)
    AmbionvCells <- read.table(file.path(mode, "siRNA_Ambion_vs_cells.txt"), header=TRUE, stringsAsFactors=FALSE)
    Dharmaconv289 <- read.table(file.path(mode, "siRNA_289_vs_Dharmacon.txt"), header=TRUE, stringsAsFactors=FALSE)
    DharmaconvCells <- read.table(file.path(mode, "siRNA_Dharmacon_vs_cells.txt"), header=TRUE, stringsAsFactors=FALSE)
    AmbionvDharmacon <- read.table(file.path(mode, "siRNA_Ambion_vs_Dharmacon.txt"), header=TRUE, stringsAsFactors=FALSE)

    common.order <- sort(Ambionv289$ENSEMBL)
    Ambionv289 <- Ambionv289[match(common.order, Ambionv289$ENSEMBL),]
    AmbionvCells <- AmbionvCells[match(common.order, AmbionvCells$ENSEMBL),]
    Dharmaconv289 <- Dharmaconv289[match(common.order, Dharmaconv289$ENSEMBL),]
    DharmaconvCells <- DharmaconvCells[match(common.order, DharmaconvCells$ENSEMBL),]
    AmbionvDharmacon <- AmbionvDharmacon[match(common.order, AmbionvDharmacon$ENSEMBL),]

    is.sig <- cbind(AmbionvCells=AmbionvCells$P.Value <= threshold, 
                    DharmaconvCells=DharmaconvCells$P.Value <= threshold) 
    pdf(sprintf("pics/venn_siRNA_Ambion_Dharmacon_vs_Cells_%s.pdf", extra))
    vennDiagram(is.sig)
    dev.off()

    is.sig <- cbind(
                    Ambionv289=Ambionv289$P.Value <= threshold, 
#                    Dharmaconv289=Dharmaconv289$P.Value <= threshold,
                    AmbionvDharmacon=AmbionvDharmacon$P.Value <= threshold
                    ) 
    pdf(sprintf("pics/venn_siRNA_289_Dharmacon_vs_Ambion_%s.pdf", extra))
    vennDiagram(is.sig)
    dev.off()

    # Doing the same for the LNAs.
    LNAvA <- read.table(file.path(mode, "LNA_289.2_vs_controlA.txt"), header=TRUE, stringsAsFactors=FALSE)
    LNAvB <- read.table(file.path(mode, "LNA_289.2_vs_controlB.txt"), header=TRUE, stringsAsFactors=FALSE)
    AvB <- read.table(file.path(mode, "LNA_controlA_vs_controlB.txt"), header=TRUE, stringsAsFactors=FALSE)

    common.order <- sort(LNAvA$ENSEMBL)
    LNAvA <- LNAvA[match(common.order, LNAvA$ENSEMBL),]
    LNAvB <- LNAvB[match(common.order, LNAvB$ENSEMBL),]
    AvB <- AvB[match(common.order, AvB$ENSEMBL),]

    is.sig <- cbind(
                    LNAvA=LNAvA$P.Value <= threshold, 
                    AvB=AvB$P.Value <= threshold
                    ) 
    pdf(sprintf("pics/venn_LNA_289_A_vs_B_%s.pdf", extra))
    vennDiagram(is.sig)
    dev.off()

    is.sig <- cbind(
                    LNAvB=LNAvB$P.Value <= threshold,
                    AvB=AvB$P.Value <= threshold
                    ) 
    pdf(sprintf("pics/venn_LNA_289_B_vs_A_%s.pdf", extra))
    vennDiagram(is.sig)
    dev.off()
}
