#################################################################################
# Making pathway analysis plots, for the top 10 pathways in each comparison.

for (mode in c("results_de", "results_lfc")) {
    if (mode=="results_de") {
        extra <- "de"
    } else {
        extra <- "lfc"
    }

    for (f in list.files(file.path("pathways", mode), full=TRUE, pattern="kegg*")) { 
        input <- read.table(f, header=TRUE, sep="\t", row.names=1, fill=TRUE, quote="", stringsAsFactors=FALSE)
        renamed <- -log10(input$P.DE)
        names(renamed) <- input$Pathway

        maxx <- ifelse(grepl("lfc", f), 8, 10)        
        pdf(file.path("pics", sub(".txt", sprintf("_%s.pdf", extra), basename(f))), width=12, height=8)
        par(mar=c(5.1, 20.1, 2.1, 2.1), cex.lab=1.4, cex.axis=1.2)
        barplot(rev(renamed[1:20]), horiz=TRUE, las=1, xlab=expression("-"~Log[10]~"p-value"), xlim=c(0, maxx))
        abline(v=-log10(0.01), col="red", lwd=2, lty=2)
        dev.off()
    }
}


