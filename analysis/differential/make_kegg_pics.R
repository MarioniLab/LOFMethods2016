#################################################################################
# Making pathway analysis plots, for the top 10 pathways in each comparison.

for (mode in c("results_de", "results_lfc")) {
    if (mode=="results_de") {
        extra <- "de"
    } else {
        extra <- "lfc"
    }
    respath <- file.path("pathways", mode)
    picpath <- file.path(respath, "pics")
    dir.create(picpath, showWarning=FALSE)

    for (f in list.files(respath, full=TRUE, pattern="kegg*")) { 
        input <- read.table(f, header=TRUE, sep="\t", row.names=1, fill=TRUE, quote="", stringsAsFactors=FALSE)
        renamed <- -log10(input$P.DE)
        names(renamed) <- input$Pathway

        maxx <- ifelse(grepl("lfc", f), 8, 10)        
        pdf(file.path(picpath, sub(".txt", sprintf("_%s.pdf", extra), basename(f))), width=12, height=8)
        par(mar=c(5.1, 20.1, 2.1, 2.1), cex.lab=1.4, cex.axis=1.2)

        chosen <- 20:1
        x <- renamed[chosen]
        y <- barplot(x, horiz=TRUE, las=1, xlab=expression("-"~Log[10]~"p-value"), xlim=c(0, maxx))
        text(x, y, pos=4, sprintf("%i/%i", input$DE[chosen], input$N[chosen]))

        legend("bottomright", sprintf("%i DE genes in total", input$Total.DE[1]), bty="n")
        dev.off()
    }
}


