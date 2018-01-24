for (dir in c("results_de", "results_lfc")) {
    all.files <- list.files(dir) # as the comparisons here are a subset in 'differential'.
    cat(sprintf(">> %s\n", dir))

    for (mode in c(".", "../differential")) { 
        total <- 0L
        for (x in all.files) {
            blah <- read.table(file.path(mode, dir, x), header=TRUE, row.names=1, stringsAsFactors=FALSE)
            nsig <- sum(blah$adj.P.Val <= 0.05)
            if (nsig > 0) {
                cat(sprintf("%i:\t%s\n", nsig, x))
                total <- total + nsig
            }
        }
        cat(sprintf("Total: %i\n\n", total)) 
    }
    cat("\n")
}

