echo '
anno.files <- file.path("../../genomes/annotation", c("hg38.gtf", "cas9_pHR_approx.gtf"))
bam.files <- list.files("bam", full=TRUE, pattern="bam$")
ispet <- TRUE
strandspec <- 2
' | cat - ../../tools/counter.R > count_me.R
R CMD BATCH --no-save count_me.R 

