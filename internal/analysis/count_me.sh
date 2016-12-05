echo '
anno.files <- file.path("/lustre/jmlab/resources/annotation", c("processed/hg38.gtf", "original/cas9_pHR_approx.gtf"))
bam.files <- list.files("../bam", full=TRUE, pattern="bam$")
ispet <- TRUE
strandspec <- 2
' | cat - ~/Code/mapping/counter.R > count_me.R
R CMD BATCH --no-save count_me.R 

