# This script links over all relevant BAM and FastQ files from the original directories.

unlink("bam", recursive=TRUE)
dir.create("bam", showWarning=FALSE)
unlink("fastq", recursive=TRUE)
dir.create("fastq", showWarning=FALSE)
to.pull <- read.table("metadata.tsv", header=TRUE, stringsAsFactor=FALSE)

for (i in seq_len(nrow(to.pull))) { 
    dir <- paste0("../real_", sub("[a-z]+$", "", to.pull$Batch[i]))

    bam.dir <- file.path(dir, "bam")
    chosen <- list.files(bam.dir, pattern=to.pull$Library[i], full=TRUE)
    stopifnot(length(chosen)>1)
    destination <- file.path("bam", paste0(to.pull$Batch[i], ".", basename(chosen)))
    file.symlink(file.path("..", chosen), destination) 

    fastq.dir <- file.path(dir, "fastq")
    chosen <- list.files(fastq.dir, pattern=to.pull$Library[i], full=TRUE)
    stopifnot(length(chosen)>1)
    destination <- file.path("fastq", paste0(to.pull$Batch[i], ".", basename(chosen)))
    file.symlink(file.path("..", chosen), destination) 
}
