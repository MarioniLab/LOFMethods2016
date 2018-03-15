# This formats the count matrix in a manner that is comparable to the sdrf.tsv file.

blah <- read.table("../genic_counts.tsv", header=TRUE, check.names=FALSE)
genedata <- blah[,1:2]
blah <- blah[,-c(1:2)]
renamed <- sub("[^\\.]+\\.", "", colnames(blah))
sdrf <- read.table("sdrf.tsv",header=TRUE, sep="\t")
all.files <- unique(sub("_[12]\\.fq\\.gz", "", sdrf$Array.Data.File))

m <- match(sub("-", ".", all.files), renamed)
stopifnot(all(!is.na(m)))
blah <- blah[,m]
colnames(blah) <- all.files
blah <- cbind(genedata, blah[,m])
write.table(file="lncRNA_counts.tsv", blah, row.names=FALSE, sep="\t", quote=FALSE)




