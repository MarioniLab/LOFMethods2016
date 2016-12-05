# This makes a plot of the H19 target genes, and how they're expressed across the CRISPRi controls.

all.counts <- read.table("../genic_counts.tsv", header=TRUE, row.names=1)
lengths <- all.counts[,1]
all.counts <- all.counts[,-1]
dim(all.counts)

library(edgeR)
lib.name <- sub("_.*", "", colnames(all.counts))
all.counts <- sumTechReps(all.counts, lib.name)
dim(all.counts)

metadata <- read.table("../metadata.tsv", stringsAsFactors=FALSE, header=TRUE)
m <- match(colnames(all.counts), paste0("X", metadata$Batch, ".", metadata$Library))
stopifnot(all(!is.na(m)))
metadata <- metadata[m,]

grouping <- paste0(metadata$Sample, ".", metadata$Batch)

# Identifying the genes of interest.
groups.of.interest <- c("Helacells.20160713b",
                        "HelaCas9clone1.20160713b",
                        "HelaCas9clone2.20160713b",
                        "HelaCas9clone4.20160713b",
                        "Negativeguide1.20160713",
                        "Negativeguide2.20160713",
                        "H19guide2.20160713")

genes.of.interest <- c("ENSG00000167244",
                       "ENSG00000197081",
                       "ENSG00000129757",
                       "ENSG00000185559",
                       "ENSG00000087460",
                       "ENSG00000254656",
                       "ENSG00000100697",
                       "ENSG00000185950",
                       "ENSG00000171105",
                       "ENSG00000149948",
                       "ENSG00000114315",
                       "ENSG00000138166",
                       "ENSG00000124762",
                       "ENSG00000139687",
                       "ENSG00000054267",
                       "ENSG00000167106",
                       "ENSG00000130600")

# Creating the barplot.

adjc <- cpm(all.counts)[genes.of.interest,]
library(org.Hs.eg.db)
anno <- select(org.Hs.eg.db, keys=genes.of.interest, keytype="ENSEMBL", column="SYMBOL")
rownames(adjc) <- anno$SYMBOL[match(genes.of.interest, anno$ENSEMBL)]

collected.mean <- collected.se <- list()
for (g in groups.of.interest) {
    in.g <- grouping==g
    current <- adjc[,in.g]
    collected.mean[[g]] <- rowMeans(current)
    collected.se[[g]] <- sqrt(apply(current, 1, var)/ncol(current))
}

collected.mean <- do.call(rbind, collected.mean)
collected.se <- do.call(rbind, collected.se)

pdf("h19_targets.pdf", width=10, height=5)
par(mar=c(5.1, 4.1, 1.1, 1.1))
col <- c("grey80", "darkblue", "blue", "dodgerblue", "darkorange", "gold1", "red")
upper <- collected.mean + collected.se
x <- barplot(collected.mean, beside=TRUE, col=col, ylab="CPM", las=2, cex.names=0.8, ylim=c(0, max(upper)))
segments(x, collected.mean, x, upper)
segments(x - 0.1, upper, x + 0.1)
legend("topright", fill=col, c("HeLa only", "dCas9-KRAB clone 1", "dCas9-KRAB clone 2", "dCas9-KRAB clone 4",
                               "Negative guide 1", "Negative guide 2", "H19 guide 2"), cex=0.8)
dev.off()


