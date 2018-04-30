# Gets the UTRs from the DEGs in the negative control vs cells comparisons, 
# split by downregulation (possibly direct) and upregulation (indirect).

path <- "../differential/results_lfc"
threshold <- readRDS("../differential/threshold_lfc_control.rds")
dir.create("motifs")

# Function to obtain UTR sequences.

library(BSgenome.Hsapiens.UCSC.hg38)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(org.Hs.eg.db)
    
threeutr <- threeUTRsByTranscript(TxDb.Hsapiens.UCSC.hg38.knownGene)
FUN <- function(ensembls) {
    entrez <- mapIds(org.Hs.eg.db, keys=ensembls, keytype="ENSEMBL", column="ENTREZID")
    tids <- mapIds(TxDb.Hsapiens.UCSC.hg38.knownGene, keys=entrez, keytype="GENEID", column="TXID")
    
    keep <- tids %in% names(threeutr)
    threeutr <- threeutr[as.character(tids[keep])]
    
    out <- getSeq(BSgenome.Hsapiens.UCSC.hg38, names=unlist(threeutr))
    names(out) <- make.unique(names(out))
    out[width(out) > 10]
}

# Ambion.

RNAi <- read.table(file.path(path, "siRNA_Ambion_vs_trans.txt"), header=TRUE, stringsAsFactors=FALSE)
up <- RNAi[RNAi$logFC > 0 & RNAi$P.Value <= threshold,]
down <- RNAi[RNAi$logFC < 0 & RNAi$P.Value <= threshold,]

writeXStringSet(FUN(rownames(up)), filepath="motifs/ambion_up.fa", format="fasta")
writeXStringSet(FUN(rownames(down)), filepath="motifs/ambion_down.fa", format="fasta")

# Control B.

LNA <- read.table(file.path(path, "LNA_controlB_vs_trans.txt"), header=TRUE, stringsAsFactors=FALSE)
up <- LNA[LNA$logFC > 0 & LNA$P.Value <= threshold,]
down <- LNA[LNA$logFC < 0 & LNA$P.Value <= threshold,]

writeXStringSet(FUN(rownames(up)), filepath="motifs/controlB_up.fa", format="fasta")
writeXStringSet(FUN(rownames(down)), filepath="motifs/controlB_down.fa", format="fasta")

# The idea is to analyze these sequences with DREME, using the upregulated UTRs as the controls and the downregulated UTRs as the search space.
# cd motifs
# CMD=~/Downloads/meme_4.12.0/scripts/dreme
# ${CMD} -o ambion_up -p ambion_up.fa -n ambion_down.fa
# ${CMD} -o ambion_down -p ambion_down.fa -n ambion_up.fa
# ${CMD} -o controlB_up -p controlB_up.fa -n controlB_down.fa
# ${CMD} -o controlB_down -p controlB_down.fa -n controlB_up.fa

