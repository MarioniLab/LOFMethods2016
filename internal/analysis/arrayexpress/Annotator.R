########################################################################################
# Print out important bits from the metadata, regarding genotype and treatment.

stuff <- read.table("../metadata.tsv", header=TRUE, stringsAsFactors=FALSE)

LOF.method <- character(nrow(stuff))
LOF.method["ConsiAmbion"==stuff$Sample | "ControlDharamaco"==stuff$Sample | "289siRNA"==stuff$Sample] <- "RNAi"
LOF.method[grepl("289LNA", stuff$Sample) | "negcontrolB"==stuff$Sample | grepl("Neg[AB]", stuff$Sample) | "CellsMax"==stuff$Sample] <- "LNA"
LOF.method[grepl("289guide", stuff$Sample) | grepl("Negativeguide2", stuff$Sample) | grepl("Hela[Cc]as9clone", stuff$Sample) |
           grepl("Negativeguide1", stuff$Sample) | grepl("H19guide2", stuff$Sample)] <- "CRISPRi"
LOF.method["cells"==stuff$Sample | "Helacells"==stuff$Sample | "Cells"==stuff$Sample] <- "None"

genotype <- rep("WT", nrow(stuff))
genotype[grepl("289guide", stuff$Sample) | grepl("Negativeguide2", stuff$Sample) | grepl("Hela[Cc]as9clone2", stuff$Sample) |
         grepl("Negativeguide1", stuff$Sample) | grepl("H19guide2", stuff$Sample)] <- "dCas9-KRAB clone 2"
genotype[stuff$Sample == "HelaCas9clone1"] <- "dCas9-KRAB clone 1"
genotype[stuff$Sample == "HelaCas9clone4"] <- "dCas9-KRAB clone 3"

compound <- character(nrow(stuff))
compound[stuff$Sample=="ConsiAmbion"] <- "Ambion control"
compound[stuff$Sample=="ControlDharamaco"] <- "Dharmacon control"
compound[stuff$Sample=="289siRNA"] <- "289 siRNA"
compound[stuff$Sample=="289LNA2"] <- "289 LNA"
compound[stuff$Sample=="negcontrolB"] <- compound[stuff$Sample=="NegB"] <- "Negative control B"
compound[stuff$Sample=="NegA"] <- "Negative control A"
compound[stuff$Sample=="CellsMax"] <- "Transfection control"
compound[stuff$Sample=="289guide1"] <- "289 guide 1"
compound[stuff$Sample=="289guide9"] <- "289 guide 9"
compound[stuff$Sample=="Negativeguide2"] <- "Negative control guide 2"
compound[stuff$Sample=="Negativeguide1"] <- "Negative control guide 1"
compound[stuff$Sample=="H19guide2"] <- "H19 guide 2"
compound[stuff$Sample %in% c("HelaCas9clone1", "HelaCas9clone4", "HelaCas9clone2", "Helacas9clone2", "Cells", "cells", "Helacells")] <- "None"

seq.date <- sub("b$", "", stuff$Batch)
batch.id <- as.integer(factor(stuff$Batch))
exp.num <- sub("^exp", "", stuff$Experiment)
exp.num[exp.num=="4b"] <- "5"

holding <- data.frame(Library=stuff$Library, Sample=stuff$Sample, LOF=LOF.method, 
                      Genotype=genotype, Compound=compound, Date=seq.date, 
                      Batch=batch.id, Experiment=exp.num, stringsAsFactors=FALSE)
write.csv(file="compressed_annotation.csv", holding, row.names=FALSE, quote=FALSE)

########################################################################################
# Getting all the file names.

relink <- "make_links.sh"
write(file=relink, c("set -e", "set -u", "mkdir fastq"), ncol=1)

fpath <- "/run/user/1753941046/gvfs/smb-share:server=jmlab-data,share=jmlab/group_folders/lun01/Odom/lncRNA_mitosis"
collected <- list()
all.md5 <- list()
for (x in seq_len(nrow(holding))) {
    cur.date <- holding$Date[x]
    cur.lib <- holding$Library[x]
    curpath <- file.path(fpath, paste0("real_", cur.date), "fastq")
    chosen <- list.files(curpath, pattern=cur.lib)

    new.names <- sub("p([12]\\.fq\\.gz)$", "_\\1", chosen)
    current <- data.frame(FileName=new.names) # Replicate=as.integer(sub(".*_CRI([0-9]+)p[12]\\.fq\\.gz$", "\\1", chosen)))
    write(file=relink, paste0("ln -s ", file.path(curpath, chosen), " ", file.path("fastq", new.names)), append=TRUE, ncol=1)

    curmd5 <- all.md5[[cur.date]]
    if (is.null(curmd5)) {
        all.md5[[cur.date]] <- read.table(file.path(curpath, "md5.all"), stringsAsFactors=FALSE)
        curmd5 <- all.md5[[cur.date]]
    }
    current$MD5sum <- curmd5[,1][match(chosen, curmd5[,2])]
    collected[[x]]  <- suppressWarnings(cbind(current, holding[x,]))
}                                                                

final <- do.call(rbind, collected)

########################################################################################
# Constructing the sdrf.tsv file.

output <- list()
output[["Source Name"]] <- final$Library
output[["Characteristics[organism]"]] <- "Homo sapiens"
output[["Characteristics[cell line]"]] <- "HeLa"
output[["Material Type"]] <- "RNA"
output[[paste0(rep(c("Protocol REF", "Performer"), 6), collapse="\t")]] <- paste0(c("Obtaining HeLa cells", "Lovorka Stojic",
                                                                                    "Culturing HeLa cells", "Lovorka Stojic",
                                                                                    "Depleting the lncRNA", "Lovorka Stojic",
                                                                                    "Reverse transcription", "Lovorka Stojic",
                                                                                    "Extracting RNA", "Lovorka Stojic",
                                                                                    "Creating libraries","Lovorka Stojic"
                                                                                    ), collapse="\t")
output[["Extract Name"]] <- final$Library
output[["Comment[LIBRARY_LAYOUT]"]] <- "PAIRED"
output[["Comment[LIBRARY_SELECTION]"]] <- "Inverse rRNA"
output[["Comment[LIBRARY_SOURCE]"]] <- "TRANSCRIPTOMIC"
output[["Comment[LIBRARY_STRAND]"]] <- "first strand"
output[["Comment[LIBRARY_STRATEGY]"]] <- "ssRNA-seq"
output[["Comment[NOMINAL_LENGTH]"]] <- 295
output[["Comment[NOMINAL_SDEV]"]] <- 25
output[["Comment[ORIENTATION]"]] <- "5'-3'-3'-5'"
output[["Protocol REF\tPerformer"]] <- "Sequencing libraries\tLovorka Stojic"
output[["Assay Name"]] <- final$Library
output[["Technology Type"]] <- "sequencing assay"
output[["Comment[batch number]"]] <- final$Batch
output[["Comment[sequencing date]"]] <- paste0(substr(final$Date, 1, 4), "-", substr(final$Date, 5, 6), "-", substr(final$Date, 7, 8))
output[["Comment[experiment number]"]] <- final$Experiment
output[["Array Data File"]] <- final$FileName
output[["Protocol REF"]] <- "Assigning reads to genes"
output[["Derived Array Data File"]] <- "lncRNA_counts.tsv"
output[["Comment[MD5]"]] <- final$MD5sum
output[["Factor Value[loss of function method]"]] <- final$LOF
output[["Factor Value[genotype]"]] <- final$Genotype
output[["Factor Value[compound]"]] <- final$Compound

output$check.names <- FALSE
sdrf <- do.call(data.frame, output)
write.table(file="sdrf.tsv", sdrf, row.names=FALSE, sep="\t", quote=FALSE)
