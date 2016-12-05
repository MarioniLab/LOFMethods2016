firstlot <- read.table("../../real_20160208/analysis/sample_metadata.txt", header=TRUE, sep ="\t")
secondlot <- read.csv("../../real_20160713/analysis/sample_metadata.csv")
thirdlot <- read.csv("../../real_20160907/analysis/sample_metadata.csv")

firstlot <- firstlot[,c(1,7)]
secondlot <- secondlot[,1:2]
secondlot[,1] <- tolower(secondlot[,1])
thirdlot[,1] <- tolower(thirdlot[,1])

# Adding batch dates.
firstlot <- cbind("20160208", firstlot)
secondlot <- cbind("20160713", secondlot)
thirdlot <- cbind("20160907", thirdlot)

# Getting rid of 'LNAold' specifier at the end of the third lot's sample names.
thirdlot <- thirdlot[grep("LNAold", thirdlot[,3]),]
thirdlot[,3] <- sub("_LNAold", "", thirdlot[,3])

# Merging lots.
colnames(firstlot) <- colnames(secondlot) <- colnames(thirdlot) <- c("Batch", "Library", "Sample")
combined <- rbind(firstlot, secondlot, thirdlot)

# For the second lot, clones prepared in a separate batch from the H19 samples.
combined$Batch <- as.character(combined$Batch)
is.batch2 <- combined$Batch=="20160713" & grepl("exp[6-9]", combined$Sample) 
combined$Batch[is.batch2] <- paste0(combined$Batch[is.batch2], "b")

# Cleaning up names and groupings.
exp.num <- sub(".*(exp[0-9]b?)$", "\\1", combined$Sample)
groupings <- sub("exp[0-9]b?$", "", combined$Sample)
groupings <- gsub("[ _-]", "", groupings)

# Discarding libraries that we don't want.
discard <- grepl("271", groupings) | grepl("3417", groupings) | grepl("centrosome", groupings) |  # Unnecessary libraries
           grepl("289LNA1", groupings) | # LNA with poor knockdown.
           (combined$Library=="do8276" & combined$Batch=="20160208") | (combined$Library=="do9614" & combined$Batch=="20160907") | # Failed due to sequencing
           (grepl("(LNA|Negcontrol*|transfection|Helacells)", groupings) & combined$Batch=="20160713") | # Failed because of Cas9 contamination
###           (combined$Library %in% c("do9295", "do9300") & combined$Batch =="20160713") | # Inconsistent with other control clones, for some reason.
           grepl("H19guide1", groupings) # CRISPRi with poor knockdown.

combined$Sample <- groupings
combined <- cbind(combined, Experiment=exp.num)
write.table(file="metadata.tsv", combined[!discard,], row.names=FALSE, sep="\t", quote=FALSE)

########################################################################################
# Print out important bits from the metadata, regarding genotype and treatment.

stuff <- read.table("metadata.tsv", header=TRUE, stringsAsFactors=FALSE)

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
genotype[stuff$Sample == "HelaCas9clone4"] <- "dCas9-KRAB clone 4"

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
replicate.num <- sub("^exp", "", stuff$Experiment)
replicate.num[replicate.num=="4b"] <- "5"

write.csv(file="annotations.csv", 
          data.frame(Library=stuff$Library, Sample=stuff$Sample, LOF=LOF.method, 
                     Genotype=genotype, Compound=compound, Date=seq.date, 
                     Batch=batch.id, Replicate=replicate.num), row.names=FALSE, quote=FALSE)


