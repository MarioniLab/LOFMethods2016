########################################################################################

# First batch of data:
firstlot <- read.table("../../real_20160208/analysis/sample_metadata.txt", header=TRUE, sep ="\t")
firstlot <- firstlot[,c(1,7)]

lib.num <- firstlot$Library
exp.num <- sub(".*exp", "\\1", firstlot$Individual)
exp.num[exp.num=="4b"] <- "5"
condition <- sub("[ _-]*exp.*", "", firstlot$Individual)
cleaned.first <- data.frame(Library=lib.num, Condition=condition, Experiment=exp.num) 

# Discarding unnecessary libraries.
discard <- grepl("^271", condition) | grepl("^3417", condition) |  # Unnecessary libraries
           grepl("289 LNA-1", condition) | # LNA with poor knockdown.
           firstlot$Library=="do8276" # Failed sequencing

cleaned.first <- cleaned.first[!discard,] 
rownames(cleaned.first) <- NULL
cleaned.first$Date <- 20160208
cleaned.first$Batch <- "20160208"

# Setting LOF mode:
LOF <- character(nrow(cleaned.first))
LOF[grepl("guide", cleaned.first$Condition) | 
    grepl("cas9", cleaned.first$Condition)] <- "CRISPRi"
LOF[grepl("Ambion", cleaned.first$Condition) |
    grepl("Dharamaco", cleaned.first$Condition) | 
    grepl("siRNA", cleaned.first$Condition)] <- "RNA interference"
LOF[cleaned.first$Condition=="neg control B" |
    cleaned.first$Condition=="289 LNA2"] <- "LNA"
LOF[cleaned.first$Condition=="cells"] <- "none"

# Setting genotype:
genotype <- rep("wild type genotype", nrow(cleaned.first))
genotype[grepl("guide", cleaned.first$Condition) | 
         grepl("cas9", cleaned.first$Condition)] <- "dCas9-KRAB clone 2"

# Setting treatment compound:
compound <- character(nrow(cleaned.first))

compound[cleaned.first$Condition=="Con si Ambion"] <- "Ambion control"
compound[cleaned.first$Condition=="Control Dharamaco"] <- "Dharmacon control"
compound[cleaned.first$Condition=="289 siRNA"] <- "289 siRNA"
 
compound[cleaned.first$Condition=="289 LNA2"] <- "289 LNA"
compound[cleaned.first$Condition=="neg control B"] <- "Negative control B"

compound[cleaned.first$Condition=="289 guide 1"] <- "289 guide 1"
compound[cleaned.first$Condition=="289 guide1"] <- "289 guide 1"
compound[cleaned.first$Condition=="289 guide 9"] <- "289 guide 9"
compound[cleaned.first$Condition=="289 guide9"] <- "289 guide 9"
compound[cleaned.first$Condition=="Negative guide 2"] <- "Negative control guide 2"

compound[cleaned.first$Condition %in% c("Hela cas9 clone 2", "cells")] <- "none"

# Creating a new group.
cleaned.first$LOF <- LOF
cleaned.first$Genotype <- genotype 
cleaned.first$Compound <- compound

########################################################################################

# Second batch of data:
secondlot <- read.csv("../../real_20160713/analysis/sample_metadata.csv")
secondlot <- secondlot[,1:2]
secondlot[,1] <- tolower(secondlot[,1])

lib.num <- secondlot$DO_name
exp.num <- sub(".*exp", "\\1", secondlot$sample_name)
condition <- sub("_exp.*", "", secondlot$sample_name)
cleaned.second <- data.frame(Library=lib.num, Condition=condition, Experiment=exp.num) 

# Sub batch within a sequencing batch.
first.sub.batch <- exp.num %in% as.character(6:9)
cleaned.second$Date <- 20160713
cleaned.second$Batch <- ifelse(first.sub.batch, "20160713b", "20160713") 

# Discarding unnecessary libraries:
discard <- grepl("centrosome", condition) |  # Unnecessary libraries
           grepl("(LNA|Neg_control*|transfection)", condition) | (condition=="Hela_cells" & !first.sub.batch) | # Failed because of Cas9 contamination
           grepl("H19guide1", condition) # CRISPRi with poor knockdown.
cleaned.second <- cleaned.second[!discard,]

# Setting LOF mode:
LOF <- character(nrow(cleaned.second))
LOF[grepl("Cas9", cleaned.second$Condition) | grepl("guide", cleaned.second$Condition)] <- "CRISPRi"
LOF["Hela_cells"==cleaned.second$Condition] <- "none"

# Setting genotype:
genotype <- rep("wild type genotype", nrow(cleaned.second))
genotype[grepl("guide", cleaned.second$Condition) | 
         cleaned.second$Condition=="HelaCas9_clone2"] <- "dCas9-KRAB clone 2"
genotype[cleaned.second$Condition == "HelaCas9_clone1"] <- "dCas9-KRAB clone 1"
genotype[cleaned.second$Condition == "HelaCas9_clone4"] <- "dCas9-KRAB clone 4"

# Setting compound:
compound <- character(nrow(cleaned.second))
compound[cleaned.second$Condition=="Negative_guide1"] <- "Negative control guide 1"
compound[cleaned.second$Condition=="Negative_guide2"] <- "Negative control guide 2"
compound[cleaned.second$Condition=="H19guide2"] <- "H19 guide 2"
compound[grepl("clone", cleaned.second$Condition) | cleaned.second$Condition=="Hela_cells"] <- "none"

# Creating a new group.
cleaned.second$LOF <- LOF
cleaned.second$Genotype <- genotype 
cleaned.second$Compound <- compound

########################################################################################

# Third batch of data:
thirdlot <- read.csv("../../real_20160907/analysis/sample_metadata.csv")
thirdlot[,1] <- tolower(thirdlot[,1])

lib.num <- thirdlot$Library
exp.num <- sub("_LNAold", "", sub(".*exp", "\\1", thirdlot$Sample))
condition <- sub("_exp.*", "", thirdlot$Sample)
cleaned.third <- data.frame(Library=lib.num, Condition=condition, Experiment=exp.num) 

# Discarding unnecessary libraries:
discard <- thirdlot$Library=="do9614" | # Failed due to sequencing
           !grepl("LNAold", thirdlot$Sample) # only looking at the samples involved in LNA.
cleaned.third <- cleaned.third[!discard,]
cleaned.third$Date <- 20160907
cleaned.third$Batch <- "20160907"

# Setting LOF mode:
LOF <- rep("LNA", nrow(cleaned.third))
LOF[cleaned.third$Condition=="Cells"] <- "none"

# Setting genotype:
genotype <- rep("wild type genotype", nrow(cleaned.third))

# Setting compound:
compound <- character(nrow(cleaned.third))
compound[cleaned.third$Condition=="Cells_Max"] <- "Transfection control"
compound[cleaned.third$Condition=="NegA"] <- "Negative control A"
compound[cleaned.third$Condition=="NegB"] <- "Negative control B"
compound[cleaned.third$Condition=="Cells"] <- "none"

# Creating a new group.
cleaned.third$LOF <- LOF
cleaned.third$Genotype <- genotype 
cleaned.third$Compound <- compound

########################################################################################

# Fourth batch of data:
fourthlot <- read.csv("../../real_20161212/analysis/metadata.csv", header=TRUE, stringsAsFactors=FALSE)
fourthlot$DO.numbers <- tolower(fourthlot$DO.numbers)

lib.num <- fourthlot$DO.numbers
exp.num <- sub("_hetero", "", sub(".*exp", "\\1", fourthlot$Description))
exp.num[exp.num=="4b"] <- "5"
condition <- sub("_exp[^_]+", "", fourthlot$Description)
cleaned.fourth <- data.frame(Library=lib.num, Condition=condition, Experiment=exp.num) 

# Discarding unnecessary libraries:
discard <- !grepl("hetero", condition) |
           grepl("271", condition) | # don't care about this one.
           grepl("289_g1", condition) # guide 1 failed.
cleaned.fourth <- cleaned.fourth[!discard,]
cleaned.fourth$Date <- 20161212
cleaned.fourth$Batch <- "20161212"

# Setting LOF mode:
LOF <- rep("CRISPRi", nrow(cleaned.fourth))
LOF[cleaned.fourth$Condition=="Hela_hetero"] <- "none"

# Setting genotype:
genotype <- rep("heterogenous dCas9-KRAB", nrow(cleaned.fourth))
genotype[cleaned.fourth$Condition=="Hela_hetero"] <- "wild type genotype"

# Setting compound:
compound <- character(nrow(cleaned.fourth))
compound[cleaned.fourth$Condition=="Hela_Nc1_hetero"] <- "Negative control guide 1"
compound[cleaned.fourth$Condition=="Hela_Nc2_hetero"] <- "Negative control guide 2"
compound[cleaned.fourth$Condition=="Hela_H19_hetero"] <- "H19 guide 2"
compound[cleaned.fourth$Condition=="Hela_289_g1_hetero"] <- "289 guide 1"
compound[cleaned.fourth$Condition=="Hela_289_g9_hetero"] <- "289 guide 9"
compound[cleaned.fourth$Condition=="Hela_BFP_hetero" | cleaned.fourth$Condition=="Hela_hetero"] <- "none"

# Creating a new group.
cleaned.fourth$LOF <- LOF
cleaned.fourth$Genotype <- genotype 
cleaned.fourth$Compound <- compound

########################################################################################

# Merging lots.
combined <- rbind(cleaned.first, cleaned.second, cleaned.third, cleaned.fourth)
combined$Batch <- as.integer(factor(combined$Batch))
new.group <- paste0(combined$LOF, ".", combined$Genotype, ".", combined$Compound, ".batch_", combined$Batch)
new.group <- gsub("[- ]", "_", new.group)
combined$Condition <- new.group
write.table(file="metadata.tsv", combined, row.names=FALSE, sep="\t", quote=FALSE)

