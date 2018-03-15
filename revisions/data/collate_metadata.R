# Inner join.
tab1 <- read.csv("contents.csv", header=TRUE, check.names=FALSE)
tab2 <- read.csv("sample_metadata.csv", header=TRUE, check.names=FALSE)
combined <- merge(tab1, tab2, by.x="Sample name", by.y="IDs")

# Removing pAs samples.
metadata <- combined[,c("Sample name", "Barcode", "Var.5")]
colnames(metadata) <- c("ID", "Barcode", "Group")
keep <- !grepl("_pAs$", metadata$Group) 
metadata <- metadata[keep,]

# Removing a sample that failed.
metadata <- metadata[metadata$Barcode!="D705-D507",]

# Breaking down the groups.
metadata$Condition <- sub("_exp[0-9]+$", "", metadata$Group)
metadata$Experiment <- sub("^.*_exp", "", metadata$Group)
metadata$Group <- NULL

# Removing uninteresting groups.
metadata <- metadata[!metadata$Condition %in% c("Malat86", "NOCO", "Monastrol", "hela_cells", "BFP_CAS9"),]

# Replacing dashes with underscores.
metadata$Barcode <- sub("-", "_", metadata$Barcode)

##########################################
# Setting up the metadata.

Date <- "20180227"
Batch <- 6

LOF <- character(nrow(metadata))
LOF[grepl("^Malat84", metadata$Condition) | grepl("^TOG94", metadata$Condition) | grepl("^NC2_CRISPRi", metadata$Condition)] <- "CRISPRi"
LOF[grepl("^TOGsi", metadata$Condition) | grepl("^Con Dharm", metadata$Condition)] <- "RNA interference"
LOF[grepl("^Malat LNA", metadata$Condition) | grepl("^LNA A", metadata$Condition)] <- "LNA"

Genotype <- rep("wild type genotype", nrow(metadata))
Genotype[LOF=="CRISPRi"] <- "heterogeneous dCas9-KRAB"

Compound <- character(nrow(metadata))
Compound[grepl("^Malat84", metadata$Condition)] <- "MALAT1 guide 1"
Compound[grepl("^TOG94", metadata$Condition)] <- "Ch-TOG guide 1"
Compound[grepl("^NC2_CRISPRi", metadata$Condition)] <- "Negative control guide 2"
Compound[grepl("^TOGsi", metadata$Condition)] <- "Ch-TOG siRNA"
Compound[grepl("^Con Dharm", metadata$Condition)] <- "Dharmacon control"
Compound[grepl("^Malat LNA", metadata$Condition)] <- "MALAT1 LNA"
Compound[grepl("^LNA A", metadata$Condition)] <- "Negative control A"

write.table(file="metadata.tsv", 
            data.frame(Library=tolower(metadata$ID),
                       Condition=gsub("[ -]", "_", paste(LOF, Genotype, Compound, paste0("batch_", Batch), sep=".")),
                       Experiment=metadata$Experiment,
                       Date=Date, Batch=Batch, LOF=LOF, Genotype=Genotype, Compound=Compound),                                   
            sep="\t", quote=FALSE, row.names=FALSE)
