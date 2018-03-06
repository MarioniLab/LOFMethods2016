# Inner join.
tab1 <- read.csv("contents.csv", header=TRUE, check.names=FALSE)
tab2 <- read.csv("sample_metadata.csv", header=TRUE, check.names=FALSE)
combined <- merge(tab1, tab2, by.x="Sample name", by.y="IDs")

# Removing pAs samples.
metadata <- combined[,c("Barcode", "Var.5")]
colnames(metadata) <- c("Library", "Group")
keep <- !grepl("_pAs$", metadata$Group) 
metadata <- metadata[keep,]

# Removing a sample that failed.
metadata <- metadata[metadata$Library!="D705-D507",]

# Breaking down the groups.
metadata$Condition <- sub("_exp[0-9]+$", "", metadata$Group)
metadata$Experiment <- sub("^.*_exp", "", metadata$Group)
metadata$Group <- NULL

# Removing uninterestin groups.
metadata <- metadata[!metadata$Condition %in% c("Malat86", "NOCO", "Monastrol", "hela_cells", "BFP_CAS9"),]

# Replacing dashes with underscores.
metadata$Library <- sub("-", "_", metadata$Library)
write.table(file="metadata.tsv", metadata, sep="\t", quote=FALSE, row.names=FALSE)

# Identifying columns to remove in the count table.
counts <- read.table("genic_counts.tsv", header=TRUE, check.names=FALSE)
re.names <- sub("SLX-[0-9]+\\.", "", colnames(counts))
re.names <- sub("\\..*", "", re.names)
counts <- counts[, re.names %in% c("GeneID", "Length", metadata$Library)]
write.table(file="subcounts.tsv", counts, sep="\t", quote=FALSE, row.names=FALSE)
