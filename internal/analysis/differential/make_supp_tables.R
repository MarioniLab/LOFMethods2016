mode <- "results_lfc"

siRNA <- read.table(file.path(mode, "combined_siRNA_289.txt"), header=TRUE, stringsAsFactors=FALSE)
LNA <- read.table(file.path(mode, "combined_LNA_289.txt"), header=TRUE, stringsAsFactors=FALSE)
CRISPRi <- read.table(file.path(mode, "combined_CRISPRi_289.txt"), header=TRUE, stringsAsFactors=FALSE)

combined.p <- cbind(siRNA$P.Value, LNA$P.Value, CRISPRi$P.Value)
adj.p <- p.adjust(combined.p, method="BH")
dim(adj.p) <- dim(combined.p)

siRNA$adj.P.Val <- adj.p[,1]
siRNA <- siRNA[siRNA$adj.P.Val <= 0.05,]
LNA$adj.P.Val <- adj.p[,2]
LNA <- LNA[LNA$adj.P.Val <= 0.05,]
CRISPRi$adj.P.Val <- adj.p[,3]
CRISPRi <- CRISPRi[CRISPRi$adj.P.Val <= 0.05,]

write.table(file=file.path(mode, "SUPP_combined_siRNA_289.txt"), siRNA, row.names=FALSE, sep="\t", quote=FALSE) 
write.table(file=file.path(mode, "SUPP_combined_LNA_289.txt"), LNA, row.names=FALSE, sep="\t", quote=FALSE) 
write.table(file=file.path(mode, "SUPP_combined_CRISPRi_289.txt"), CRISPRi, row.names=FALSE, sep="\t", quote=FALSE) 

