library(edgeR)
yold <- readRDS("../../analysis/differential/object.rds")
ynew <- readRDS("object.rds")

shared <- intersect(rownames(yold), rownames(ynew))
yold <- yold[shared,]
ynew <- ynew[shared,]
y <- cbind(yold, ynew)

#############################
# Removing the "block" effect. 

g <- factor(y$samples$group)
b <- factor(y$samples$block)
design <- model.matrix(~0 + g + b)
design <- design[,!colnames(design) %in% c("b24", "b39", "b44", "b56", "b62", "b67")]

v.all <- voomWithQualityWeights(y, design)
fit <- lmFit(v.all, design)

coefs <- fit$coefficients
resids <- v.all$E - fitted(fit)
keep <- grepl("^g", colnames(coefs))
coefs[,!keep] <- 0

# Removing the batch effect, by centering all "cells-only" groups.
is.cell <- grepl("^gnone.wild_type_genotype.none", colnames(coefs))
for (b in seq_len(5)) {
    in.batch <- grepl(paste0("batch_", b, "$"), colnames(coefs))
    cur.cell <- in.batch & is.cell
    if (any(cur.cell)) { 
        coefs[,in.batch] <- coefs[,in.batch] - coefs[,cur.cell]
    }
}

# Centering the clone groups for batch 2.
is.clone.b1 <- "gCRISPRi.dCas9_KRAB_clone_2.none.batch_1"==colnames(coefs)
b1 <- grepl("batch_1$", colnames(coefs))
is.clone.b2 <- "gCRISPRi.dCas9_KRAB_clone_2.none.batch_2"==colnames(coefs)
b2 <- grepl("batch_2$", colnames(coefs))
coefs[,b2] <- coefs[,b2] - coefs[,is.clone.b2] + coefs[,is.clone.b1]

# Centering the Dharmacon groups for batch 6.
is.dharam.b1 <- "gRNA_interference.wild_type_genotype.Dharmacon_control.batch_1"==colnames(coefs)
b1 <- grepl("batch_1$", colnames(coefs))
is.dharam.b6 <- "gRNA_interference.wild_type_genotype.Dharmacon_control.batch_6"==colnames(coefs)
b6 <- grepl("batch_6$", colnames(coefs))
coefs[,b6] <- coefs[,b6] - coefs[,is.dharam.b6] + coefs[,is.dharam.b1]

# Centering the LNA controls.
is.lnaa.b4 <- "gLNA.wild_type_genotype.Negative_control_A.batch_4"==colnames(coefs)
b4 <- grepl("batch_4$", colnames(coefs))
is.lnaa.b6 <- "gLNA.wild_type_genotype.Negative_control_A.batch_6"==colnames(coefs)
b6 <- grepl("batch_6$", colnames(coefs))
coefs[,b6] <- coefs[,b6] - coefs[,is.lnaa.b6] + coefs[,is.lnaa.b4]

# Centering the CRISPRi negative guide 2.
is.nc2.b5 <- "gCRISPRi.heterogeneous_dCas9_KRAB.Negative_control_guide_2.batch_5"==colnames(coefs)
b5 <- grepl("batch_5$", colnames(coefs))
is.nc2.b6 <- "gCRISPRi.heterogeneous_dCas9_KRAB.Negative_control_guide_2.batch_6"==colnames(coefs)
b6 <- grepl("batch_6$", colnames(coefs))
coefs[,b6] <- coefs[,b6] - coefs[,is.nc2.b6] + coefs[,is.nc2.b5]

# Computing "corrected" expression values.
corrected <- coefs %*% t(design) + resids

#############################
# Color by tech.

# Running the PCA after an ANODEV to detect top 1000 genes.
con <- matrix(0, ncol(design), sum(keep))
diag(con) <- -1
con[1,] <- 1
con <- con[,-1]

fitx <- contrasts.fit(fit, con)
fitx <- eBayes(fitx, robust=TRUE)
top <- topTable(fitx, n=1000)

X <- prcomp(t(corrected[rownames(top),]))
var.exp <- (X$sdev)^2
pct.exp <- var.exp/sum(var.exp) * 100

# Plotting.
pdf("pics/all_pca.pdf", width=8, height=5)
tech.colors <- list(RNAi="#7F4098", LNA="#28A8E0", CRISPRi="#F79420", CRISPRi.het="#6D3F19")

# Color by technology.
color <- rep("grey", ncol(y))
color[grepl("^CRISPRi.dCas9", y$samples$group)] <- tech.colors$CRISPRi
color[grepl("^CRISPRi.het", y$samples$group) | grepl("Malat84|TOG94", y$samples$group)] <- tech.colors$CRISPRi.het
color[grepl("^LNA", y$samples$group) | grepl("MALAT_LNA", y$samples$group)] <- tech.colors$LNA
color[grepl("^RNA", y$samples$group) | grepl("Con_Dharm", y$samples$group) | grepl("TOGsi", y$samples$group)] <- tech.colors$RNAi
color[grepl("Transfection", y$samples$group)] <- "black"

layout(cbind(1,2), width=c(5,2.5))
par(mar=c(5.1, 4.1, 0.1, 1.0))
plot(X$x[,1], X$x[,2], col=color, pch=16, 
     xlab=sprintf("PC1 (%.1f%%)", pct.exp[1]), 
     ylab=sprintf("PC2 (%.1f%%)", pct.exp[2])) 

par(mar=c(5.1, 0.1, 0.1, 0.0))
plot.new()
legend("topleft", 
       col=c("grey", "black", unlist(tech.colors)),
       pch=16,
       legend=c("HeLa cells", 
                "Cells + transfection reagent", 
                "RNAi",
                "LNA",
                "CRISPRi (clonal)", 
                "CRISPRi (non-clonal)"
                )
       )

# Shape by knockdown.
target.colors <- list(`289`="black", H19="gold", TOG="cadetblue", MALAT1="indianred")

color <- rep("grey", ncol(y))
color[grepl("289", y$samples$group)] <- target.colors$`289`
color[grepl("H19", y$samples$group)] <- target.colors$H19
color[grepl("TOG", y$samples$group)] <- target.colors$TOG
color[grepl("MALAT", ignore.case=TRUE, y$samples$group)] <- target.colors$MALAT1

layout(cbind(1,2), width=c(5,2.5))
par(mar=c(5.1, 4.1, 0.1, 1.0))
plot(X$x[,1], X$x[,2], col=color, pch=16, 
     xlab=sprintf("PC1 (%.1f%%)", pct.exp[1]), 
     ylab=sprintf("PC2 (%.1f%%)", pct.exp[2])) 

par(mar=c(5.1, 0.1, 0.1, 0.0))
plot.new()
legend("topleft", 
       col=c("grey", unlist(target.colors)),
       pch=16,
       legend=c("No knockdown", "289", "H19", "Ch-TOG", "MALAT1"))

dev.off()
