library(edgeR)
y <- readRDS("../differential/object.rds")

#############################
# Removing the "block" effect. 
v.all <- voomWithQualityWeights(y, y$design)
fit <- lmFit(v.all, y$design)

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

# Computing "corrected" expression values.
corrected <- coefs %*% t(y$design) + resids

#############################
# Color by tech.
color <- rep("grey", ncol(y))
color[grepl("^CRISPRi.dCas9", y$samples$group)] <- "orange"
color[grepl("^CRISPRi.het", y$samples$group)] <- "brown"
color[grepl("^LNA", y$samples$group)] <- "dodgerblue"
color[grepl("^RNA", y$samples$group)] <- "purple"
color[grepl("Transfection", y$samples$group)] <- "black"

# Shape by knockdown.
pch <- rep(1, ncol(y))
pch[grepl("289", y$samples$group)] <- 4
pch[grepl("H19", y$samples$group)] <- 6
pch[grepl("Negative|Control", y$samples$group, ignore.case=TRUE) & !grepl("Transfection", y$samples$group)] <- 3

# Running the PCA after an ANODEV to detect top 1000 genes.
con <- matrix(0, ncol(y$design), sum(keep))
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
pdf("all_pca.pdf", width=8, height=5)
layout(cbind(1,2), width=c(5,2.5))

par(mar=c(5.1, 4.1, 0.1, 1.0))
plot(X$x[,1], X$x[,2], col=color, pch=pch, 
     xlab=sprintf("PC1 (%.1f%%)", pct.exp[1]), 
     ylab=sprintf("PC2 (%.1f%%)", pct.exp[2])) 

par(mar=c(5.1, 0.1, 0.1, 0.0))
plot.new()
legend("topleft", 
       col=c("grey", 
             "black",
             "purple",
             "dodgerblue", 
             "orange", 
             "brown",
             rep("black", 3)
             ) ,
       pch=c(rep(16, 6),
             3,
             4, 
             6
             ),
       legend=c("HeLa cells", 
                "Cells + transfection reagent", 
                "RNAi",
                "LNA",
                "CRISPRi (clonal)", 
                "CRISPRi (non-clonal)",
                "Control sequences",
                "SLC25A25-AS1 depletion",
                "H19 depletion"
                )
       )
dev.off()
