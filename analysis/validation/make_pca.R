library(edgeR)
y <- readRDS("../differential/object.rds")

# Color by tech.
color <- rep("black", ncol(y))
color[grepl("^CRISPRi.dCas9", y$samples$group)] <- "orange"
color[grepl("^CRISPRi.het", y$samples$group)] <- "brown"
color[grepl("^LNA", y$samples$group)] <- "dodgerblue"
color[grepl("^RNA", y$samples$group)] <- "purple"
color[grepl("Transfection", y$samples$group)] <- "grey"

# Shape by knockdown.
pch <- rep(1, ncol(y))
pch[grepl("289", y$samples$group)] <- 4
pch[grepl("H19", y$samples$group)] <- 6

# Running the PCA.
lcpm <- cpm(y, log=TRUE, prior=3)
X <- prcomp(t(lcpm))
var.exp <- (X$sdev)^2
pct.exp <- var.exp/sum(var.exp) * 100

# Plotting.
pdf("all_pca.pdf")
par(mfrow=c(2,2), mar=c(5.1, 4.1, 0.1, 1.0))
plot(X$x[,1], X$x[,2], col=color, pch=pch, 
     xlab=sprintf("PC1 (%.1f%%)", pct.exp[1]), 
     ylab=sprintf("PC2 (%.1f%%)", pct.exp[2])) 

plot.new()
legend("topleft", 
       col=c("black", 
             "grey",
             "purple",
             "dodgerblue", 
             "orange", 
             "brown",
             rep("black", 3)
             ) ,
       pch=c(rep(16, 6),
             1,
             4, 
             6
             ),
       legend=c("cells", 
                "transfection", 
                "RNAi",
                "LNA",
                "CRISPRi (clonal)", 
                "CRISPRi (non-clonal)",
                "No depletion",
                "289 depletion",
                "H19 depletion"
                )
       )

plot(X$x[,1], X$x[,3], col=color, pch=pch, 
     xlab=sprintf("PC1 (%.1f%%)", pct.exp[1]), 
     ylab=sprintf("PC3 (%.1f%%)", pct.exp[3])) 

plot(X$x[,2], X$x[,3], col=color, pch=pch, 
     xlab=sprintf("PC2 (%.1f%%)", pct.exp[2]), 
     ylab=sprintf("PC3 (%.1f%%)", pct.exp[3])) 

dev.off()
