# Evaluating the specificity of LOF methods

## Overview

This repository provides the code for the paper **Specificity of RNAi, LNA and CRISPRi as loss-of-function methods in transcriptional analysis**, by Stojic et al., currently in biorXiv (https://doi.org/10.1101/234930).

## Generating the results

This repository contains code for Lovorka and Aaron's lncRNA loss-of-function project.
To reproduce the results in Figures 1-4, and most of the supplementary figures:

1. Run `download.sh`, which will download the SDRF file from ArrayExpress (E-MTAB-5308).
It will also download the processed data and place it in the relevant subdirectories.
2. The `run_order.R` script in `analysis/differential/` executes all the Rmarkdown scripts to perform the DE (`de_analysis.Rmd`) and pathway analyses (`kegg.Rmd`).
It also executes the various R scripts to create the necessary plots (`make_*.R`).
This must be done in the correct order.
3. The `fdr_check.Rmd` script in `analysis/validation` performs a series of DE analyses between replicates of the same group, to check that the error rate is properly controlled.
This should be followed by `run_checker.R`, which will summarize the false positives across all comparisons.

To reproduce the results from depletion of _MALAT1_ or _Ch-TOG_ (Figures 5-6):

1. Enter the `revisions/analysis` directory and run `de_anaysis.Rmd` to perform DE analyses.
This assumes that `download.sh` has already been run.
2. Run `make_pics.R` to create the figures for the depletion results.
This includes MA plots, volcano plots and the log-fold change plots.
3. The `make_pca.R` script in `analysis/validation` will generate a PCA plot of all samples, colored and shaped according to the group in which each sample belongs.

To integrate our results with other data:

- The `genomic_tracker.R` script in `analysis/genomic/` looks at the genomic context of the DE genes for the CRISPRi comparisons.
The Bash scripts in `Thakore_peaks` align and call peaks from the Thakore _et al._ data.
- The `make_h19.R` script in `analysis/h19/` makes a barplot of expression for known target genes of H19.

## Realigning the sequence data 

If you want to re-align the sequence data and regenerate the count matrix:

1. Enter `genomes/builds` and follow the instructions in the README to create a `subread` index for hg38 with the Cas9 plasmid.
Also enter `genomes/annotation` and follow the instructions in the README to create a GTF file using Ensembl GRCm38 annotation version 90.
2. Enter the `data/` subdirectory and download all FASTQ files from ArrayExpress (E-MTAB-5308) into `data/fastq/`.
Move all FASTQ files corresponding to libraries in batch 6 to `revisions/data/fastq/`.
3. Running  `data/mapme.sh` aligns Gzipped FASTQ files using the `subread` aligner.
Running `data/count_me.sh` generates the `count_me.R` script and counts the number of reads mapped to each gene.
4. Enter `revisions/data/` and run `mapme.sh` followed by `count_me.sh`, to obtain a separate count matrix for the final batch.

## For internal use:

- In `data/`, the `combine_metadata.R` script collates metadata across various batches and standardizes their labels.
Running `pulldown.sh` then pulls the relevant FASTQ files from tier II storage.
- In `revisions/data/`, the `collate_metadata.R` script collates metadata for the libraries generated during the revision.
Running `pulldown.sh` then pulls the relevant FASTQ files from tier II storage.
- The various R scripts in `submitAE/` prepare the data for submission to ArrayExpress.
