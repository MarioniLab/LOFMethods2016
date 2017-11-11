# LOFMethods2016

Code for Lovorka and Aaron's lncRNA loss-of-function project.
To reproduce the results:

- Enter the `analysis/arrayexpress/` directory and download the SDRF file from ArrayExpress (E-MTAB-5308).
Also download the processed data and unzip it in the same directory.
- The `run_order.R` script in `analysis/differential/` executes all the Rmarkdown scripts to perform the DE (`de_analysis.Rmd`) and pathway analyses (`kegg.Rmd`).
It also executes the various R scripts to create the necessary plots (`make_*.R`).
This must be done in the correct order.
- The `genomic_tracker.R` script in `analysis/genomic/` looks at the genomic context of the DE genes for the CRISPRi comparisons.
The Bash scripts in `Thakore_peaks` align and call peaks from the Thakore _et al._ data.
- The `make_h19.R` script in `analysis/h19/` makes a barplot of expression for known target genes of H19.

If you want to re-align the sequence data and regenerate the count matrix:

- Enter the `data/` subdirectory and download all FASTQ files from ArrayExpress (E-MTAB-5308) into `data/fastq`.
Running  `mapme.sh` aligns Gzipped FASTQ files to the hg38 genome, using the code in the [`CRUKtools` repository](https://github.com/LTLA/CRUKtools).
- Running `count_me.sh` generates the `count_me.R` script (using the `CRUKtools` scripts again) and counts the number of reads mapped to each gene.
This was done using Ensembl GRCm38 annotation version 90, see https://github.com/MarioniLab/CommonResources for details on how the GTF files were generated.
- Move the resulting `data/genic_counts.tsv` file into `analysis/arrayexpress/`.
Make sure that the SDRF file is also available in this subdirectory before running `de_analysis.Rmd`.

For internal use:

- In `data/`, the `combine_metadata.R` script collates libraries that were used in this project across various batches.
Running `pullover.sh` then pulls the relevant FASTQ files from tier II storage.
- The various R scripts in `analysis/arrayexpress/` prepare the data for submission to ArrayExpress.
