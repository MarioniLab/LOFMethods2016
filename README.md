# LOFMethods2016
Code for Lovorka and Aaron's lncRNA loss-of-function project.

The `internal` subdirectory contains the code used to produce the analyses in the manuscript.

- The `mapme.sh` script performs alignment of Gzipped FASTQ files to the hg38 genome, using the code in the `cruktools` directory of the `OddsAndEnds` [repository](https://github.com/LTLA/OddsAndEnds).
For internal use, the `pullover.R` script creates links to FASTQ and BAM files in their original directories, which avoids having to regenerate all the files.
- The `count_me.R` script in `analysis/` is generated from `count_me.sh` (using the `OddsAndEnds` scripts again) and counts the number of reads mapped to each gene.
For internal use, the `combine_metadata.R` script pools metadata information from the original directories.
- The `run_order.R` script in `analysis/differential/` executes all the Rmarkdown scripts to perform the DE (`de_analysis.Rmd`) and pathway analyses (`kegg.Rmd`).
It also executes the various R scripts to create the necessary plots (`make_*.R`).
This must be done in the correct order.
- The `genomic_tracker.R` script in `analysis/genomic/` looks at the genomic context of the DE genes for the CRISPRi comparisons.
The Bash scripts in `Thakore_peaks` align and call peaks from the Thakore _et al._ data.
- The `make_h19.R` script in `analysis/h19/` makes a barplot of expression for known target genes of H19.
- The various R scripts in `analysis/arrayexpress/` prepare the data for submission to ArrayExpress.



