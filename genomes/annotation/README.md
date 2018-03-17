# Annotation sources

- cas9_pHR_approx.gtf was manually constructed by examining the Cas9-coding domain in https://www.addgene.org/46911/
- Homo_sapiens.GRCh38.90.gtf.gz was obtained from http://ftp.ensembl.org/pub/release-90/gtf/homo_sapiens/ and processed as described below.

```sh
zcat Homo_sapiens.GRCh38.90.gtf.gz | \
    sed -r "s/^([0-9MXY])/chr\1/" | \
    sed "s/^chrMT/chrM/g" | \
    awk '$3 == "exon"' > hg38.gtf
```
