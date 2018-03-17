# Sequence sources

- pHR-SFFV-dCas9-BFP-KRAB.fa was obtained from https://www.addgene.org/46911/
- hsa.hg38.fa was obtained from UCSC (via `/lustre/reference_data/mib-cri/reference_genomes/homo_sapiens/hg38/fasta/hsa.hg38.fa`)

# Index builds

This combines the hg38 genome with a plasmid sequence expressing dCas9-KRAB.

```sh
subread-buildindex -o hg38_cas9 pHR-SFFV-dCas9-BFP-KRAB.fa \
    /lustre/reference_data/mib-cri/reference_genomes/homo_sapiens/hg38/fasta/hsa.hg38.fa
```
