# Requires the mapping file in https://github.com/LTLA/CRUKtools.
# Requires the genome builds in https://github.com/MarioniLab/CommonResources.
ispet=1
fastq=($(ls fastq/*.fq.gz))
genome=/lustre/jmlab/resources/genomes/subread/hg38_cas9

source ${HOME}/Code/mapping/multi_align.sh
