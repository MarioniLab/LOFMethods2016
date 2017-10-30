if [ ! -d fastq ]
then
    echo "No fastq folder, make it striped"
    exit
fi

allbatches=$(cat metadata.tsv | cut -f4 | tail -n +2 | sed "s/[^0-9]//" | sort | uniq)

for b in ${allbatches[@]}
do
    keep=$(cat metadata.tsv | awk -v VAR=$b '$4 == VAR { print $1 }')

    for f in ${keep[@]}
    do
        scp hpcgate:/mnt/nas-srv002/jmlab/group_folders/lun01/Odom/lncRNA_mitosis/real_${b}/fastq/${f}* fastq
    done
done

