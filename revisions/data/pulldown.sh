# This Bash script pulls FASTQ files and MD5 sums from tier II storage.
# This is necessary as crinode108 cannot use R, thus requiring a complete
# implementation of the code in Bash. Note: this will not be applicable to 
# any user who is not the author!

set -e
set -u

# Requesting a striped FASTQ folder.
if [ ! -d fastq ]
then
    echo "No fastq folder, make it striped"
    exit
fi

# Figuring out where to place the MD5 sums.
md5dest=fastq/md5.all
if [ -e $md5dest ]
then
    rm $md5dest
fi

tmp_meta=fastq/contents.csv
scp hpcgate:/mnt/nas-srv002/jmlab/group_folders/lun01/Odom/lncRNA_mitosis/revisions_20180227/contents.csv ${tmp_meta}

tmp_md5=fastq/md5.all
scp hpcgate:/mnt/nas-srv002/jmlab/group_folders/lun01/Odom/lncRNA_mitosis/revisions_20180227/fastq/md5.all ${tmp_md5}

keep=$(cat metadata.tsv | cut -f1 | tail -n +2)

# Pulling down the files, and renaming them to fit the "do" convention.
for f in ${keep[@]}
do
    ref=\"${f}\"
    index=$(cat ${tmp_meta} | awk -F',' -v VAR=$ref 'tolower($4) == VAR { print $2 }' | sed "s/-/_/" | sed "s/\"//g")
    
    collected=$(ls fastq/ | grep ${f} | wc -l)
    if [ ${collected} -eq 0 ]
    then
        scp hpcgate:/mnt/nas-srv002/jmlab/group_folders/lun01/Odom/lncRNA_mitosis/revisions_20180227/fastq/SLX-15604.${index}* fastq
        
        for x in $(ls fastq/SLX-15604.${index}*)
        do 
            prefix=$(basename ${x})
            mv $x fastq/${f}_${prefix}
        done
    fi

    # Also grabbing their MD5 sums.
    cat ${tmp_md5} | grep ${index} | sed "s/SLX-15604/${f}_SLX-15604/" >> ${md5dest}
done
rm ${tmp_md5}
rm ${tmp_meta}
