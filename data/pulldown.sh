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

# Iterating over the batches and figuring out which files to keep.
allbatches=$(cat metadata.tsv | cut -f4 | tail -n +2 | sed "s/[^0-9]//" | sort | uniq)
for b in ${allbatches}
do
    keep=$(cat metadata.tsv | awk -v VAR=$b '$4 ~ VAR { print $1 }')

    if [ $b != "20161212" ]
    then
        # Pulling down the files.
        for f in ${keep[@]}
        do
            collected=$(ls fastq/ | grep ${f} | wc -l)
            if [ ${collected} -eq 0 ]
            then
                scp hpcgate:/mnt/nas-srv002/jmlab/group_folders/lun01/Odom/lncRNA_mitosis/real_${b}/fastq/${f}* fastq
            fi
        done
    
        # Pulling down the MD5 sums.
        tmp_md5=fastq/${b}.md5
        scp hpcgate:/mnt/nas-srv002/jmlab/group_folders/lun01/Odom/lncRNA_mitosis/real_${b}/fastq/md5.all ${tmp_md5}
        for f in ${keep[@]}
        do
            cat ${tmp_md5} | grep "${f}" >> ${md5dest}
        done
        rm ${tmp_md5}

    else 
        # Special behaviour for the final batch, which doesn't follow consistent naming, dammit.
        tmp_meta=fastq/metadata.csv
        scp hpcgate:/mnt/nas-srv002/jmlab/group_folders/lun01/Odom/lncRNA_mitosis/real_${b}/metadata.csv ${tmp_meta}

        tmp_md5=fastq/${b}.md5
        scp hpcgate:/mnt/nas-srv002/jmlab/group_folders/lun01/Odom/lncRNA_mitosis/real_${b}/fastq/md5.all ${tmp_md5}

        # Pulling down the files, and renaming them to fit the "do" convention.
        for f in ${keep[@]}
        do
            collected=$(ls fastq/ | grep ${f} | wc -l)
            if [ ${collected} -eq 0 ]
            then
                index=$(cat ${tmp_meta} | awk -F',' -v VAR=$f 'tolower($1) == VAR { print $3 }' | sed "s/-/_/")
                scp hpcgate:/mnt/nas-srv002/jmlab/group_folders/lun01/Odom/lncRNA_mitosis/real_${b}/fastq/SLX-12829.${index}* fastq
                
                for x in $(ls fastq/SLX-12829.${index}*)
                do 
                    prefix=$(basename ${x})
                    mv $x fastq/${f}_${prefix}
                done
            fi

            # Also grabbing their MD5 sums.
            cat ${tmp_md5} | grep ${index} | sed "s/SLX-12829/${f}_SLX-12829/" >> ${md5dest}
        done
        rm ${tmp_md5}
        rm ${tmp_meta}
    fi

done

