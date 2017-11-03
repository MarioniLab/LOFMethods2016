cd ../../data/fastq
for x in $(ls)
do 
    echo $x
ftp -n -v ftp-private.ebi.ac.uk << EOT
ascii
pass
user aexpress aexpress1
prompt
cd E-MTAB-5308
put $x
bye
EOT
done

cd ../
ln -s genic_counts.tsv lncRNA_counts.tsv
ftp -n -v ftp-private.ebi.ac.uk << EOT
ascii
pass
user aexpress aexpress1
prompt
cd E-MTAB-5308
put lncRNA_counts.tsv
bye
EOT
rm lncRNA_counts.tsv
