cd fastq
for x in $(ls)
do 
    echo $x
ftp -n -v ftp-private.ebi.ac.uk << EOT
binary
pass
user aexpress aexpress1
prompt
cd E-MTAB-5308
put $x
bye
EOT
done

cd -
ftp -n -v ftp-private.ebi.ac.uk << EOT
binary
pass
user aexpress aexpress1
prompt
cd E-MTAB-5308
put lncRNA_counts.tsv
bye
EOT

