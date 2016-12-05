cd fastq
for x in $(ls)
do 
    echo $x
ftp -n -v ftp-private.ebi.ac.uk << EOT
ascii
pass
user aexpress aexpress1
prompt
cd ivnpoy2h-3d6jcuf6fpuxi
put $x
bye
EOT
done

cd -
ftp -n -v ftp-private.ebi.ac.uk << EOT
ascii
pass
user aexpress aexpress1
prompt
cd ivnpoy2h-3d6jcuf6fpuxi
put lncRNA_counts.tsv
bye
EOT

