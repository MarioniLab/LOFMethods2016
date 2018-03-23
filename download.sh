cd data/
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-5308/E-MTAB-5308.sdrf.txt
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-5308/E-MTAB-5308.processed.1.zip
unzip E-MTAB-5308.processed.1.zip
rm E-MTAB-5308.processed.1.zip
cd -

cd ../revisions/data
ln -s ../../data/E-MTAB-5308.sdrf.txt
cd -

# Also cloning the tools (using the last commit known to work).
git clone https://github.com/LTLA/CRUKTools tools
cd tools
git reset --hard 915706ac016f6f59e74b4ec266f06349293b997f
