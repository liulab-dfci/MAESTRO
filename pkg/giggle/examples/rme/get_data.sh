mkdir data
cd data
wget http://egg2.wustl.edu/roadmap/data/byFileType/chromhmmSegmentations/ChmmModels/coreMarks/jointModel/final/all.mnemonics.bedFiles.tgz
tar zxvf all.mnemonics.bedFiles.tgz
cd ..
mkdir split
python rename.py states.txt EDACC_NAME.txt "data/*gz" "split/"
