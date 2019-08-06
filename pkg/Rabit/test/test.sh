# /bin/sh

CMP=$PWD/test/compare.py
DATA_PATH=$PWD/data/
CMD=$PWD/src/Rabit 

cd ${DATA_PATH}

$CMD -x RBP_motif.hg19_3UTR.targets.mat -y TCGA.all_cancers.Expression -b hg19.background.PWM.RBP -c CNA:TCGA.all_cancers.CNA -c Methylation:TCGA.all_cancers.DNA_Methylation -o test.output

if [ $? -ne 0 ]         # Test exit status of "Rabit"
then
  echo "Rabit failure."
  exit 99
fi

python $CMP output test.output
STATUS=$?

rm test.output*

exit $STATUS
