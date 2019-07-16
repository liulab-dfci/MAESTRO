import sys
import os.path
import gzip

if len (sys.argv) != 3:
    print "usage:\t" + sys.argv[0] + " <sql file> <output dir>"
    exit(1)

in_file = sys.argv[1]
out_dir = sys.argv[2]

chrom = -1
chromEnd = -1
chromStart = -1
tName = -1
tStart = -1
tEnd = -1
txStart = -1
txEnd = -1
genoName = -1
genoStart = -1
genoEnd = -1

base_name = os.path.basename(in_file)[:-4]
data_name = in_file[:-4] + '.txt.gz'
out_file_name = out_dir + '/' + base_name + '.bed.gz'


look = 0
col = 0
for l in open(in_file):
    if "KEY" in l:
        break
    elif look == 1:
        col += 1
        if "`chrom`" in l:
            chrom = col
        elif "`tName`" in l:
            tName = col
        elif "`chromEnd`" in l:
            chromEnd = col
        elif "`chromStart`" in l:
            chromStart = col
        elif "`tStart`" in l:
            tStart = col
        elif "`tEnd`" in l:
            tEnd = col
        elif "`txStart`" in l:
            txStart = col
        elif "`txEnd`" in l:
            txEnd = col
        elif "`genoName`" in l:
            genoName = col
        elif "`genoStart`" in l:
            genoStart = col
        elif "`genoEnd`" in l:
            genoEnd = col
    elif 'CREATE TABLE' in l:
        look = 1

cols=[]

if chrom != -1:
    cols.append(chrom)
elif  tName != -1:
    cols.append(tName)
elif  genoName != -1:
    cols.append(genoName)
else:
    cols.append(-1)

if (chromStart != -1) and (chromEnd != -1):
    cols.append(chromStart)
    cols.append(chromEnd)
elif (txStart != -1) and (txEnd != -1):
    cols.append(txStart)
    cols.append(txEnd)
elif (tStart != -1) and (tEnd != -1):
    cols.append(tStart)
    cols.append(tEnd)
elif (genoStart != -1) and (genoEnd != -1):
    cols.append(genoStart)
    cols.append(genoEnd)
else:
    cols.append(-1)
    cols.append(-1)


if -1 not in cols:
    out_file = gzip.open(out_file_name, 'wb')
    for l in gzip.open(data_name, 'rb'):
        A = l.rstrip().split('\t')
        out = []
        for col in cols:
            out.append(A[col - 1])
        out_file.write('\t'.join(out) + '\n')
    out_file.close()
