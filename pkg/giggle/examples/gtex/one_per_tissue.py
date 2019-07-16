import sys

trans_to_cord = {}

for l in open("genes.bed", "r"):
	A = l.rstrip().split('\t')
	trans_to_cord[A[3]] = '\t'.join(A[:3])

sampid_to_smts = {}
smts_to_sampid = {}

for l in open("SAMPID_to_SMTS.txt", "r"):
	A = l.rstrip().split('\t')
	if len(A) == 2:
		sampid_to_smts[A[0]] = A[1]
		if A[1] not in smts_to_sampid:
			smts_to_sampid[A[1]] = []
		smts_to_sampid[A[1]].append(A[0])


cols = []
go = 0
files = {}
for l in sys.stdin:
	A = l.rstrip().split('\t')
	if A[0] == 'Name':
		cols = A[2:]
		go = 1
	elif go == 1:
		gene = A[0]		
		rpkms = A[2:]
		tissue_rpkm = {}
		assert len(rpkms) == len (cols)
		for i in range(len(rpkms)):
			sampleid = cols[i]
			tissue = sampid_to_smts[sampleid]
			if tissue not in tissue_rpkm:
				tissue_rpkm[tissue] = {}
			tissue_rpkm[tissue][sampleid] = rpkms[i]
		for tissue in tissue_rpkm:
			if tissue not in files:
				files[tissue] = open(tissue + ".bed", "w")
			files[tissue].write(trans_to_cord[gene] + '\t' + 
					    '\t'.join([tissue_rpkm[tissue][sampleid] for sampleid in tissue_rpkm[tissue]]) +
					    '\n')
for f in files:
	f.close()
