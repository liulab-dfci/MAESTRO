import sys
import re

N = int(sys.argv[1])

region = ''
tissue = ''

BP = []
tissues = {}
regions = []

for l in sys.stdin:
    if l[:2] == '##':
        A = l[2:].rstrip().split('\t')
        region = A[0] + ':' + A[1] + '-' + str(int(A[1]) + 1)
        regions.append(region)
        BP.append({})
    elif l[0] == '#':
        tissue = l.rstrip().split('/')[1].split('.')[0]
        tissues[tissue] = 1
    else:
        A = l.rstrip().split('\t')
        BP[-1][tissue] = A[3:]

print ' '.join(sorted(tissues.keys()))
i = 0;
for bp in BP:
    if len(bp) == 0:
        continue
    O = []
    for tissue in sorted(tissues.keys()):
        if tissue in bp:
            O.append(','.join([str(x) for x in bp[tissue]]))
        else:
            O.append(0)
    print ' '.join([regions[i]] + O)
    i+=1
    if i == N:
        break
