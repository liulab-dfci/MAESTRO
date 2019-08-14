import sys

if len(sys.argv) != 4:
    sys.stderr.write('usage:\t' + \
                     sys.argv[0] + \
                     ' <name2library file>' + \
                     ' <expression count matrix file>' + \
                     ' <out dir>\n')
    sys.exit(1)

name2library_file=sys.argv[1]
expression_count_matrix_file=sys.argv[2]
out_dir=sys.argv[3]

files={}

names = {}
for l in open(name2library_file, 'r'):
    A = l.rstrip().split('\t')
    names[A[1]] = A[0]

header = []
for l in open(expression_count_matrix_file, 'r'):
    A = l.rstrip().split()
    if A[0] == 'Id':
        header = A[1:]
        print len(header)
        0/1
    else:
        i = 0
        for a in A[1:]:
            if a != '0':
                if names[header[i]] not in files:
                    files[names[header[i]]] = \
                            open(out_dir + \
                                 '/' + \
                                 names[header[i]] + \
                                 '.bed',
                                 'w')

                
                files[names[header[i]]].write( \
                        A[0].replace(':','\t').replace('-','\t') + \
                        '\t' + a + '\n')
            i+=1
