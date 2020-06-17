import sys
import os

fail_flag = 99
tol = 5e-2


def fail_msg(msg):
    sys.stderr.write(msg + '\n')
    sys.exit(fail_flag)


def isfloat(value):
  try:
    float(value)
    return True
  except ValueError:
    return False


def diff_float(v1, v2):
    denom = abs(v1) + abs(v2)
    if denom < tol: return False
    
    diff = abs(v1-v2)/denom
    return (diff > tol)



def compare_file(f1, f2):
    
    msg = 'File not equal between \"' + os.path.basename(f1) + '\" \"' + os.path.basename(f2) + '\"'
    
    if not os.path.exists(f1): fail_msg('File not exists \"' + f1 + '\"')
    if not os.path.exists(f2): fail_msg('File not exists \"' + f2 + '\"')
    
    fin1 = open(f1)
    fin2 = open(f2)
    
    for l1 in fin1:
        fields1 = l1.strip().split('\t')
        N = len(fields1)
        
        fields2 = fin2.readline().strip().split('\t')
        if N != len(fields2): fail_msg(msg)
        
        for i in range(N):
            v1 = fields1[i]
            v2 = fields2[i]
            
            float_flag = isfloat(v1)
            
            if float_flag != isfloat(v2): fail_msg(msg)
            
            if float_flag:
                if diff_float(float(v1), float(v2)):
                    fail_msg(msg)
            else:
                if v1 != v2:
                    fail_msg(msg)
    
    fin2.close()
    fin1.close()


def main():
    argv = sys.argv
    
    if len(argv) < 3: fail_msg('Not enough parameters')
    
    f1 = argv[1]
    f2 = argv[2]
    
    compare_file(f1, f2)
    compare_file(f1 + '.t', f2 + '.t')
    #compare_file(f1 + '.FDR', f2 + '.FDR')
    
    sys.stdout.write('\"' + os.path.basename(f1) + '\" \"' + os.path.basename(f2) + '\" are the same.\n')

    sys.exit(0)

if __name__ == '__main__':
    main()
