import os

def get_fastqfile(fastqpath):
    files = os.listdir(fastqpath)
    fastq_1 = []
    fastq_2 = []
    fastqs = []
    for f in files:
        if f.endswith('.fastq'):
            fastqs.append(f)
            if f.endswith('_1.fastq'):
                fastq_1.append(fastqpath + f)
            elif f.endswith('_2.fastq'):
                fastq_2.append(fastqpath + f)

    fastq_1 = sorted(fastq_1)
    fastq_2 = sorted(fastq_2)
    fastqs = sorted(fastqs)
    
    fastqstr = ''
    if len(fastq_1) != 0:
        if len(fastq_1) == len(fastq_2) and len(fastqs) == len(fastq_1) + len(fastq_2):
            fastqstr = ','.join(fastq_1) + ' ' + ','.join(fastq_2)
        else:
            print("Invalid fastq files!")
    else:
        if len(fastq_1) == len(fastq_2) and len(fastqs) != 0:
            fastqstr = ','.join(fastqs)
        else:
            print("Invalid fastq files!")

    return(fastqstr)

def get_fastqid(fastqpath):
    files = os.listdir(fastqpath)
    fastq_1 = []
    fastq_2 = []
    fastqs = []
    for f in files:
        if f.endswith('.fastq'):
            fastqs.append(f)
            if f.endswith('_1.fastq'):
                fastq_1.append(f)
            elif f.endswith('_2.fastq'):
                fastq_2.append(f)

    fastq_1 = sorted(fastq_1)
    fastq_2 = sorted(fastq_2)
    fastqs = sorted(fastqs)
    
    fastqidstr = ''
    if len(fastq_1) != 0:
        if len(fastq_1) == len(fastq_2) and len(fastqs) == len(fastq_1) + len(fastq_2):
            samplelist = ['ID:'+ i[0:len(i)-8] for i in fastq_1]
            fastqidstr = ' , '.join(samplelist)
        else:
            print("Invalid fastq files!")
    else:
        if len(fastq_1) == len(fastq_2) and len(fastqs) != 0:
            samplelist = ['ID:' + i[0:len(i)-6] for i in fastqs]
            fastqidstr = ' , '.join(fastqs)
        else:
            print("Invalid fastq files!")

    return(fastqidstr)

def get_fastqlist(fastqpath):
    files = os.listdir(fastqpath)
    fastq_1 = []
    fastq_2 = []
    fastqs = []
    for f in files:
        if f.endswith('.fastq'):
            fastqs.append(f)
            if f.endswith('_1.fastq'):
                fastq_1.append(f)
            elif f.endswith('_2.fastq'):
                fastq_2.append(f)

    fastq_1 = sorted(fastq_1)
    fastq_2 = sorted(fastq_2)
    fastqs = sorted(fastqs)
    
    if len(fastq_1) != 0:
        if len(fastq_1) == len(fastq_2) and len(fastqs) == len(fastq_1) + len(fastq_2):
            samplelist = [i[0:len(i)-8] for i in fastq_1]
        else:
            print("Invalid fastq files!")
    else:
        if len(fastq_1) == len(fastq_2) and len(fastqs) != 0:
            samplelist = [i[0:len(i)-6] for i in fastqs]
        else:
            print("Invalid fastq files!")

    return(samplelist)

def get_bamfile(bampath):
    files = os.listdir(bampath)
    bams = []
    for f in files:
        if f.endswith('Aligned.sortedByReads.out.bam'):
            bams.append(bampath + f)
    bams = sorted(bams)

    bam_str = ' '.join(bams)
    return(bam_str)
