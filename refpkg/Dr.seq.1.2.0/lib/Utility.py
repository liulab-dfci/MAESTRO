#!/usr/bin/env python
"""

Function declare:
 
def CMD                  (cmd)
def sp                   (cmd)
def sperr                (cmd)
def raise_error          ()
def detect_memory        ()
def pdf_name             (input_name)
def wlog                 (message,logfile)
def ewlog                (message,logfile)
def rwlog                (message,logfile)
def readAnnotation       (annotation)
def textformat           (inp)
def createDIR            (dirname)
def strlatexformat       (instr)
def strdis               (str1,str2)
def sample_down_transform_sam (samfile,outbed,sampledownsam,sampledownbed,sample_reads)
def transform_refgene    (refgene,ttsdis,outname)
def reform_barcode_fastq (fq,reformtxt,cbL,umiL)
def combine_reads        (barcodeF,cdsF,utr3F,utr5F,symbolF,ttsdisF,outF,dup_measure)
def generate_matrix      (refgene,inputbed,ttsdis,qcmatfull,qcmat,expmat,coverGNcutoff,umidis1)
"""
import subprocess
import sys
import os
import math
import random
import string
def CMD(cmd):
    os.system(cmd)

def sp(cmd):
    '''
    Call shell cmd or software and return its stdout
    '''
    a=subprocess.Popen(cmd, stdout=subprocess.PIPE, shell='TRUE')
    ac = a.communicate()
    return ac
def sperr(cmd):
    '''
    Call shell cmd or software and return its stdout
    '''
    a=subprocess.Popen(cmd, stderr=subprocess.PIPE, shell='TRUE')
    ac = a.communicate()
    return ac
   
def raise_error():
    '''
    Raise an error messgae and exit
    '''
    print('error occurs, check log file~!')
    sys.exit(1)


def detect_memory():
    meminfo={}#OrderedDict()
    try:
        with open('/proc/meminfo') as f:
            for line in f:
                meminfo[line.split(':')[0].strip()] = line.split(':')[1].strip()
        totalM = meminfo['MemTotal'].split()
        #freeM = meminfo['MemFree'].split()
        if totalM[1].lower() == "kb":
            try:
                totalM_G = int(totalM[0])/1e6
                return totalM_G
            except:
                return 'NA'
        else:
            return 'NA'    
    except:
        return 'NA'


def pdf_name(input_name):
    '''
    Change filename to pdf file name
    '''
    outputname = "_".join(input_name.split('.')[:-1])+".pdf"
    return outputname
    
def wlog(message,logfile):
    '''
    print a message and write the message to logfile
    '''
    print(message)
    os.system('echo "%s " >> %s'%(message,logfile))
    
def ewlog(message,logfile):
    '''
    print an error message and write the error message to logfile
    then exit Dr.seq
    error messages start with [ERROR]
    '''
    print("[ERROR] %s "%(message))
    os.system('echo "[ERROR] %s " >> %s'%(message,logfile))
    raise_error()
    
def rwlog(cmd,logfile) :
    '''
    print an (shell) command line and write the command line to logfile
    then conduct the command line
    command lines start with [CMD]
    '''
    print("[CMD] %s "%(cmd))
    os.system('echo "[CMD] %s " >> %s'%(cmd,logfile))
    CMD(cmd)
   
    
def readAnnotation(annotation):
    '''
    read full annotation file and output as a dictionary 
    file format is fixed to UCSC full annotation format
    '''
    inf = open(annotation)
    outdict = {}
    for line in inf:
        ll = line.split()
        outdict[ll[1]] = ll[12]
    return outdict
    
     
def textformat(inp):
    '''
    transfer 1000000 to  1,000,000 for better visualization in output report
    '''
    o = ''
    comma = 0
    for i in (inp[::-1]):
        comma += 1
        o += i
        if comma%3 == 0:
            o += ','
    return o[::-1].strip(',')

def createDIR(dirname):
    '''
    check dir name and create new dir
    '''
    if not os.path.isdir(dirname):
        os.system('mkdir %s'%(dirname))

def strlatexformat(instr):
    outstr = instr.replace('_','\_')
    return(outstr)
    
def sample_down_transform_sam(samfile,outBed,sampledown_sam,sampledown_bed,sample_reads,q30filter):
    '''
    transform aligned samfile to bed file, and randomly sample N reads to a small samfile (for reads/bulk QC)
    '''
    q30 = int(q30filter)
    inf = open(samfile)
    outbed = open(outBed,'w')
    outSDbed = open(sampledown_bed,'w')
    outSDsam = open(sampledown_sam,'w')

   # totalN = int(sp('wc -l %s'%(samfile))[0].split()[0]) - headline_N
   # p = 10000 * (int(sample_reads)*1.0/totalN)
    totalN = 0
    for line in inf:
        if line.startswith('@'):
            continue
        ll = line.strip().split("\t")
        chrom = ll[2]
        start = int(ll[3])-1
        seqlen = len(ll[9])
        end = start + seqlen
        txname = ll[0]
        mapQ = ll[4]
        flag = bin(int(ll[1]))[2:]
        if len(flag) < 5:
            if len(flag) < 3:
                strand = "+"
            else:
                if flag[-3] == "1":
                    continue
                else:
                    strand = "+"
        else:
            if flag[-3] == "1":
                continue
            elif flag[-5] == "1":
                strand = "-"
            else:
                strand = "+"
        if q30 == 1 and int(mapQ) < 30:
            continue
        totalN += 1
        newll = [chrom,start,end,txname,'255',strand]
        outbed.write("\t".join(map(str,newll))+"\n")
    outbed.close()

    inf.seek(0)
    p = 10000 * (int(sample_reads)*1.0/totalN)

    for line in inf:
        if line.startswith('@'):
            outSDsam.write(line)
            continue
        ll = line.strip().split("\t")
        chrom = ll[2]
        start = int(ll[3])-1
        seqlen = len(ll[9])
        end = start + seqlen
        txname = ll[0]
        mapQ = ll[4]
        flag = bin(int(ll[1]))[2:]
        if len(flag) < 5:
            if len(flag) < 3:
                strand = "+"
            else:
                if flag[-3] == "1":
                    continue
                else:
                    strand = "+"
        else:
            if flag[-3] == "1":
                continue
            elif flag[-5] == "1":
                strand = "-"
            else:
                strand = "+"
        if q30 == 1 and int(mapQ) < 30:
            continue
        if random.randint(1,10000) <= p:
            outSDsam.write(line)
            newll = [chrom,start,end,txname,'255',strand]
            outSDbed.write("\t".join(map(str,newll))+"\n")
    outSDbed.close()
    outSDsam.close()
    inf.close()

    
def transform_refgene(refgene,ttsdis,outname):
    '''
    transform full gene annotation downloaded from UCSC to different meterials used in Dr.seq
    including
    1. CDS exon 
    2. 5' utr exon
    3. 3' utr exon
    4. transcript bed
    5. gene symbol bed(merge transcripts belonging to same gene)
    6. TTS +- 400(default) bed
    7. full bed for RseQC input 
    '''
    ret_lst=[]
    inf = open(refgene)
    
    utr5_list = []
    utr3_list = []    
    CDSexon_list = []
    Exon_list = [] ### only consider exon, do not divide utr and cds, for genebody coverage
    transcript_list = []
    symbol_dict = {}    
    TTS400_list = []
    GBbin_list = []
    for line in inf:
        if line.startswith('#'):
            continue
        f = line.strip().split()
        chrom = f[2]
        strand = f[3]
        symbol = f[12]
        txStart = int(f[4])
        txEnd = int(f[5])
        txName = f[1]
        cdsStart = int(f[6])
        cdsEnd = int(f[7])
        exonCount = int(f[8])
        exonStarts = [ int(i) for i in f[9].strip(',').split(',') ]
        exonEnds = [ int(i) for i in f[10].strip(',').split(',') ]
        exonlength = []
        exondisTss = []
        for i in range(len(exonStarts)):
            exonlength.append(exonEnds[i] - exonStarts[i])
            exondisTss.append(exonStarts[i] - txStart)
        
        total_exon_length = sum(exonlength)*1.0
        bin_length = total_exon_length/100
        
        if strand == "+" :
            TSS = txStart
            TTS = txEnd
        else:
            TSS = txEnd
            TTS = txStart
        ### collect transcript info
        #refbed_list.append([chrom,txStart,txEnd,txName,'0',strand,cdsStart,cdsEnd,'0',exonCount,",".join(map(str,exonlength))+",",",".join(map(str,exondisTss))+","])
        transcript_list.append([chrom,txStart,txEnd,txName,symbol,strand])
        TTS400_list.append([chrom,max(0,TTS-int(ttsdis)),TTS+int(ttsdis)])
        mergeYes = 0
        if symbol not in symbol_dict:
            symbol_dict[symbol] = [[chrom,txStart,txEnd]]
        else:
            for mergeTX in symbol_dict[symbol]:

                if mergeTX[0] == chrom and mergeTX[1] < txEnd and mergeTX[2] > txStart:
                    mergeTX[1] = min(mergeTX[1],txStart)
                    mergeTX[2] = max(mergeTX[2],txEnd)
                    mergeYes = 1
                    break
                else:
                    mergeYes = 0
            if mergeYes == 0:
                symbol_dict[symbol].append([chrom,txStart,txEnd])
                    
        for st,end in zip(exonStarts,exonEnds):
            Exon_list.append([st,end])
            if st < cdsStart:
                utr_st = st
                utr_end = min(end,cdsStart)
                utr5_list.append([chrom,utr_st,utr_end])
                if cdsStart < end:
                    CDSexon_list.append([chrom,cdsStart,end])
            if end > cdsEnd:
                utr_st = max(st, cdsEnd)
                utr_end = end
                utr3_list.append([chrom,utr_st,utr_end])
                if st < cdsEnd:
                    CDSexon_list.append([chrom,st,cdsEnd])
            if st >= cdsStart and end <= cdsEnd:
                CDSexon_list.append([chrom,st,end])
        
        ### output 100 bins for exons
        exonnum = 0
        current_position = exonStarts[exonnum]
        for N in range(100):
            Nlen = math.floor(bin_length*(N+1))-math.floor(bin_length*N)
            if Nlen < exonEnds[exonnum] - current_position:
                # this exon is enough to support this bin
                newll = [chrom,int(current_position),int(current_position+Nlen),txName,symbol,strand,N]
                GBbin_list.append(newll)
                current_position = current_position + Nlen
            elif Nlen == exonEnds[exonnum] - current_position:
                newll = [chrom,int(current_position),int(current_position+Nlen),txName,symbol,strand,N]
                GBbin_list.append(newll)
                if N < 99:
                    exonnum += 1
                    current_position = exonStarts[exonnum]
            else:
                # this exon is not enough, add parts of next exon
                newll = [chrom,int(current_position),int(exonEnds[exonnum]),txName,symbol,strand,N]
                GBbin_list.append(newll)
                leftlen = Nlen - (exonEnds[exonnum] - current_position)
                while 1:
                    exonnum += 1
                    current_position = exonStarts[exonnum]
                    if leftlen < exonEnds[exonnum] - current_position:
                        newll = [chrom,int(current_position),int(current_position+leftlen),txName,symbol,strand,N]
                        GBbin_list.append(newll)
                        current_position = current_position + leftlen
                        break
                    elif leftlen == exonEnds[exonnum] - current_position:
                        newll = [chrom,int(current_position),int(current_position+leftlen),txName,symbol,strand,N]
                        GBbin_list.append(newll)
                        if N < 99:
                            exonnum += 1
                            current_position = exonStarts[exonnum]
                        break
                    else:
                        newll = [chrom,int(current_position),int(exonEnds[exonnum]),txName,symbol,strand,N]
                        GBbin_list.append(newll)
                        leftlen = leftlen - (exonEnds[exonnum] - current_position)
                    
                    
    ### output tx information
    outf = open(outname + '_gene_anno_binexon.bed','w')
    for i in GBbin_list:
        outf.write("\t".join(map(str,i))+"\n")
    outf.close()        
    
    outf = open(outname+'_gene_anno_transcript.bed','w')
    for i in transcript_list:
        outf.write("\t".join(map(str,i))+"\n")
    outf.close()
    
    outf = open(outname+'_gene_anno_TTSdis.bed','w')
    for i in TTS400_list:
        outf.write("\t".join(map(str,i))+"\n")
    outf.close()

    
    outf = open(outname+'_gene_anno_5utr.bed','w')
    for i in utr5_list:
        outf.write("\t".join(map(str,i))+"\n")
    outf.close()
    
    outf = open(outname+'_gene_anno_3utr.bed','w')
    for i in utr3_list:
        outf.write("\t".join(map(str,i))+"\n")
    outf.close()

    outf = open(outname+'_gene_anno_cds.bed','w')
    for i in CDSexon_list:
        outf.write("\t".join(map(str,i))+"\n")
    outf.close()    
        
    outf = open(outname+'_gene_anno_symbol.bed','w')
    for sym in symbol_dict:
        count=1
        for region in symbol_dict[sym]:
            newll = region + [sym]
            outf.write("\t".join(map(str,newll))+"\n")
            count += 1 
    outf.close()
           
def reform_barcode_fastq(fq,reformtxt,cbL,umiL):
    '''
    transform barcode fastq to another txt format [name,cell_barcode,umi] for following usage
    '''
    lastL = cbL+umiL
    inf = open(fq)
    outf = open(reformtxt,'w')
    count = 0
    for line in inf:
        count += 1
        if count%4 == 1:
            head = line.split()[0][1:]
        if count%4 == 2:
            seq = line.strip()
        if count%4 == 3:
            pass
        if count%4 == 0:
            if cbL == lastL:
                newll = [head,seq[:cbL],"NA"]
            else:
                newll = [head, seq[:cbL], seq[cbL:lastL]]
            outf.write("\t".join(newll)+"\n")
    
    outf.close()
    inf.close()    


def combine_reads(barcodeF,cdsF,utr3F,utr5F,symbolF,ttsdisF,outF,dup_measure):
    '''
    combine annotation information of all reads which is generate in previous step with bedtools
    '''
    barcode_file = open(barcodeF)
    cds_file = open(cdsF)
    utr3_file = open(utr3F)
    utr5_file = open(utr5F)
    symbol_file = open(symbolF)
    TTS400_file = open(ttsdisF)
    outf = open(outF,'w')

    last_read = ["NA"]*8
    next_read_sym = symbol_file.readline().split()
    for line in cds_file:
        current_read_cds = line.split()
        current_read_utr3 = utr3_file.readline().split()
        current_read_utr5 = utr5_file.readline().split()
        current_read_TTS400 = TTS400_file.readline().split()

        read_info = current_read_cds[:6]
        read_name = current_read_cds[3]
        cdsinfo = current_read_cds[6]
        utr3info = current_read_utr3[6]
        utr5info = current_read_utr5[6]
        TTS400info = current_read_TTS400[6]
        newll = []
        if read_name == last_read[3] :
            newll = last_read
        else:
            while(1):
                current_barcode = barcode_file.readline().split()
                if current_barcode[0] == read_name:
                    newll = read_info  + current_barcode[1:] #+ [cdsinfo,utr3info,utr5info]                 
                    break
        if newll == []:
            print('error in match barcode')
            sys.exit(1)
        last_read = newll[:8]
        if len(next_read_sym) > 3 and newll[3] == next_read_sym[3] :
            addsym_list = [next_read_sym[9]]
            while(1):
                next_read_sym = symbol_file.readline().split()
                if len(next_read_sym)>3 and newll[3] == next_read_sym[3]   :
                    if not next_read_sym[9] in addsym_list:
                        addsym_list.append(next_read_sym[9])
                else:
                    break
            addsym = ",".join(addsym_list)
        else:
            addsym = "NA"
                
        newll += [cdsinfo,utr3info,utr5info,TTS400info,addsym] 
        if int(dup_measure) == 1:
            newll[4] = "_".join([newll[7],newll[0],newll[5],newll[1]])
        elif int(dup_measure) == 2:
            newll[4] = newll[7]
        elif int(dup_measure) == 3:
            newll[4] = "_".join([newll[0],newll[5],newll[1]])
        else:
            newll[4] = "NA"                
            
        outf.write("\t".join(newll)+"\n")

    outf.close()

def strdis(str1,str2):
    if len(str1) != len(str2):
        print('umi have different length in different reads')
        sys.exit(1)
    diff = 0
    for i in range(len(str1)):
        if str1[i] != str2[i]:
            diff += 1
    return diff

def generate_matrix(refgene,inputbed,ttsdis,qcmatfull,qcmat,expmat,coverGNcutoff,umidis1):
    '''
    generate two matrix
    1. expression matrix whose row/column is corresponded to genes and cell_barcodes
    2. QC matrix whose row/column is corresponded to cell_barcodes and measurement(including total reads, #umi, #cds exon reads ,..., )
    '''

    inf = open(refgene)
    allgenes = []
    for line in inf:
        if line.startswith('#'):
            continue
        ll = line.strip().split()
        if not ll[12] in allgenes:
            allgenes.append(ll[12])
    inf.close()

    inf = open(inputbed)

    QCmat = {}
    Expmat = {}

    last_cell = "NA"
    last_umi = "NA"
    last_gname = []
    for line in inf:
        ll = line.strip().split()
        cell = ll[6]
        umi = ll[4]
        gname = ll[12].split(",") #['a','b']
    
        if cell not in Expmat:
            Expmat[cell] = {}
            QCmat[cell] = [0]*9# allreads,umi,cds,3utr,5utr,intron,intergenic,#greatertts400,covered gene number
    
        QCmat[cell][0] += 1

        share_gene = 0
        for g in gname:
            if g in last_gname:        
                share_gene = 1
        
        if cell == last_cell and umi == last_umi and umi != "NA" and  share_gene == 1 : 
            pass
        elif cell == last_cell and umi != "NA" and int(umidis1) == 1 and strdis(umi,lastumi) == 1 :
            pass
        else:
            if ttsdis ==1 and int(ll[11]) == 0:
                pass
            else:
                if int(ll[8]) > 0 or int(ll[9]) > 0 or int(ll[10]) > 0:
                    if ll[12] != "NA":
                        for target_gene in ll[12].split(","):
                            if target_gene not in Expmat[cell]:
                                Expmat[cell][target_gene] = 0
                            Expmat[cell][target_gene] += 1           
        
            QCmat[cell][1] += 1
            if ll[12] == "NA":
                QCmat[cell][6] += 1
            elif int(ll[8]) > 0:
                QCmat[cell][2] += 1
            elif int(ll[9]) > 0:
                QCmat[cell][3] += 1
            elif int(ll[10]) > 0:
                QCmat[cell][4] += 1
            else:
                QCmat[cell][5] += 1
            if int(ll[11]) > 0:
                QCmat[cell][7] += 1
        
        last_cell =  cell#ll[6]
        last_umi = umi#ll[4] 
        last_gname = gname
    
    inf.close()
    outf0 = open(qcmatfull,'w')
    outf1 = open(qcmat,'w')
    outf2 = open(expmat,'w')

    newll = ['cellname','allreads','umi','cds','utr3','utr5','intron','intergenic','awayTTS','coveredGN']
    outf0.write("\t".join(newll)+"\n")
    outf1.write("\t".join(newll)+"\n")

    newll = ['cellname'] + allgenes
    EXPmat = []
    EXPmat.append(newll)
#    outf2.write("\t".join(newll)+"\n")
    for cell in sorted(Expmat.keys()):
        coverN = len(Expmat[cell].keys())
        QCmat[cell][8] = coverN
        newllqc = [cell] + QCmat[cell]
        outf0.write("\t".join(map(str,newllqc))+"\n")
        if coverN >= int(coverGNcutoff):
            outf1.write("\t".join(map(str,newllqc))+"\n")
            explist = []
            for g in allgenes:
                if g in Expmat[cell]:
                    explist.append(Expmat[cell][g])
                else:
                    explist.append(0)
            newll = [cell] + explist
            EXPmat.append(newll)
#            outf2.write("\t".join(map(str,newll))+"\n")
    for i in range(len(EXPmat[0])):
        newll  =[]
        for j in range(len(EXPmat)):
            newll.append(EXPmat[j][i])
        outf2.write("\t".join(map(str,newll))+"\n")
    outf0.close()  
    outf1.close()
    outf2.close()
  
def readsqc(SDsamfile,outname):

    inf = open(SDsamfile)
    test_line_number = 0
    while 1:
        line = inf.readline()
        if line.strip() == "":
            continue
        if line.startswith("@"):
            continue
        ll = line.strip().split("\t")
        if len(ll) < 11:
            continue
        seq = ll[9]
        seqlen = len(seq)
        test_line_number += 1
        if test_line_number == 10:
            break
    inf.seek(0)    

    readsQuality = {}
    NVC = {}
    NVC['A'] = [0]*seqlen
    NVC['C'] = [0]*seqlen
    NVC['G'] = [0]*seqlen
    NVC['T'] = [0]*seqlen
    GCsummary = {}
    for i in range(seqlen+1):
        GCsummary[i] = 0
    for line in inf:
        if line.startswith("@"):
            continue
        if line.strip() == "":
            continue
        ll = line.split()
        if len(ll) < 11:
            continue
        try:
            flag = bin(int(ll[1]))[2:]
        except:
            continue
        if len(flag) < 5:
            if len(flag) < 3:
                strand = "+"
            else:
                if flag[-3] == "1":
                    continue
                else:
                    strand = "+"
        else:
            if flag[-3] == "1":
                continue
            elif flag[-5] == "1":
                strand = "-"
            else:
                strand = "+"

        if strand == "+":
            thisseq = ll[9].upper()
            thisQuality = ll[10]
        
        else:
            transtab = str.maketrans("ACGTNX","TGCANX")
            thisseq = ll[9].upper().translate(transtab)[::-1]
            thisQuality = ll[10][::-1]
        if len(thisseq) != seqlen:
            continue
        ### GC 
        GCnumber = 0
        for position in range(len(thisseq)):
            bp = thisseq[position].upper()
            qul = ord(thisQuality[position])-33
            if qul not in readsQuality:
                readsQuality[qul] = [0]*seqlen
            readsQuality[qul][position] += 1
            if bp in ['A','C','G','T']:
                NVC[bp][position] += 1
            if bp == "G" or bp == "C":
                GCnumber += 1
        
        #if not GCsummary.has_key(GCnumber):
        #    GCsummary[GCnumber] = 0
        GCsummary[GCnumber] += 1
        
    inf.close()
    
    GCoutf = open(outname+'_qcGC.txt','w')
    for i in range(seqlen+1):
        newll = [i,GCsummary[i],seqlen]
        GCoutf.write("\t".join(map(str,newll))+"\n")
    GCoutf.close()
    
    QULoutf = open(outname+'_qcQul.txt','w')
    for i in range(max(readsQuality.keys())+1):
        if i in readsQuality:
            newll = [i]+readsQuality[i]
        else:
            newll = [i] + [0]*seqlen
        QULoutf.write("\t".join(map(str,newll))+"\n")
    QULoutf.close()
    
    NVCoutf = open(outname+'_qcNVC.txt','w')
    for i in ['A','C','G','T']:
        newll = [i]+NVC[i]
        NVCoutf.write("\t".join(map(str,newll))+"\n")
    NVCoutf.close()
      

def GBcover(SDreads_on_gbbin,outname):
    inf= open(SDreads_on_gbbin)
    gbcount = [0]*100
    for line in inf:
        ll = line.split()
        if ll[5] == "-":
            gbcount[99-int(ll[6])]+=int(ll[7])
        else:
            gbcount[int(ll[6])]+=int(ll[7])
    inf.close()
    outf = open(outname+'_qcGBcover.txt','w')
    for i in range(len(gbcount)):
        newll = [i+1,gbcount[i]]
        outf.write("\t".join(map(str,newll))+"\n")
    outf.close()
            
       
