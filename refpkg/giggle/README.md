[![Build Status](https://travis-ci.org/ryanlayer/giggle.svg?branch=master)](https://travis-ci.org/ryanlayer/giggle)

<img src="https://raw.githubusercontent.com/ryanlayer/giggle/master/img/logo.png" width="300"/>

GIGGLE is a genomics search engine that identifies and ranks the
significance of shared genomic loci between query features and thousands of
genome interval files.

For questions and discussion about GIGGLE please visit/join the mailing list:
https://groups.google.com/d/forum/giggle-discuss

For more information about how GIGGLE works, please read the manuscript in Nature Methods:
https://www.nature.com/articles/nmeth.4556

Or watch a presentation about GIGGLE on YouTube (14m 37s)
[![GIGGLE](https://img.youtube.com/vi/yw8H7PhtZoA/0.jpg)](https://www.youtube.com/watch?v=yw8H7PhtZoA)


## Usage

GIGGLE has two high-level functions:
* `index` creates an index from a directory of bgzipped annotations (BED files
  or VCF files)
* `search` takes a region or a file of regions and searches them against an
  index

```
giggle, v0.6.3
usage:   giggle <command> [options]
     index     Create an index
     search    Search an index
```

### Indexing

    giggle, v0.6.3
    usage:   giggle index -i <input files> -o <output dir> -f
             options:
                 -s  Files are sorted
                 -i  Files to index (e.g. data/*.gz)
                 -o  Index output directory
                 -f  For reindex if output directory exists

### Searching

    giggle, v0.6.3
    usage:   giggle search -i <index directory> [options]
         options:
             -i giggle index directory
             -r <regions (CSV)>
             -q <query file>
             -o give results per record in the query file (omits empty results)
             -c give counts by indexed file
             -s give significance by indexed file (requires query file)
             -v give full record results
             -f print results for files that match a pattern (regex CSV)
             -g genome size for significance testing (default 3095677412)
             -l list the files in the index  

### Example

To demonstrate GIGGLE indexing and searching, we will curate and query a genome
repeat reference dataset. 

This example will use `gargs` from `https://github.com/brentp/gargs`

    wget -O gargs https://github.com/brentp/gargs/releases/download/v0.3.8/gargs_linux
    chmod +x gargs

or if you are on a Mac

    wget -O gargs https://github.com/brentp/gargs/releases/download/v0.3.8/gargs_darwin
    chmod +x gargs

This reference will be based on the following annotation from the UCSC genome
browser:
* Repeat Masker
* Segmental Duplications 
* Microsatellites 
* Simple Repeats

UCSC stores data as tables, and the relevant columns vary between files, so we
must take some care in curating the data.
    
    mkdir repeat
    url="http://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/rmsk.txt.gz"
    curl -s $url | gunzip -c | cut -f 6,7,8,11,12,13 > repeat/rmsk.bed

    url="http://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/simpleRepeat.txt.gz"
    curl -s $url | gunzip -c | cut -f 2,3,4,17 > repeat/simpleRepeat.bed

    url="http://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/microsat.txt.gz"
    curl -s $url | gunzip -c | cut -f 2,3,4 > repeat/microsat.bed

    url="http://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/genomicSuperDups.txt.gz"
    curl -s $url | gunzip -c | cut -f 2,3,4,5 > repeat/genomicSuperDups.bed

    url="http://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/chainSelf.txt.gz"
    curl -s $url | gunzip -c  | cut -f 3,5,6,7,10,11 > repeat/chainSelf.bed
 
Once all of the data are in bed files in the `repeat` directory, we sort and bgzip the
files then index.

    mkdir repeat_sort
    giggle/scripts/sort_bed "repeat/*.bed" repeat_sort 4
    giggle index -i "repeat_sort/*gz" -o repeat_sort_b -f -s

A GIGGLE index can be queried and the output formatted in a variety of ways.
The most basic is to search a single interval and get the number of overlaps for 
each database file.


    giggle search -i repeat_sort_b -r 1:200457776-200457776

    #repeat_sort/chainSelf.bed.gz   size:1058543    overlaps:0
    #repeat_sort/genomicSuperDups.bed.gz    size:51599  overlaps:0
    #repeat_sort/microsat.bed.gz    size:41572  overlaps:0
    #repeat_sort/rmsk.bed.gz    size:5298130    overlaps:1
    #repeat_sort/simpleRepeat.bed.gz    size:962714 overlaps:0

To search only a subset of database files use the `-f` option, which takes a 
comma separated list of regular expressions. Only those database files that
match one of the regular expressions will be considered.

    giggle search -i repeat_sort_b -r 1:200457776-200457776 -f rmsk,simple

    #repeat_sort/rmsk.bed.gz    size:5298130    overlaps:1
    #repeat_sort/simpleRepeat.bed.gz    size:962714 overlaps:0

To retrieve the original records for each overlap, use the `-v` option.  This
is useful for detailed filtering and summaries.

    giggle search -i repeat_sort_b -r 1:200457776-200457776 -f rmsk,simple -v

    chr1    200457488   200457811   L2a LINE    L2  repeat_sort/rmsk.bed.gz

GIGGLE also accepts query files in either `bed.gz` or `vcf.gz` formats. When a 
query file is given all of the above options are valid. In addition, GIGGLE can
perform statistical tests between the query file and each database file using the 
`-s` option. These tests include the:
* odds ratio that estimates the enrichment of observed v. expected
* the Fisher's two tailed, left tailed, and right tailed tests that estimate p-values
* the GIGGLE combo score that combines the odds ratio and Fisher's two tailed tests

```
giggle search -i repeat_sort_b -q bed.bed.gz -s
#file                             file_size   overlaps              odds_ratio         fishers_two_tail       fishers_left_tail   fishers_rigth_tail       combo_score
repeat_sort/chainSelf.bed.gz        1058543     410556  1.4753028915472957e-10  9.9350208733579337e-201 9.9350208733579337e-201   1                        0
repeat_sort/genomicSuperDups.bed.gz   51599      31434  1.2707669134541957      7.63990775297919e-91    1                         3.8346413983550421e-91  31.153365312726858
repeat_sort/microsat.bed.gz           41572      12320  1.0772570268110226      3.7685521408880773e-10  0.99999999982079036       1.9339512340673677e-10   1.0117655469639908
repeat_sort/rmsk.bed.gz             5298130    1599051  1.0068476528982573e-10  3.9822462776218267e-200 2.2545997326354446e-200   1                        0
repeat_sort/simpleRepeat.bed.gz      962714     308460  1.0977857352008243e-10  4.4719559780220419e-201 4.4719559780220419e-201   1                        0
```

Original records can also be retrieved and grouped by query interval with the `-v -o` options.

    giggle search -i repeat_sort_b -q ALL.wgs.integrated_sv_map_v2.20130502.svs.genotypes.vcf.gz -v -o | head

    ####1     645710  ALU_umary_ALU_2 A       <INS:ME:ALU>    .       .       TSD=null;SVTYPE=ALU;MEINFO=AluYa4_5,1,223,-;SVLEN=222;CS=ALU_umary;AC=35;AF=0.00698882;NS=2504;AN=5008;EAS_AF=0.0069;EUR_AF=0.0189;AFR_AF=0;AMR_AF=0.0072;SAS_AF=0.0041;SITEPOST=0.9998
    chr1    18392   804926  chr19   60000   244029  repeat_sort/chainSelf.bed.gz
    chr1    70007   667633  chr6    60000   147635  repeat_sort/chainSelf.bed.gz
    chr1    114546  672699  chr1    260861  676115  repeat_sort/chainSelf.bed.gz
    chr1    126221  713929  chr7    56425213        56942114        repeat_sort/chainSelf.bed.gz
    chr1    130987  703609  chr7    55804351        56479927        repeat_sort/chainSelf.bed.gz
    chr1    131111  751608  chr1    222641449       224228809       repeat_sort/chainSelf.bed.gz
    chr1    227417  706973  chr5    14937   199599  repeat_sort/chainSelf.bed.gz
    chr1    228316  648877  chr6    60000   193267  repeat_sort/chainSelf.bed.gz
    chr1    230805  804926  chr16   60000   164433  repeat_sort/chainSelf.bed.gz

## Building

### Dependencies
From a fresh install of Ubuntu, the following steps should provide all the
required dependencies.

    sudo apt install gcc make autoconf zlib1g-dev libbz2-dev libcurl4-openssl-dev libssl-dev ruby
    
### Giggle command line interface

    git clone https://github.com/ryanlayer/giggle.git
    cd giggle
    make
    export GIGGLE_ROOT=`pwd`
    cd ..

### Run tests

The first set of tests require bedtools to be in your path.

    sudo apt install g++ python

    git clone https://github.com/arq5x/bedtools2.git
    cd bedtools2
    make
    cd bin
    export PATH=$PATH:`pwd`
    cd ../..
    
Now run the tests

    cd $GIGGLE_ROOT/test/func
    ./giggle_tests.sh
    cd ../unit
    make
    cd ../../..
    
### Hosted data and services

#### Data
Roadmap Epigenomics:  https://s3.amazonaws.com/layerlab/giggle/roadmap/roadmap_sort.tar.gz

UCSC Genome browser:  https://s3.amazonaws.com/layerlab/giggle/ucsc/ucscweb_sort.tar.gz

Fantom5:  https://s3.amazonaws.com/layerlab/giggle/fantom/fantom_sort.tar.gz

#### Interactive heatmap

http://ryanlayer.github.io/giggle/index.html?primary_index=stix.colorado.edu/rme&ucsc_index=stix.colorado.edu/ucsc

### Web server (optional)
This is based on [libmicrohttpd](http://www.gnu.org/software/libmicrohttpd/)

    mkdir -p $HOME/usr/local/
    wget http://ftpmirror.gnu.org/libmicrohttpd/libmicrohttpd-0.9.46.tar.gz
    tar zxvf libmicrohttpd-0.9.46.tar.gz
    cd libmicrohttpd-0.9.46
    ./configure --prefix=$HOME/usr/local/
    make
    make install

    export LD_LIBRARY_PATH=$HOME/usr/local/lib/

    cd ..

    sudo apt install libtool

    wget https://github.com/json-c/json-c/archive/json-c-0.12.1-20160607.tar.gz
    tar xvf json-c-0.12.1-20160607.tar.gz  
    cd json-c-json-c-0.12.1-20160607
    ./configure --prefix=$HOME/usr/local/
    make
    make install

    cd $GIGGLE_ROOT
    make
    make server
    cd ..
    
To host the site shown in Supplemental Figure 3, you will need host web servers for both 
the Roadmap Epigenomics data and the UCSC data. Here we will run both servers on the same
host from ports `8080` and `8081` and access the web services using `localhost`, but these 
are general steps and apply to many other configurations including hosting the data sets 
on different servers.
     
    wget https://s3.amazonaws.com/layerlab/giggle/roadmap/roadmap_sort.tar.gz
    tar -zxvf roadmap_sort.tar.gz
    
    # NOTE, if the following command gives "Too many open files" try:
    ulimit -Sn 16384
    $GIGGLE_ROOT/bin/giggle index -s -f \
        -i "roadmap_sort/*gz" \
        -o roadmap_sort_b 
        
    wget https://s3.amazonaws.com/layerlab/giggle/ucsc/ucscweb_sort.tar.gz
    tar -zxvf ucscweb_sort.tar.gz
    
    $GIGGLE_ROOT/bin/giggle index -s -f \
        -i "ucscweb_sort/*gz" \
        -o ucscweb_sort_b
    
Start a web server for each index. 

    $GIGGLE_ROOT/bin/server_enrichment -i roadmap_sort_b/ -u /tmp/ -d $GIGGLE_ROOT/examples/rme/data_def.json -p 8080 &
    $GIGGLE_ROOT/bin/server_enrichment -i ucscweb_sort_b/ -u /tmp/ -d $GIGGLE_ROOT/examples/ucsc/data_def.json -p 8081 &

    If you get Access-Control-Allow-Origin errors, then pass the `-a` option to `server_enrichment`

Pass these two services to the web interface through URL arguments:

    http://ryanlayer.github.io/giggle/index.html?primary_index=localhost:8080&ucsc_index=localhost:8081
   
These data are also being served here:

http://ryanlayer.github.io/giggle/index.html?primary_index=ec2-54-227-176-15.compute-1.amazonaws.com/rme&ucsc_index=ec2-54-227-176-15.compute-1.amazonaws.com/ucsc

## Example analysis

**NOTE:** Index files and query files MUST be bgzipped (https://github.com/samtools/htslib, https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3042176/).

### Roadmap Epigenomics

    # details of how to recreate the data at 
    # https://github.com/ryanlayer/giggle/blob/master/examples/rme/README.md
    wget https://s3.amazonaws.com/layerlab/giggle/roadmap/roadmap_sort.tar.gz
    tar -zxvf roadmap_sort.tar.gz
    
    # NOTE, if the following command gives "Too many open files" try:
    # ulimit -Sn 16384
    $GIGGLE_ROOT/bin/giggle index -s -f \
        -i "roadmap_sort/*gz" \
        -o roadmap_sort_b 

    wget ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM1218nnn/GSM1218850/suppl/GSM1218850_MB135DMMD.peak.txt.gz
    # take the just the top peaks
    zcat GSM1218850_MB135DMMD.peak.txt.gz \
    | awk '$8>100' \
    | cut -f1,2,3 \
    | $GIGGLE_ROOT/lib/htslib/bgzip -c \
    > GSM1218850_MB135DMMD.peak.q100.bed.gz

    # List files in the index
    $GIGGLE_ROOT/bin/giggle search -l \
        -i roadmap_sort_b/ 

    # Search
    $GIGGLE_ROOT/bin/giggle search -s \
        -i roadmap_sort_b/ \
        -q GSM1218850_MB135DMMD.peak.q100.bed.gz \
    > GSM1218850_MB135DMMD.peak.q100.bed.gz.result

    
    # Plot
    sudo apt install python python-pip python-tk
    pip install matplotlib
    $GIGGLE_ROOT/scripts/giggle_heat_map.py \
        -s $GIGGLE_ROOT/examples/rme/states.txt \
        -c $GIGGLE_ROOT/examples/rme/EDACC_NAME.txt \
        -i GSM1218850_MB135DMMD.peak.q100.bed.gz.result \
        -o GSM1218850_MB135DMMD.peak.q100.bed.gz.result.3x11.pdf \
        -n $GIGGLE_ROOT/examples/rme/new_groups.txt \
        --x_size 3 \
        --y_size 11 \
        --stat combo \
        --ablines 15,26,31,43,52,60,72,82,87,89,93,101,103,116,120,122,127 \
        --state_names $GIGGLE_ROOT/examples/rme/short_states.txt \
        --group_names $GIGGLE_ROOT/examples/rme/new_groups_names.txt


## APIs

### [Python](https://github.com/brentp/python-giggle) by Brent Pedersen

```
from giggle import Giggle
index = Giggle('existing-index-dir') # or Giggle.create('new-index-dir', 'files/*.bed')
print(index.files)

result = index.query('chr1', 9999, 20000)
print(result.n_files)
print(result.n_total_hits) # integer number sum of hits across all files

print(result.n_hits(0)) # integer number of hits for the 0th file...

for hit in result[0]:
    print(hit) # hit is a string
```

#### Installation

make sure you have `liz`, `libcurl`, `libcrypto`, `libbz2` and `liblzma` installed in the appropriate
place on your system.
```
git clone --recursive https://github.com/brentp/python-giggle
cd python-giggle
python setup.py test
python setup.py install
```

### [Go](https://github.com/brentp/go-giggle) by Brent Pedersen

[![GoDoc](https://godoc.org/github.com/brentp/go-giggle?status.png)](https://godoc.org/github.com/brentp/go-giggle)

```Go

import (
    giggle "github.com/brentp/go-giggle"
    "fmt"
) 

func main() {

    index := giggle.Open("/path/to/index")
    res := index.Query("1", 565657, 567999)

    // all files in the index
    index.Files()

    // int showing total count
    res.TotalHits()

    // []uint32 giving number of hits for each file
    res.Hits()

    var lines []string
    # access results by index of file.
    lines = res.Of(0)
    fmt.Println(strings.Join(lines, "\n"))
    lines = res.Of(1)
}
```

## Docker

### [giggle-docker](https://github.com/kubor/giggle-docker) by Ryuichi Kubo

```
docker run kubor/giggle-docker giggle -h
```

```
Unknown command
giggle, v0.6.3
usage:   giggle <command> [options]
         index     Create an index
                  search    Search an index
```
