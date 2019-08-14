Get data

    wget ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM1218nnn/GSM1218850/suppl/GSM1218850_MB135DMMD.peak.txt.gz

    zcat GSM1218850_MB135DMMD.peak.txt.gz \
    | awk '$8 >= 100'
    | $HTSLIB_ROOT/bgzip -c > GSM1218850_MB135DMMD.peak.q100.bed.gz
