rule scrna_merge:
    input:
        rawmtx = expand("Result/STAR/{sample}/{sample}Solo.out/Gene/raw/matrix.mtx", sample=ALL_SAMPLES),
        features = expand("Result/STAR/{sample}/{sample}Solo.out/Gene/raw/features.tsv", sample=ALL_SAMPLES)[1],
        barcodes = expand("Result/STAR/{sample}/{sample}Solo.out/Gene/raw/barcodes.tsv", sample=ALL_SAMPLES)
    output:
        mergedmtx = "Result/%s/rawmatrix/matrix.mtx" % config["mergedname"],
        mergedfeatures = "Result/%s/rawmatrix/features.tsv" % config["mergedname"],
        mergedbarcodes = "Result/%s/rawmatrix/barcodes.tsv" % config["mergedname"],
        sampleannotation = "Result/%s/BarcodeAnnotation.tsv" % config["mergedname"]
    params:
        samplenames = ALL_SAMPLES
    shell:
        """
        numbarcodes=$(wc -l {input.barcodes} | awk 'NR==1{{printf $1}} NR>1{{printf ","$1}}')
        totalmtx=$(wc -l {input.rawmtx} | tail -1 | awk '{{print $1}}')
        awk -v numbars=$numbarcodes -v totalmtx=$totalmtx '\
          BEGIN {{lenbars=split(numbars,bararray,","); fileindex=1; extrabarcount=0}} \
          NR==1 || NR==2{{print; next}} \
          NR==3{{print $1,bararray[lenbars],totalmtx-(3*(lenbars-1)); next}} \
          FNR==1{{extrabarcount+=bararray[fileindex]; fileindex++; next}} \
          FNR==2 || FNR==3{{next}} \
          {{print $1,$2+extrabarcount,$3}}' {input.rawmtx} > {output.mergedmtx}
        cp {input.features} {output.mergedfeatures}
        awk -v names=$(echo {params.samplenames} | tr " " ",") '\
          BEGIN {{split(names,samplenames,","); sampleindex=0}} \
          FNR==1 {{sampleindex++}} \
          {{print samplenames[sampleindex]"_"$1}}
          {{print samplenames[sampleindex]"\t"samplenames[sampleindex]"_"$1 > "{output.sampleannotation}"}}' \
          {input.barcodes} > {output.mergedbarcodes}
        """
