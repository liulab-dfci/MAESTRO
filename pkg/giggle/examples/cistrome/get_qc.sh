#!/bin/bash

echo -e "id\tmap\tpeaks\tfastqc\tfrip\tpbc\tmotif_judge\tdhs"

for i in `cat TF_human_data_information.txt | cut -f 1`; do
    echo -en "$i\t"
    curl -s http://dc2.cistrome.org/api/inspector?id=$i \
    | jq -cr '.qc.judge | [.map, .peaks, .fastqc, .frip, .pbc, .motif_judge, .dhs ]' \
    | jq -cr 'map(tostring)' \
    | jq -cr 'join("\t")'
done
