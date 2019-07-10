rsync -a -P rsync://hgdownload.cse.ucsc.edu/goldenPath/hg19/database ./

mkdir parsed_tracks

ls database/*sql \
| xargs \
    -I{} \
    -P 30 \
    bash -c "python parse_sql.py {} parsed_tracks/"

mkdir parsed_tracks_sorted

~/src/giggle/scripts/sort_bed "parsed_tracks/*gz" parsed_tracks_sorted 30

time ~/src/giggle/bin/giggle \
    index \
    -i "parsed_tracks_sorted/*gz" \
    -o parsed_tracks_sorted_b \
    -s \
    -f
Indexed 6980993757 intervals.

real    268m46.033s
user    245m44.262s
sys     11m51.567s
