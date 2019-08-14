mkdir orig
cd orig
wget -r -nH -nd -np -R index.html* http://noble.gs.washington.edu/proj/encyclopedia/interpreted/

cd ..
mkdir split

python ~/src/giggle/examples/segway/rename.py "orig/*gz" split

ls *bed | xargs -I{} -P 20 bash -c "bgzip {}"

cd ..

mkdir split_sort

~/src/giggle/scripts/sort_bed "split/*gz" split_sort 30

time ~/src/giggle/bin/giggle index -i "split_sort/*gz" -o split_sort_b -f -s

time ~/src/giggle/bin/giggle index -i "split_sort/*gz" -o split_sort_b -f -s

Indexed 1013379899 intervals.

real    23m28.112s
user    22m50.339s
sys     0m36.730s
