#include <stdlib.h>
#include <stdio.h>
#include <err.h>
#include <string.h>

#include "giggle_index.h"
#include "wah.h"
#include "cache.h"
#include "file_read.h"
#include "kfunc.h"
#include "util.h"
#include "ll.h"

#define MAX(X, Y) (((X) > (Y)) ? (X) : (Y))




int main(int argc, char **argv)
{
    WAH_SIZE = 32;
    WAH_MAX_FILL_WORDS = (1<<(WAH_SIZE-1)) - 1;

    uint32_t num_chrms = 100;

    if ((argc != 4)) {
        errx(1,
             "usage:\t%s <input file> <index dir> <w|i>",
             argv[0]);
    }

    double genome_size =  3095677412.0;

    char *input_file = argv[1];
    char *index_dir = argv[2];
    char *i_type = argv[3];

    struct input_file *in_f = input_file_init(input_file);

    int chrm_len = 50;
    char *chrm = (char *)malloc(chrm_len*sizeof(char));
    uint32_t start, end;
    long offset;
    kstring_t line = {0, 0, NULL};

    struct giggle_index *gi;

    gi = giggle_load(index_dir,
                     uint64_t_ll_giggle_set_data_handler);

    uint32_t *file_counts = (uint32_t *)
            calloc(gi->file_idx->index->num, sizeof(uint32_t));

    uint32_t num_intervals = 0;
    double mean_interval_size = 0.0;
    while ( in_f->input_file_get_next_interval(in_f, 
                                               &chrm,
                                               &chrm_len,
                                               &start,
                                               &end,
                                               &offset,
                                               &line) >= 0 ) {
        num_intervals += 1;
        mean_interval_size += end - start;

        struct uint64_t_ll *R =
                (struct uint64_t_ll *)giggle_query_region(gi,
                                                          chrm,
                                                          start,
                                                          end);
        if (R != NULL) {
            struct uint64_t_ll_node *curr = R->head;

            while (curr != NULL) {
                /*
                struct file_id_offset_pair *fid_off = 
                    (struct file_id_offset_pair *)
                    unordered_list_get(gi->offset_index, curr->val);
                */
                struct file_id_offset_pair fid_off = 
                        offset_index_get(gi->offset_idx, curr->val);
                    //gi->offset_idx->index->vals[curr->val];
                struct file_data *fd = file_index_get(gi->file_idx,
                                                      fid_off.file_id);

                file_counts[fid_off.file_id] += 1;

                curr = curr->next;
            }
            uint64_t_ll_free((void **)&R);
        }
    }

    if (line.s != NULL)
        free(line.s);

    mean_interval_size = mean_interval_size/num_intervals;

    struct doubles_uint32_t_tuple *sig = (struct doubles_uint32_t_tuple *)
            calloc(gi->file_idx->index->num,
                   sizeof(struct doubles_uint32_t_tuple));

    uint32_t i;
    for (i = 0; i < gi->file_idx->index->num; ++i) {
        struct file_data *fd = file_index_get(gi->file_idx, i);

        long long n11 = (long long)(file_counts[i]);
        long long n12 = (long long)(MAX(0,num_intervals - file_counts[i]));
        long long n21 = (long long)(MAX(0,fd->num_intervals - file_counts[i]));
        double comp_mean = ((fd->mean_interval_size+mean_interval_size));
        long long n22_full = (long long)
            MAX(n11 + n12 + n21, genome_size/comp_mean);
        long long n22 = MAX(0, n22_full - (n11 + n12 + n21));
        double left, right, two;
        double r = kt_fisher_exact(n11, n12, n21, n22, &left, &right, &two);

        double ratio = (((double)n11/(double)n12) / ((double)n21/(double)n22));

        //fprintf(stderr, "%s\t%f\n", fd->file_name, two);
        sig[i].d1 = right;
        sig[i].d2 = ratio;
        sig[i].u1 = i;
        sig[i].u2 = file_counts[i];
    }

    qsort(sig,
          gi->file_idx->index->num,
          sizeof(struct doubles_uint32_t_tuple), 
          doubles_uint32_t_tuple_cmp);

    for (i = 0; i < gi->file_idx->index->num; ++i) {
        struct file_data *fd = file_index_get(gi->file_idx,  sig[i].u1);
        /*
        printf("%s\t"
               "right:%f\t"
               "%f\n", fd->file_name, sig[i].d1, sig[i].d2);
        */
        printf( "sig:%f\t"
                "size:%u\t"
                "overlap:%u\t"
                "ratio:%f\t"
                "%s\n",
                sig[i].d1,
                fd->num_intervals,
                sig[i].u2,
                sig[i].d2,
                fd->file_name);
    }

    giggle_index_destroy(&gi);
    cache.destroy();
}
