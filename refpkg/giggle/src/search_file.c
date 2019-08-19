#include <stdlib.h>
#include <stdio.h>
#include <err.h>
#include <string.h>
#include <htslib/kstring.h>

#include "util.h"
#include "giggle_index.h"
#include "wah.h"
#include "cache.h"
#include "file_read.h"
#include "kfunc.h"
#include "ll.h"


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

    struct long_ll **offsets = (struct long_ll **)
            calloc(gi->file_idx->index->num, sizeof(struct long_ll *));

    uint32_t i,j;

    for (i = 0; i < gi->file_idx->index->num; ++i)
        offsets[i] = NULL;


    uint32_t num_intervals = 0;
    while ( in_f->input_file_get_next_interval(in_f, 
                                               &chrm,
                                               &chrm_len,
                                               &start,
                                               &end,
                                               &offset,
                                               &line) >= 0 ) {
        struct uint64_t_ll *R =
                (struct uint64_t_ll *)giggle_query_region(gi,
                                                          chrm,
                                                          start,
                                                          end);
        if (R != NULL) {
            struct uint64_t_ll_node *curr = R->head;

            uint32_t count = 0;
            while (curr != NULL) {
                struct file_id_offset_pair fid_off = 
                    offset_index_get(gi->offset_idx, curr->val);
                    //gi->offset_idx->index->vals[curr->val];

                long_ll_append(&(offsets[fid_off.file_id]),fid_off.offset);

                curr = curr->next;
            }
            uint64_t_ll_free((void **)&R);
            R=NULL;
        } 
        num_intervals += 1;
    }

    if (line.s != NULL)
        free(line.s);

    uint32_t sorted_offsets_num = num_intervals * 2;
    long *sorted_offsets = (long *)malloc(sorted_offsets_num*sizeof(long));
    for (i = 0; i < gi->file_idx->index->num; ++i) {

        struct file_data *fd = file_index_get(gi->file_idx, i);

        printf("#\t%s\t", fd->file_name);

        if (offsets[i] == NULL) {
            printf("0\n");
            continue;
        } else {
            printf("%llu\n", offsets[i]->len);
        }

       
        if (sorted_offsets_num < offsets[i]->len) {
            sorted_offsets_num = offsets[i]->len * 2;
            sorted_offsets = (long *)realloc(sorted_offsets,
                                             sorted_offsets_num*sizeof(long));
        }
    
        j = 0;
        struct long_ll_node *curr = offsets[i]->head;
        while (curr != NULL) {
            sorted_offsets[j++] = curr->val;
            curr = curr->next;
        }

        qsort(sorted_offsets, offsets[i]->len, sizeof(long), long_cmp);

        struct input_file *ipf = input_file_init(fd->file_name);


        char *str;
        for (j = 0; j < offsets[i]->len; ++j) {
            ipf->input_file_seek(ipf, sorted_offsets[j]);
            ipf->input_file_get_next_line(ipf,
                                          &str);
            printf("%s\n", str);
        }

        input_file_destroy(&ipf);
    }

    giggle_index_destroy(&gi);
    cache.destroy();
}
