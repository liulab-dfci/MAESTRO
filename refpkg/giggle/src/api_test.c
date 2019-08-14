#include <stdlib.h>
#include <stdio.h>
#include <err.h>
#include <string.h>
#include <htslib/kstring.h>

#include "giggle_index.h"
#include "wah.h"
#include "cache.h"
#include "file_read.h"
#include "kfunc.h"
#include "ll.h"

struct uint_flt_pair
{
    uint32_t uint;
    float flt;
};

void offset_data_append_uint_flt(uint8_t *dest, kstring_t *line)
{
    int n;
    int *fields = ksplit(line, '\t', &n);

    struct uint_flt_pair a;
    if (n > 7) {
        a.uint = atoi(line->s + fields[4]);
        a.flt = atof(line->s + fields[6]);
    } else {
        a.uint = 0;
        a.flt = 0;
    }

    memcpy(dest, &a, sizeof(struct uint_flt_pair));

    free(fields);
}


int main(int argc, char **argv)
{

    if ((argc != 6)) {
        errx(1,
             "usage:\t%s <input path> <output path> <chrom> <start> <end>\n",
             argv[0]);
    }

    char *input_path = argv[1];
    char *index_path = argv[2];
    char *chrom = argv[3];
    uint32_t start = atoi(argv[4]);
    uint32_t end = atoi(argv[5]);

    offset_data_size =
            sizeof(struct uint_flt_pair);
    offset_data_append_data = offset_data_append_uint_flt;

    uint64_t num_intervals = giggle_bulk_insert(input_path,
                                                index_path,
                                                1);

    struct giggle_index *gi = giggle_load(index_path,
                                          uint64_t_ll_giggle_set_data_handler);

    giggle_data_handler.giggle_collect_intersection =
        giggle_collect_intersection_data_in_block;

    giggle_data_handler.map_intersection_to_offset_list =
            leaf_data_map_intersection_to_offset_list;

    struct giggle_query_result *gqr = giggle_query(gi, chrom, start, end, NULL);

    uint32_t i;
    for(i = 0; i < gqr->num_files; i++) {
        char *result = NULL;
        struct giggle_query_iter *gqi = giggle_get_query_itr(gqr, i);
        while (giggle_query_next(gqi, &result) == 0) {
            printf("line\t%s\n", result);
        }
        giggle_iter_destroy(&gqi);
    }

    giggle_query_result_destroy(&gqr);

    gqr = giggle_query(gi, chrom, start, end, NULL);

    for(i = 0; i < gqr->num_files; i++) {
        struct uint_flt_pair *pair = NULL;
        struct giggle_query_iter *gqi = giggle_get_query_data_itr(gqr, i);
        while (giggle_query_next_data(gqi, (void **)&pair) == 0) {
            printf("data\t%u\t%f\n", pair->uint, pair->flt);
        }
        giggle_iter_destroy(&gqi);
    }

    giggle_index_destroy(&gi);

    cache.destroy();
}
