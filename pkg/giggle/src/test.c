#include <stdlib.h>
#include <stdio.h>
#include <err.h>
#include <string.h>

#include "giggle_index.h"
#include "wah.h"
#include "cache.h"
#include "ll.h"

int main(int argc, char **argv)
{
    uint32_t num_chrms = 100;

    if ((argc != 4)) {
        errx(1,
             "usage:\t%s <index dir> <region> <w|i>",
             argv[0]);
    }

    char *index_dir = argv[1];
    char *region_s = argv[2];
    char *i_type = argv[3];

    struct giggle_index *gi;
    gi = giggle_load(index_dir,
                     uint64_t_ll_giggle_set_data_handler);


#if 0
    char *chrm = region_s;
    uint32_t start = 0, end = 0;
    uint32_t i, len = strlen(region_s);
    
    for (i = 0; i < len; ++i) {
        if (region_s[i] == ':') {
            region_s[i] = '\0';
            start = atoi(region_s + i + 1);
        } else if (region_s[i] == '-') {
            region_s[i] = '\0';
            end = atoi(region_s + i + 1);
            break;
        }
    }

    struct giggle_index *gi;
    if (i_type[0] == 'i') {
        gi = giggle_load(index_dir,
                         uint64_t_ll_giggle_set_data_handler);

        struct uint64_t_ll *R =
                (struct uint64_t_ll *)giggle_query_region(gi,
                                                          chrm,
                                                          start,
                                                          end);

        if (R != NULL)
            printf("Hits:%u\n", R->len);
        else
            printf("Hits:0\n");

    } else {
        gi = giggle_load(index_dir,
                         wah_giggle_set_data_handler);

        uint32_t chr_id = giggle_get_chrm_id(gi, chrm);
        //return giggle_search(chr_id, gi->root_ids[chr_id], start, end);
        
        uint32_t domain = chr_id;
        uint32_t root_id = gi->root_ids[chr_id];

        uint32_t leaf_start_id;
        int pos_start_id;

        uint32_t nld_start_id = bpt_find(domain,
                                         root_id,
                                         &leaf_start_id,
                                         &pos_start_id,
                                         start);
        fprintf(stderr,
                "nld_start_id:%u\t"
                "leaf_start_id:%u\t"
                "pos_start_id:%u\n",
                nld_start_id,
                leaf_start_id,
                pos_start_id);

        struct bpt_node *leaf_start = cache.get(domain,
                                                leaf_start_id - 1,
                                                &bpt_node_cache_handler);
        bpt_print_node(leaf_start);

        
        struct wah_bpt_non_leading_data *nld = 
                cache.get(domain,
                          BPT_POINTERS(leaf_start)[0] - 1,
                          &wah_non_leading_cache_handler);

        fprintf(stderr,
                "WAH_LEN:%u\t"
                "wah_get_ints_count:%u\t"
                "\n",
                WAH_LEN(nld->SA),
                wah_get_ints_count(nld->SA));
            
        uint32_t *R = NULL;
        uint32_t R_len = wah_get_ints(nld->SA, &R);

        uint32_t i;
        for (i = 0; i < R_len; ++i) {
            fprintf(stderr, "%u:%u\n", i, R[i]);
        }

        /*
        uint8_t *R = (uint8_t *)giggle_query_region(gi,
                                                    chrm,
                                                    start,
                                                    end);
        if (R != NULL)
            printf("Hits:%u\n", wah_get_ints_count(R));
        else
            printf("Hits:0\n");
        */

    }
#endif
    giggle_index_destroy(&gi);
    cache.destroy();
}
