#define _GNU_SOURCE

#include "util.h"

#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <inttypes.h>
#include <string.h>
#include <err.h>
#include <sysexits.h>

#include "unity.h"
#include "fastlz.h"
#include "bpt.h"
#include "giggle_index.h"
#include "lists.h"
#include "file_read.h"
#include "wah.h"
#include "cache.h"
#include "ll.h"

#define MIN(X, Y) (((X) < (Y)) ? (X) : (Y))

void setUp(void) { }
void tearDown(void) { }

//
////{{{ void test_get_leaf_data(void)
//void test_get_leaf_data(void)
//{
//    ORDER = 10;
//    struct simple_cache *sc = simple_cache_init(1000, 30, NULL);
//    uint64_t_ll_giggle_set_data_handler();
//    struct giggle_index *gi =
//            giggle_init_index(30, "test_get_leaf_data_offset.idx");
//    char *file_name = "../data/1k.unsort.bed.gz";
//    uint32_t ret = giggle_index_file(gi, file_name);
//
//    uint32_t i;
//    for (i = 0; i < gi->num; ++i) {
//        struct bpt_node *curr = cache.get(i,
//                                          gi->root_ids[i] - 1,
//                                          &bpt_node_cache_handler);
//
//        while (BPT_IS_LEAF(curr) != 1) {
//            curr = cache.get(i,
//                             BPT_POINTERS(curr)[0] - 1,
//                             &bpt_node_cache_handler);
//        }
//
//        uint32_t id = BPT_ID(curr);
//        uint32_t len = 0;
//        while (id != 0) {
//            curr = cache.get(i,
//                             id - 1,
//                             &bpt_node_cache_handler);
//
//            struct leaf_data *lf = NULL;
//            uint16_t *starts_ends_offsets = NULL;
//                uint32_t leaf_data_size = 
//                        giggle_get_leaf_data(gi,
//                                             i,
//                                             BPT_ID(curr),
//                                             &lf,
//                                             &starts_ends_offsets);
// 
//            if (BPT_LEADING(curr) != 0) {
//                struct uint64_t_ll_bpt_leading_data *ld =
//                    cache.get(i,
//                              BPT_LEADING(curr) - 1,
//                              &uint64_t_ll_leading_cache_handler);
//                TEST_ASSERT_EQUAL(ld->B->len, lf->num_leading);
//            } else {
//                TEST_ASSERT_EQUAL(0, lf->num_leading);
//            }
//
//            uint32_t s_len = 0;
//            uint32_t e_len = 0;
//            uint32_t j;
//            for (j = 0; j < BPT_NUM_KEYS(curr) ; ++j) {
//                if (BPT_POINTERS(curr)[j] != 0) {
//                    struct uint64_t_ll_bpt_non_leading_data *nld =
//                            cache.get(i,
//                                      BPT_POINTERS(curr)[j] - 1,
//                                      &uint64_t_ll_non_leading_cache_handler);
//                    s_len += (nld->SA == NULL) ? 0 : nld->SA->len;
//                    e_len += (nld->SE == NULL) ? 0 : nld->SE->len;
//
//                    TEST_ASSERT_EQUAL((nld->SA == NULL) ? 0 : nld->SA->len,
//                                      starts_ends_offsets[j*2] -
//                                      ((j==0) ? 0 
//                                              : starts_ends_offsets[j*2-2]));
//                    uint32_t k;
//                    for (k = (j==0)? 0 : starts_ends_offsets[j*2-2];
//                         k < starts_ends_offsets[j*2];
//                         ++k) {
//                        TEST_ASSERT_EQUAL(1,
//                                          uint64_t_ll_contains(nld->SA,
//                                                               lf->starts[k]));
//                    }
//
//                    TEST_ASSERT_EQUAL((nld->SE == NULL) ? 0 : nld->SE->len,
//                                      starts_ends_offsets[j*2+1] -
//                                      ((j==0) ? 0 
//                                              : starts_ends_offsets[j*2-1]));
//
//                    for (k = (j==0)? 0 : starts_ends_offsets[j*2-1];
//                         k < starts_ends_offsets[j*2+1];
//                         ++k) {
//                        TEST_ASSERT_EQUAL(1,
//                                          uint64_t_ll_contains(nld->SE,
//                                                               lf->ends[k]));
//                    }
//                }
//            }
//            TEST_ASSERT_EQUAL(s_len, lf->num_starts);
//            TEST_ASSERT_EQUAL(e_len, lf->num_ends);
//
//            id = BPT_NEXT(curr);
//            leaf_data_free_mem((void **)&lf);
//            free(starts_ends_offsets);
//        }
//    }
//    giggle_index_destroy(&gi);
//    cache.destroy();
//    remove("test_get_leaf_data_offset.idx");
//}
////}}}
//
////{{{ void test_get_leaf_data(void)
//void test_leaf_data_cache_handler(void)
//{
//    ORDER = 10;
//    struct simple_cache *sc = simple_cache_init(1000, 30, NULL);
//    uint64_t_ll_giggle_set_data_handler();
//    struct giggle_index *gi = 
//            giggle_init_index(30,
//                              "test_leaf_data_cache_handler_offset.idx");
//    char *file_name = "../data/1k.unsort.bed.gz";
//    uint32_t ret = giggle_index_file(gi, file_name);
//
//    uint32_t i;
//    for (i = 0; i < gi->num; ++i) {
//        struct bpt_node *curr = cache.get(i,
//                                          gi->root_ids[i] - 1,
//                                          &bpt_node_cache_handler);
//
//        while (BPT_IS_LEAF(curr) != 1) {
//            curr = cache.get(i,
//                             BPT_POINTERS(curr)[0] - 1,
//                             &bpt_node_cache_handler);
//        }
//
//        uint32_t id = BPT_ID(curr);
//        uint32_t len = 0;
//        while (id != 0) {
//            curr = cache.get(i,
//                             id - 1,
//                             &bpt_node_cache_handler);
//
//            struct leaf_data *lf = NULL;
//            uint16_t *starts_ends_offsets = NULL;
//                uint32_t leaf_data_size = 
//                        giggle_get_leaf_data(gi,
//                                             i,
//                                             BPT_ID(curr),
//                                             &lf,
//                                             &starts_ends_offsets);
//
//            /*
//            leaf_data_deserialize(void *serialized,
//                               uint64_t serialized_size,
//                               void **deserialized);
//            leaf_data_serialize(void *deserialized, void **serialized);
//
//
//            */
//            id = BPT_NEXT(curr);
//            leaf_data_free_mem((void **)&lf);
//            free(starts_ends_offsets);
//        }
//    }
//
//    giggle_index_destroy(&gi);
//    cache.destroy();
//    remove("test_leaf_data_cache_handler_offset.idx");
//}
//
////}}}

//{{{void test_leaf_data_ops(void)
void test_leaf_data_ops(void)
{
    ORDER = 5;
    struct giggle_index *gi = giggle_init(
                23,
                "tmp",
                1,
                uint64_t_ll_giggle_set_data_handler);

    giggle_data_handler.write_tree = &giggle_write_tree_leaf_data;

    uint32_t domain = 0;

    //  1    2    3    4    5    6    7    8   9   10   11   12   13   14   15
    //1:|--------------------------------------|
    //2:     |---------|
    //3:          |--------------|
    //4:               |-----------------------------------------------|
    //5:                    |----------------------|
    //6:                              |---------------------------|
    //7:                              |------------|
    //8:                                                     |--------------|
    //
    //Q:                              |------------------|

    //   L   1    2    3    4
    //   --------------------
    // S     1    2    3    4
    // E 
    //
    //   L        5    7    10
    //   ---------------------
    // S 1,2,3,4  5    6,7   
    // E          2    3    1
    //
    //   L     11   12   14   15    16
    //   ------------------------
    // S 4,6,7      8
    // E       5    7    6    4     8
     

    uint32_t root_id = 0;
    char *chrm = "0";
    uint32_t chrm_id = giggle_get_chrm_id(gi, chrm);
    giggle_insert(chrm_id, &root_id, 1, 9, 1);
    giggle_insert(chrm_id, &root_id, 2, 4, 2);
    giggle_insert(chrm_id, &root_id, 3, 6, 3);
    giggle_insert(chrm_id, &root_id, 4, 14, 4);
    giggle_insert(chrm_id, &root_id, 5, 10, 5);
    giggle_insert(chrm_id, &root_id, 7, 13, 6);
    giggle_insert(chrm_id, &root_id, 7, 11, 7);
    giggle_insert(chrm_id, &root_id, 12, 15, 8);

    gi->root_ids[chrm_id] = root_id;

    giggle_store(gi);
    giggle_index_destroy(&gi);
    cache.destroy();

    gi = giggle_load("tmp",
        uint64_t_ll_giggle_set_data_handler);

    struct leaf_data_result *R = NULL;

    uint32_t start = 2, end = 11;

    //{{{ search
    uint32_t leaf_start_id = 0;
    int pos_start_id = 0;

    uint32_t nld_start_id = bpt_find(domain,
                                     gi->root_ids[domain],
                                     &leaf_start_id, 
                                     &pos_start_id,
                                     start);
    struct bpt_node *leaf_start = cache.get(domain,
                                            leaf_start_id - 1,
                                            &bpt_node_cache_handler);
    if ((pos_start_id == 0) && (BPT_KEYS(leaf_start)[0] != start))
        pos_start_id = -1;
    else if ( (pos_start_id >=0) && 
              (pos_start_id < BPT_NUM_KEYS(leaf_start)) &&
              (BPT_KEYS(leaf_start)[pos_start_id] > start))
        pos_start_id -= 1;


    uint32_t leaf_end_id = 0;
    int pos_end_id = 0;

    uint32_t nld_end_id = bpt_find(domain,
                                  gi->root_ids[domain],
                                  &leaf_end_id, 
                                  &pos_end_id,
                                  end);
    struct bpt_node *leaf_end = cache.get(domain,
                                          leaf_end_id - 1,
                                          &bpt_node_cache_handler);

    if ((pos_end_id == 0) && (BPT_KEYS(leaf_end)[0] != end))
        pos_end_id = -1;
    else if ( (pos_end_id >=0) && 
              (pos_end_id < BPT_NUM_KEYS(leaf_end)) &&
              (BPT_KEYS(leaf_end)[pos_end_id] > end))
        pos_end_id -= 1;
    //}}}
    
    R = giggle_collect_intersection_data_in_block(leaf_start_id,
                                                  pos_start_id,
                                                  leaf_end_id,
                                                  pos_end_id,
                                                  domain,
                                                  NULL);
    // 2,11 : 1,2,3,4,5,6,7
    TEST_ASSERT_EQUAL(7, R->len);
    uint32_t A_2_11[7] = {1,2,3,4,5,6,7};

    uint32_t i;
    for (i = 0; i < R->len; ++i) 
        TEST_ASSERT_EQUAL(A_2_11[i], R->data[i]);

    free(R->data);
    free(R);

    start = 6;
    end = 8;

    //{{{ search
    nld_start_id = bpt_find(domain,
                            gi->root_ids[domain],
                            &leaf_start_id, 
                            &pos_start_id,
                            start);
    leaf_start = cache.get(domain,
                           leaf_start_id - 1,
                           &bpt_node_cache_handler);
    if ((pos_start_id == 0) && (BPT_KEYS(leaf_start)[0] != start))
        pos_start_id = -1;
    else if ( (pos_start_id >=0) && 
              (pos_start_id < BPT_NUM_KEYS(leaf_start)) &&
              (BPT_KEYS(leaf_start)[pos_start_id] > start))
        pos_start_id -= 1;


    nld_end_id = bpt_find(domain,
                          gi->root_ids[domain],
                          &leaf_end_id, 
                          &pos_end_id,
                          end);

    leaf_end = cache.get(domain,
                         leaf_end_id - 1,
                         &bpt_node_cache_handler);
    if ((pos_end_id == 0) && (BPT_KEYS(leaf_end)[0] != end))
        pos_end_id = -1;
    else if ( (pos_end_id >=0) && 
              (pos_end_id < BPT_NUM_KEYS(leaf_end)) &&
              (BPT_KEYS(leaf_end)[pos_end_id] > end))
        pos_end_id -= 1;
    //}}}
    
    R = giggle_collect_intersection_data_in_block(leaf_start_id,
                                                  pos_start_id,
                                                  leaf_end_id,
                                                  pos_end_id,
                                                  domain,
                                                  NULL);
    // 6,8 : 1,3,4,5,6,7
    TEST_ASSERT_EQUAL(6, R->len);
    uint32_t A_6_8[6] = {1,3,4,5,6,7};
    for (i = 0; i < R->len; ++i) 
        TEST_ASSERT_EQUAL(A_6_8[i], R->data[i]);

    free(R->data);
    free(R);

    start = 10;
    end = 11;

    //{{{ search
    nld_start_id = bpt_find(domain,
                            gi->root_ids[domain],
                            &leaf_start_id, 
                            &pos_start_id,
                            start);
    leaf_start = cache.get(domain,
                           leaf_start_id - 1,
                           &bpt_node_cache_handler);
    if ((pos_start_id == 0) && (BPT_KEYS(leaf_start)[0] != start))
        pos_start_id = -1;
    else if ( (pos_start_id >=0) && 
              (pos_start_id < BPT_NUM_KEYS(leaf_start)) &&
              (BPT_KEYS(leaf_start)[pos_start_id] > start))
        pos_start_id -= 1;


    nld_end_id = bpt_find(domain,
                          gi->root_ids[domain],
                          &leaf_end_id, 
                          &pos_end_id,
                          end);

    leaf_end = cache.get(domain,
                         leaf_end_id - 1,
                         &bpt_node_cache_handler);
    if ((pos_end_id == 0) && (BPT_KEYS(leaf_end)[0] != end))
        pos_end_id = -1;
    else if ( (pos_end_id >=0) && 
              (pos_end_id < BPT_NUM_KEYS(leaf_end)) &&
              (BPT_KEYS(leaf_end)[pos_end_id] > end))
        pos_end_id -= 1;
    //}}}
 
    R = giggle_collect_intersection_data_in_block(leaf_start_id,
                                                  pos_start_id,
                                                  leaf_end_id,
                                                  pos_end_id,
                                                  domain,
                                                  NULL);


    // 10,12 : 4,5,6,7
    TEST_ASSERT_EQUAL(4, R->len);
    uint32_t A_10_12[4] = {4,5,6,7};
    for (i = 0; i < R->len; ++i) 
        TEST_ASSERT_EQUAL(A_10_12[i], R->data[i]);

    free(R->data);
    free(R);
    giggle_index_destroy(&gi);
    cache.destroy();
    rmrf("tmp");
}
//}}}

