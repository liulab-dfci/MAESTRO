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
#include "bpt.h"
#include "giggle_index.h"
#include "lists.h"
#include "file_read.h"
#include "wah.h"
#include "ll.h"

void valid_giggle_index(struct giggle_index *gi);

void setUp(void) { }
void tearDown(void) { }

//{{{void test_uint64_t_ll_giggle_insert(void)
void test_uint64_t_ll_giggle_insert(void)
{
    ORDER = 4;
    struct simple_cache *sc = simple_cache_init(5, 1, NULL);
    uint32_t domain = 0;
    uint64_t_ll_giggle_set_data_handler();
    //wah_giggle_set_data_handler();

    uint32_t root_id = 0;
    uint32_t r = giggle_insert(domain, &root_id, 1, 3, 0);

    // 1(SA:0),4(SE:0)
    struct bpt_node *root = cache.get(domain,
                                      root_id - 1,
                                      &bpt_node_cache_handler);
    TEST_ASSERT_EQUAL(2, BPT_NUM_KEYS(root));
    TEST_ASSERT_EQUAL(1, BPT_KEYS(root)[0]);
    TEST_ASSERT_EQUAL(4, BPT_KEYS(root)[1]);

    r = giggle_insert(domain, &root_id, 2, 4, 1);

    // 1(SA:0),2(SA:1),4(SE:0),5(SE:1)
    root = cache.get(domain, root_id - 1, &bpt_node_cache_handler);
    TEST_ASSERT_EQUAL(4, BPT_NUM_KEYS(root));
    TEST_ASSERT_EQUAL(1, BPT_KEYS(root)[0]);
    TEST_ASSERT_EQUAL(2, BPT_KEYS(root)[1]);
    TEST_ASSERT_EQUAL(4, BPT_KEYS(root)[2]);
    TEST_ASSERT_EQUAL(5, BPT_KEYS(root)[3]);

    r = giggle_insert(domain, &root_id, 4, 6, 2);
    // 4
    // 1(SA:0),2(SA:1) (0 1)4(SA:2 SE:0),5(SE:1),7(SE:2)
    root = cache.get(domain, root_id - 1, &bpt_node_cache_handler);
    TEST_ASSERT_EQUAL(1, BPT_NUM_KEYS(root));
    TEST_ASSERT_EQUAL(4, BPT_KEYS(root)[0]);

    struct bpt_node *first_leaf = cache.get(domain,
                                            BPT_POINTERS(root)[0] - 1,
                                            &bpt_node_cache_handler);
    TEST_ASSERT_EQUAL(2, BPT_NUM_KEYS(first_leaf));
    TEST_ASSERT_EQUAL(1, BPT_KEYS(first_leaf)[0]);
    TEST_ASSERT_EQUAL(2, BPT_KEYS(first_leaf)[1]);

    struct bpt_node *second_leaf = cache.get(domain,
                                             BPT_POINTERS(root)[1] - 1,
                                             &bpt_node_cache_handler);
    TEST_ASSERT_EQUAL(BPT_ID(second_leaf), BPT_NEXT(first_leaf));

    TEST_ASSERT_EQUAL(3, BPT_NUM_KEYS(second_leaf));
    TEST_ASSERT_EQUAL(4, BPT_KEYS(second_leaf)[0]);
    TEST_ASSERT_EQUAL(5, BPT_KEYS(second_leaf)[1]);
    TEST_ASSERT_EQUAL(7, BPT_KEYS(second_leaf)[2]);
    TEST_ASSERT_EQUAL(0, BPT_LEADING(first_leaf));
    TEST_ASSERT_TRUE(BPT_LEADING(second_leaf) != 0);

    struct uint64_t_ll_bpt_leading_data *ld =
            (struct uint64_t_ll_bpt_leading_data *)
            cache.get(domain,
                      BPT_LEADING(second_leaf) - 1,
                      &giggle_data_handler.leading_cache_handler);
    TEST_ASSERT_EQUAL(2, ld->B->len);
    TEST_ASSERT_EQUAL(0, ld->B->head->val);
    TEST_ASSERT_EQUAL(1, ld->B->head->next->val);

    r = giggle_insert(domain, &root_id, 4, 5, 3);
    // 4
    // 1(SA:0),2(SA:1) (0 1)4(SA:2,3 SE:0),5(SE:1),6(SE:3),7(SE:2)
    first_leaf = cache.get(domain,
                           BPT_POINTERS(root)[0] - 1,
                           &bpt_node_cache_handler);
    TEST_ASSERT_EQUAL(0, BPT_LEADING(first_leaf));
    TEST_ASSERT_EQUAL(2, BPT_NUM_KEYS(first_leaf));
    struct uint64_t_ll_bpt_non_leading_data *nld = 
        cache.get(domain,
                  BPT_POINTERS(first_leaf)[0] - 1,
                  &giggle_data_handler.non_leading_cache_handler);
    TEST_ASSERT_EQUAL(NULL, nld->SE);
    TEST_ASSERT_EQUAL(1, nld->SA->len);
    TEST_ASSERT_EQUAL(0, nld->SA->head->val);
    nld = cache.get(domain,
                    BPT_POINTERS(first_leaf)[1] - 1,
                    &giggle_data_handler.non_leading_cache_handler);
    TEST_ASSERT_EQUAL(NULL, nld->SE);
    TEST_ASSERT_EQUAL(1, nld->SA->len);
    TEST_ASSERT_EQUAL(1, nld->SA->head->val);

    struct bpt_node *next_leaf = 
            cache.get(domain,
                      BPT_NEXT(first_leaf) - 1,
                      &bpt_node_cache_handler);
    ld = cache.get(domain,
                   BPT_LEADING(next_leaf) - 1,
                   &giggle_data_handler.leading_cache_handler);

    TEST_ASSERT_EQUAL(2, ld->B->len);
    TEST_ASSERT_EQUAL(0, ld->B->head->val);
    TEST_ASSERT_EQUAL(1, ld->B->head->next->val);
    TEST_ASSERT_EQUAL(4, BPT_NUM_KEYS(next_leaf));
    TEST_ASSERT_EQUAL(4, BPT_KEYS(next_leaf)[0]);
    TEST_ASSERT_EQUAL(5, BPT_KEYS(next_leaf)[1]);
    TEST_ASSERT_EQUAL(6, BPT_KEYS(next_leaf)[2]);
    TEST_ASSERT_EQUAL(7, BPT_KEYS(next_leaf)[3]);
    nld = cache.get(domain,
                    BPT_POINTERS(next_leaf)[0] - 1,
                    &giggle_data_handler.non_leading_cache_handler);
    TEST_ASSERT_EQUAL(2, nld->SA->len);
    TEST_ASSERT_EQUAL(2, nld->SA->head->val);
    TEST_ASSERT_EQUAL(3, nld->SA->head->next->val);
    TEST_ASSERT_EQUAL(1, nld->SE->len);
    TEST_ASSERT_EQUAL(0, nld->SE->head->val);

    nld = cache.get(domain,
                    BPT_POINTERS(next_leaf)[1] - 1,
                    &giggle_data_handler.non_leading_cache_handler);
    TEST_ASSERT_EQUAL(NULL, nld->SA);
    TEST_ASSERT_EQUAL(1, nld->SE->len);
    TEST_ASSERT_EQUAL(1, nld->SE->head->val);

    nld = cache.get(domain,
                    BPT_POINTERS(next_leaf)[2] - 1,
                    &giggle_data_handler.non_leading_cache_handler);
    TEST_ASSERT_EQUAL(NULL, nld->SA);
    TEST_ASSERT_EQUAL(1, nld->SE->len);
    TEST_ASSERT_EQUAL(3, nld->SE->head->val);

    nld = cache.get(domain,
                    BPT_POINTERS(next_leaf)[3] - 1,
                    &giggle_data_handler.non_leading_cache_handler);
    TEST_ASSERT_EQUAL(NULL, nld->SA);
    TEST_ASSERT_EQUAL(1, nld->SE->len);
    TEST_ASSERT_EQUAL(2, nld->SE->head->val);

    r = giggle_insert(domain, &root_id, 1, 6, 4);
    // 4
    // 1(SA:0,4),2(SA:1) (0 1 4)4(SA:2,3 SE:0),5(SE:1),6(SE:3),7(SE:2,4)
    root = cache.get(domain,
                     root_id - 1,
                     &bpt_node_cache_handler);
    first_leaf = cache.get(domain,
                           BPT_POINTERS(root)[0] - 1,
                           &bpt_node_cache_handler);
    TEST_ASSERT_EQUAL(0, BPT_LEADING(first_leaf));
    TEST_ASSERT_EQUAL(2, BPT_NUM_KEYS(first_leaf));

    nld = cache.get(domain,
                    BPT_POINTERS(first_leaf)[0] - 1,
                    &giggle_data_handler.non_leading_cache_handler);
    TEST_ASSERT_EQUAL(NULL, nld->SE);
    TEST_ASSERT_EQUAL(2, nld->SA->len);
    TEST_ASSERT_EQUAL(0, nld->SA->head->val);
    TEST_ASSERT_EQUAL(4, nld->SA->head->next->val);

    nld = cache.get(domain,
                    BPT_POINTERS(first_leaf)[1] - 1,
                    &giggle_data_handler.non_leading_cache_handler);
    TEST_ASSERT_EQUAL(NULL, nld->SE);
    TEST_ASSERT_EQUAL(1, nld->SA->len);
    TEST_ASSERT_EQUAL(1, nld->SA->head->val);

    next_leaf = cache.get(domain,
                          BPT_NEXT(first_leaf) - 1,
                          &bpt_node_cache_handler);
    ld = cache.get(domain,
                   BPT_LEADING(next_leaf) - 1,
                   &giggle_data_handler.leading_cache_handler);
    TEST_ASSERT_EQUAL(3, ld->B->len);
    TEST_ASSERT_EQUAL(0, ld->B->head->val);
    TEST_ASSERT_EQUAL(1, ld->B->head->next->val);
    TEST_ASSERT_EQUAL(4, ld->B->head->next->next->val);

    TEST_ASSERT_EQUAL(4, BPT_NUM_KEYS(next_leaf));
    TEST_ASSERT_EQUAL(4, BPT_KEYS(next_leaf)[0]);
    TEST_ASSERT_EQUAL(5, BPT_KEYS(next_leaf)[1]);
    TEST_ASSERT_EQUAL(6, BPT_KEYS(next_leaf)[2]);
    TEST_ASSERT_EQUAL(7, BPT_KEYS(next_leaf)[3]);

    nld = cache.get(domain,
                    BPT_POINTERS(next_leaf)[0] - 1,
                    &giggle_data_handler.non_leading_cache_handler);
    TEST_ASSERT_EQUAL(2, nld->SA->len);
    TEST_ASSERT_EQUAL(2, nld->SA->head->val);
    TEST_ASSERT_EQUAL(3, nld->SA->head->next->val);
    TEST_ASSERT_EQUAL(1, nld->SE->len);
    TEST_ASSERT_EQUAL(0, nld->SE->head->val);

    nld = cache.get(domain,
                    BPT_POINTERS(next_leaf)[1] - 1,
                    &giggle_data_handler.non_leading_cache_handler);
    TEST_ASSERT_EQUAL(NULL, nld->SA);
    TEST_ASSERT_EQUAL(1, nld->SE->len);
    TEST_ASSERT_EQUAL(1, nld->SE->head->val);

    nld = cache.get(domain,
                    BPT_POINTERS(next_leaf)[2] - 1,
                    &giggle_data_handler.non_leading_cache_handler);
    TEST_ASSERT_EQUAL(NULL, nld->SA);
    TEST_ASSERT_EQUAL(1, nld->SE->len);
    TEST_ASSERT_EQUAL(3, nld->SE->head->val);

    nld = cache.get(domain,
                    BPT_POINTERS(next_leaf)[3] - 1,
                    &giggle_data_handler.non_leading_cache_handler);
    TEST_ASSERT_EQUAL(NULL, nld->SA);
    TEST_ASSERT_EQUAL(2, nld->SE->len);
    TEST_ASSERT_EQUAL(2, nld->SE->head->val);
    TEST_ASSERT_EQUAL(4, nld->SE->head->next->val);

    r = giggle_insert(domain, &root_id, 1, 7, 5);
    // 4,6
    //   (NULL)    1(SA:0,4,5),2(SA:1)
    //   (0 1 4 5) 4(SA:2,3 SE:0),5(SE:1)
    //   (2 3 4 5) 6(SE:3),7(SE:2,4),8(SE:5)
    root = cache.get(domain,
                     root_id - 1,
                     &bpt_node_cache_handler);
    TEST_ASSERT_EQUAL(2, BPT_NUM_KEYS(root));
    TEST_ASSERT_EQUAL(4, BPT_KEYS(root)[0]);
    TEST_ASSERT_EQUAL(6, BPT_KEYS(root)[1]);
    first_leaf = cache.get(domain,
                           BPT_POINTERS(root)[0] - 1,
                           &bpt_node_cache_handler);
    TEST_ASSERT_EQUAL(0, BPT_LEADING(first_leaf));

    next_leaf = cache.get(domain,
                          BPT_NEXT(first_leaf) - 1,
                          &bpt_node_cache_handler);
    ld = cache.get(domain,
                   BPT_LEADING(next_leaf) - 1,
                   &giggle_data_handler.leading_cache_handler);
    TEST_ASSERT_EQUAL(4, ld->B->len);
    TEST_ASSERT_EQUAL(0, ld->B->head->val);
    TEST_ASSERT_EQUAL(1, ld->B->head->next->val);
    TEST_ASSERT_EQUAL(4, ld->B->head->next->next->val);
    TEST_ASSERT_EQUAL(5, ld->B->head->next->next->next->val);

    next_leaf = cache.get(domain,
                          BPT_NEXT(next_leaf) - 1,
                          &bpt_node_cache_handler);
    ld = cache.get(domain,
                   BPT_LEADING(next_leaf) - 1,
                   &giggle_data_handler.leading_cache_handler);
    TEST_ASSERT_EQUAL(4, ld->B->len);
    TEST_ASSERT_EQUAL(4, ld->B->head->val);
    TEST_ASSERT_EQUAL(2, ld->B->head->next->val);
    TEST_ASSERT_EQUAL(3, ld->B->head->next->next->val);
    TEST_ASSERT_EQUAL(5, ld->B->head->next->next->next->val);

    cache.destroy();
}
//}}}

//{{{void test_uint64_t_ll_giggle_insert(void)
void test_wah_giggle_insert(void)
{
    uint32_t *R = NULL, R_size;
    struct simple_cache *sc = simple_cache_init(5, 1, NULL);
    uint32_t domain = 0;
    wah_giggle_set_data_handler();

    uint32_t root_id = 0;
    uint32_t r = giggle_insert(domain, &root_id, 1, 3, 0);

    // 1(SA:0),4(SE:0)
    struct bpt_node *root = cache.get(domain,
                                      root_id - 1,
                                      &bpt_node_cache_handler);
    TEST_ASSERT_EQUAL(2, BPT_NUM_KEYS(root));
    TEST_ASSERT_EQUAL(1, BPT_KEYS(root)[0]);
    TEST_ASSERT_EQUAL(4, BPT_KEYS(root)[1]);

    r = giggle_insert(domain, &root_id, 2, 4, 1);

    // 1(SA:0),2(SA:1),4(SE:0),5(SE:1)
    root = cache.get(domain, root_id - 1, &bpt_node_cache_handler);
    TEST_ASSERT_EQUAL(4, BPT_NUM_KEYS(root));
    TEST_ASSERT_EQUAL(1, BPT_KEYS(root)[0]);
    TEST_ASSERT_EQUAL(2, BPT_KEYS(root)[1]);
    TEST_ASSERT_EQUAL(4, BPT_KEYS(root)[2]);
    TEST_ASSERT_EQUAL(5, BPT_KEYS(root)[3]);

    struct wah_bpt_non_leading_data *nld = 
        cache.get(domain,
                  BPT_POINTERS(root)[0] - 1,
                  &giggle_data_handler.non_leading_cache_handler);
    TEST_ASSERT_EQUAL(NULL, nld->SE);
    R_size = wah_get_ints(nld->SA, &R);
    TEST_ASSERT_EQUAL(1, R_size);
    TEST_ASSERT_EQUAL(0, R[0] - 1);
    free(R);
    R = NULL;

    nld = cache.get(domain,
                    BPT_POINTERS(root)[1] - 1,
                    &giggle_data_handler.non_leading_cache_handler);
    TEST_ASSERT_EQUAL(NULL, nld->SE);
    R_size = wah_get_ints(nld->SA, &R);
    TEST_ASSERT_EQUAL(1, R_size);
    TEST_ASSERT_EQUAL(1, R[0] - 1);
    free(R);
    R = NULL;

    nld = cache.get(domain,
                    BPT_POINTERS(root)[2] - 1,
                    &giggle_data_handler.non_leading_cache_handler);
    TEST_ASSERT_EQUAL(NULL, nld->SA);
    R_size = wah_get_ints(nld->SE, &R);
    TEST_ASSERT_EQUAL(1, R_size);
    TEST_ASSERT_EQUAL(0, R[0] - 1);
    free(R);
    R = NULL;

    nld = cache.get(domain,
                    BPT_POINTERS(root)[3] - 1,
                    &giggle_data_handler.non_leading_cache_handler);
    TEST_ASSERT_EQUAL(NULL, nld->SA);
    R_size = wah_get_ints(nld->SE, &R);
    TEST_ASSERT_EQUAL(1, R_size);
    TEST_ASSERT_EQUAL(1, R[0] - 1);
    free(R);
    R = NULL;

    r = giggle_insert(domain, &root_id, 4, 6, 2);
    // 4
    // 1(SA:0),2(SA:1) (0 1)4(SA:2 SE:0),5(SE:1),7(SE:2)
    root = cache.get(domain, root_id - 1, &bpt_node_cache_handler);
    TEST_ASSERT_EQUAL(1, BPT_NUM_KEYS(root));
    TEST_ASSERT_EQUAL(4, BPT_KEYS(root)[0]);

    struct bpt_node *first_leaf = cache.get(domain,
                                            BPT_POINTERS(root)[0] - 1,
                                            &bpt_node_cache_handler);
    TEST_ASSERT_EQUAL(2, BPT_NUM_KEYS(first_leaf));
    TEST_ASSERT_EQUAL(1, BPT_KEYS(first_leaf)[0]);
    TEST_ASSERT_EQUAL(2, BPT_KEYS(first_leaf)[1]);

    struct bpt_node *second_leaf = cache.get(domain,
                                             BPT_POINTERS(root)[1] - 1,
                                             &bpt_node_cache_handler);
    TEST_ASSERT_EQUAL(BPT_ID(second_leaf), BPT_NEXT(first_leaf));

    TEST_ASSERT_EQUAL(3, BPT_NUM_KEYS(second_leaf));
    TEST_ASSERT_EQUAL(4, BPT_KEYS(second_leaf)[0]);
    TEST_ASSERT_EQUAL(5, BPT_KEYS(second_leaf)[1]);
    TEST_ASSERT_EQUAL(7, BPT_KEYS(second_leaf)[2]);
    TEST_ASSERT_EQUAL(0, BPT_LEADING(first_leaf));
    TEST_ASSERT_TRUE(BPT_LEADING(second_leaf) != 0);

    struct wah_bpt_leading_data *ld =
            (struct wah_bpt_leading_data *)
            cache.get(domain,
                      BPT_LEADING(second_leaf) - 1,
                      &giggle_data_handler.leading_cache_handler);
    R_size = wah_get_ints(ld->B, &R);
    TEST_ASSERT_EQUAL(2, R_size);
    TEST_ASSERT_EQUAL(0, R[0] - 1);
    TEST_ASSERT_EQUAL(1, R[1] - 1);
    free(R);
    R = NULL;

    r = giggle_insert(domain, &root_id, 4, 5, 3);
    // 4
    // 1(SA:0),2(SA:1) (0 1)4(SA:2,3 SE:0),5(SE:1),6(SE:3),7(SE:2)
    first_leaf = cache.get(domain,
                           BPT_POINTERS(root)[0] - 1,
                           &bpt_node_cache_handler);
    TEST_ASSERT_EQUAL(0, BPT_LEADING(first_leaf));
    TEST_ASSERT_EQUAL(2, BPT_NUM_KEYS(first_leaf));

    nld = cache.get(domain,
                    BPT_POINTERS(first_leaf)[0] - 1,
                    &giggle_data_handler.non_leading_cache_handler);
    TEST_ASSERT_EQUAL(NULL, nld->SE);
    R_size = wah_get_ints(nld->SA, &R);
    TEST_ASSERT_EQUAL(1, R_size);
    TEST_ASSERT_EQUAL(0, R[0] - 1);
    free(R);
    R = NULL;

    nld = cache.get(domain,
                    BPT_POINTERS(first_leaf)[1] - 1,
                    &giggle_data_handler.non_leading_cache_handler);
    TEST_ASSERT_EQUAL(NULL, nld->SE);
    R_size = wah_get_ints(nld->SA, &R);
    TEST_ASSERT_EQUAL(1, R_size);
    TEST_ASSERT_EQUAL(1, R[0] - 1);
    free(R);
    R = NULL;

    struct bpt_node *next_leaf = 
            cache.get(domain,
                      BPT_NEXT(first_leaf) - 1,
                      &bpt_node_cache_handler);
    ld = cache.get(domain,
                   BPT_LEADING(next_leaf) - 1,
                   &giggle_data_handler.leading_cache_handler);
    R_size = wah_get_ints(ld->B, &R);
    TEST_ASSERT_EQUAL(2, R_size);
    TEST_ASSERT_EQUAL(0, R[0] - 1);
    TEST_ASSERT_EQUAL(1, R[1] - 1);
    free(R);
    R=NULL;
    TEST_ASSERT_EQUAL(4, BPT_NUM_KEYS(next_leaf));
    TEST_ASSERT_EQUAL(4, BPT_KEYS(next_leaf)[0]);
    TEST_ASSERT_EQUAL(5, BPT_KEYS(next_leaf)[1]);
    TEST_ASSERT_EQUAL(6, BPT_KEYS(next_leaf)[2]);
    TEST_ASSERT_EQUAL(7, BPT_KEYS(next_leaf)[3]);

    nld = cache.get(domain,
                    BPT_POINTERS(next_leaf)[0] - 1,
                    &giggle_data_handler.non_leading_cache_handler);
    R_size = wah_get_ints(nld->SA, &R);
    TEST_ASSERT_EQUAL(2, R_size);
    TEST_ASSERT_EQUAL(2, R[0] - 1);
    TEST_ASSERT_EQUAL(3, R[1] - 1);
    free(R);
    R = NULL;
    R_size = wah_get_ints(nld->SE, &R);
    TEST_ASSERT_EQUAL(1, R_size);
    TEST_ASSERT_EQUAL(0, R[0] - 1);
    free(R);
    R = NULL;

    nld = cache.get(domain,
                    BPT_POINTERS(next_leaf)[1] - 1,
                    &giggle_data_handler.non_leading_cache_handler);
    TEST_ASSERT_EQUAL(NULL, nld->SA);
    R_size = wah_get_ints(nld->SE, &R);
    TEST_ASSERT_EQUAL(1, R_size);
    TEST_ASSERT_EQUAL(1, R[0] - 1);
    free(R);
    R = NULL;

    nld = cache.get(domain,
                    BPT_POINTERS(next_leaf)[2] - 1,
                    &giggle_data_handler.non_leading_cache_handler);
    TEST_ASSERT_EQUAL(NULL, nld->SA);
    R_size = wah_get_ints(nld->SE, &R);
    TEST_ASSERT_EQUAL(1, R_size);
    TEST_ASSERT_EQUAL(3, R[0] - 1);
    free(R);
    R = NULL;

    nld = cache.get(domain,
                    BPT_POINTERS(next_leaf)[3] - 1,
                    &giggle_data_handler.non_leading_cache_handler);
    TEST_ASSERT_EQUAL(NULL, nld->SA);
    R_size = wah_get_ints(nld->SE, &R);
    TEST_ASSERT_EQUAL(1, R_size);
    TEST_ASSERT_EQUAL(2, R[0] - 1);
    free(R);
    R = NULL;

    r = giggle_insert(domain, &root_id, 1, 6, 4);
    // 4
    // 1(SA:0,4),2(SA:1) (0 1 4)4(SA:2,3 SE:0),5(SE:1),6(SE:3),7(SE:2,4)
    root = cache.get(domain,
                     root_id - 1,
                     &bpt_node_cache_handler);
    first_leaf = cache.get(domain,
                           BPT_POINTERS(root)[0] - 1,
                           &bpt_node_cache_handler);
    TEST_ASSERT_EQUAL(0, BPT_LEADING(first_leaf));
    TEST_ASSERT_EQUAL(2, BPT_NUM_KEYS(first_leaf));

    nld = cache.get(domain,
                    BPT_POINTERS(first_leaf)[0] - 1,
                    &giggle_data_handler.non_leading_cache_handler);
    TEST_ASSERT_EQUAL(NULL, nld->SE);
    R_size = wah_get_ints(nld->SA, &R);
    TEST_ASSERT_EQUAL(2, R_size);
    TEST_ASSERT_EQUAL(0, R[0] - 1);
    TEST_ASSERT_EQUAL(4, R[1] - 1);
    free(R);
    R = NULL;

    nld = cache.get(domain,
                    BPT_POINTERS(first_leaf)[1] - 1,
                    &giggle_data_handler.non_leading_cache_handler);
    TEST_ASSERT_EQUAL(NULL, nld->SE);
    R_size = wah_get_ints(nld->SA, &R);
    TEST_ASSERT_EQUAL(1, R_size);
    TEST_ASSERT_EQUAL(1, R[0] - 1);
    free(R);
    R = NULL;

    next_leaf = cache.get(domain,
                          BPT_NEXT(first_leaf) - 1,
                          &bpt_node_cache_handler);
    ld = cache.get(domain,
                   BPT_LEADING(next_leaf) - 1,
                   &giggle_data_handler.leading_cache_handler);
    R_size = wah_get_ints(ld->B, &R);
    TEST_ASSERT_EQUAL(3, R_size);
    TEST_ASSERT_EQUAL(0, R[0] - 1);
    TEST_ASSERT_EQUAL(1, R[1] - 1);
    TEST_ASSERT_EQUAL(4, R[2] - 1);
    free(R);
    R = NULL;

    TEST_ASSERT_EQUAL(4, BPT_NUM_KEYS(next_leaf));
    TEST_ASSERT_EQUAL(4, BPT_KEYS(next_leaf)[0]);
    TEST_ASSERT_EQUAL(5, BPT_KEYS(next_leaf)[1]);
    TEST_ASSERT_EQUAL(6, BPT_KEYS(next_leaf)[2]);
    TEST_ASSERT_EQUAL(7, BPT_KEYS(next_leaf)[3]);

    nld = cache.get(domain,
                    BPT_POINTERS(next_leaf)[0] - 1,
                    &giggle_data_handler.non_leading_cache_handler);
    R_size = wah_get_ints(nld->SA, &R);
    TEST_ASSERT_EQUAL(2, R_size);
    TEST_ASSERT_EQUAL(2, R[0] - 1);
    TEST_ASSERT_EQUAL(3, R[1] - 1);
    free(R);
    R = NULL;
    R_size = wah_get_ints(nld->SE, &R);
    TEST_ASSERT_EQUAL(1, R_size);
    TEST_ASSERT_EQUAL(0, R[0] - 1);
    free(R);
    R = NULL;

    nld = cache.get(domain,
                    BPT_POINTERS(next_leaf)[1] - 1,
                    &giggle_data_handler.non_leading_cache_handler);
    TEST_ASSERT_EQUAL(NULL, nld->SA);
    R_size = wah_get_ints(nld->SE, &R);
    TEST_ASSERT_EQUAL(1, R_size);
    TEST_ASSERT_EQUAL(1, R[0] - 1);
    free(R);
    R = NULL;

    nld = cache.get(domain,
                    BPT_POINTERS(next_leaf)[2] - 1,
                    &giggle_data_handler.non_leading_cache_handler);
    TEST_ASSERT_EQUAL(NULL, nld->SA);
    R_size = wah_get_ints(nld->SE, &R);
    TEST_ASSERT_EQUAL(1, R_size);
    TEST_ASSERT_EQUAL(3, R[0] - 1);
    free(R);
    R = NULL;

    nld = cache.get(domain,
                    BPT_POINTERS(next_leaf)[3] - 1,
                    &giggle_data_handler.non_leading_cache_handler);
    TEST_ASSERT_EQUAL(NULL, nld->SA);
    R_size = wah_get_ints(nld->SE, &R);
    TEST_ASSERT_EQUAL(2, R_size);
    TEST_ASSERT_EQUAL(2, R[0] - 1);
    TEST_ASSERT_EQUAL(4, R[1] - 1);
    free(R);
    R = NULL;

    r = giggle_insert(domain, &root_id, 1, 7, 5);
    // 4,6
    //   (NULL)    1(SA:0,4,5),2(SA:1)
    //   (0 1 4 5) 4(SA:2,3 SE:0),5(SE:1)
    //   (2 3 4 5) 6(SE:3),7(SE:2,4),8(SE:5)
    root = cache.get(domain,
                     root_id - 1,
                     &bpt_node_cache_handler);
    TEST_ASSERT_EQUAL(2, BPT_NUM_KEYS(root));
    TEST_ASSERT_EQUAL(4, BPT_KEYS(root)[0]);
    TEST_ASSERT_EQUAL(6, BPT_KEYS(root)[1]);
    first_leaf = cache.get(domain,
                           BPT_POINTERS(root)[0] - 1,
                           &bpt_node_cache_handler);
    TEST_ASSERT_EQUAL(0, BPT_LEADING(first_leaf));

    next_leaf = cache.get(domain,
                          BPT_NEXT(first_leaf) - 1,
                          &bpt_node_cache_handler);
    ld = cache.get(domain,
                   BPT_LEADING(next_leaf) - 1,
                   &giggle_data_handler.leading_cache_handler);
    R_size = wah_get_ints(ld->B, &R);
    TEST_ASSERT_EQUAL(4, R_size);
    TEST_ASSERT_EQUAL(0, R[0] - 1);
    TEST_ASSERT_EQUAL(1, R[1] - 1);
    TEST_ASSERT_EQUAL(4, R[2] - 1);
    TEST_ASSERT_EQUAL(5, R[3] - 1);
    free(R);
    R = NULL;

    next_leaf = cache.get(domain,
                          BPT_NEXT(next_leaf) - 1,
                          &bpt_node_cache_handler);
    ld = cache.get(domain,
                   BPT_LEADING(next_leaf) - 1,
                   &giggle_data_handler.leading_cache_handler);
    R_size = wah_get_ints(ld->B, &R);
    TEST_ASSERT_EQUAL(4, R_size);
    TEST_ASSERT_EQUAL(2, R[0] - 1);
    TEST_ASSERT_EQUAL(3, R[1] - 1);
    TEST_ASSERT_EQUAL(4, R[2] - 1);
    TEST_ASSERT_EQUAL(5, R[3] - 1);
    free(R);
    R = NULL;

    cache.destroy();
}
//}}}

//{{{ void test_uint64_t_ll_giggle_search(void)
void test_uint64_t_ll_giggle_search(void)
{
    ORDER = 7;
    struct simple_cache *sc = simple_cache_init(5, 1, NULL);
    uint32_t domain = 0;
    uint64_t_ll_giggle_set_data_handler();

    /*
     *  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19 20
     *  |0-------------------------|  |1----|     |9------| 
     *     |2----|  |3-------|
     *        |4-------|  |5----|     |7-------------------|
     *                    |6------------|
     *  |8----------------------------------------------|
     *
     */
    uint32_t root_id = 0;

    uint32_t r = giggle_insert(domain, &root_id, 1, 10, 0);
    r = giggle_insert(domain, &root_id, 11, 13, 1);
    r = giggle_insert(domain, &root_id, 2, 4, 2);
    r = giggle_insert(domain, &root_id, 5, 8, 3);
    r = giggle_insert(domain, &root_id, 3, 6, 4);
    r = giggle_insert(domain, &root_id, 7, 9, 5);
    r = giggle_insert(domain, &root_id, 7, 12, 6);

    r = giggle_insert(domain, &root_id, 11, 18, 7);
    r = giggle_insert(domain, &root_id, 1, 17, 8);
    r = giggle_insert(domain, &root_id, 15, 18, 9);

    /*
    7 13
      (NULL)    1(SA:0,8) 2(SA:2) 3(SA:4) 5(SA:3 SE:2) 
      (0 3 4 8) 7(SA:5,6 SE:4) 9(SE: 3) 10(SE:5) 11(SA:1,7 SE:0) 
      (8 6 1 7) 13(SE:6) 14(SE:1) 15(SA:9) 18(SE:8) 19(SE:7,9)
    */
#if 0
    //{{{
    struct bpt_node *root = cache.get(cache.cache, root_id);
    TEST_ASSERT_EQUAL(2, root->num_keys);
    TEST_ASSERT_EQUAL(7, root->keys[0]);
    TEST_ASSERT_EQUAL(13, root->keys[1]);

    struct bpt_node *first_leaf = (struct bpt_node *) root->pointers[0];
    TEST_ASSERT_EQUAL(4, first_leaf->num_keys);
    TEST_ASSERT_EQUAL(1, first_leaf->keys[0]);
    TEST_ASSERT_EQUAL(2, first_leaf->keys[1]);
    TEST_ASSERT_EQUAL(3, first_leaf->keys[2]);
    TEST_ASSERT_EQUAL(5, first_leaf->keys[3]);

    struct uint64_t_ll_bpt_leading_data *ld =
            (struct uint64_t_ll_bpt_leading_data *)
            first_leaf->leading;
    TEST_ASSERT_EQUAL(NULL, ld);

    TEST_ASSERT_EQUAL(4, first_leaf->next->num_keys);
    TEST_ASSERT_EQUAL(7, first_leaf->next->keys[0]);
    TEST_ASSERT_EQUAL(9, first_leaf->next->keys[1]);
    TEST_ASSERT_EQUAL(10, first_leaf->next->keys[2]);
    TEST_ASSERT_EQUAL(11, first_leaf->next->keys[3]);

    ld = (struct uint64_t_ll_bpt_leading_data *) first_leaf->next->leading;
    TEST_ASSERT_EQUAL(4, ld->B->len);
    TEST_ASSERT_EQUAL(1, uint64_t_ll_contains(ld->B, 0));
    TEST_ASSERT_EQUAL(1, uint64_t_ll_contains(ld->B, 3));
    TEST_ASSERT_EQUAL(1, uint64_t_ll_contains(ld->B, 4));
    TEST_ASSERT_EQUAL(1, uint64_t_ll_contains(ld->B, 8));

    TEST_ASSERT_EQUAL(5, first_leaf->next->next->num_keys);
    TEST_ASSERT_EQUAL(13, first_leaf->next->next->keys[0]);
    TEST_ASSERT_EQUAL(14, first_leaf->next->next->keys[1]);
    TEST_ASSERT_EQUAL(15, first_leaf->next->next->keys[2]);
    TEST_ASSERT_EQUAL(18, first_leaf->next->next->keys[3]);
    TEST_ASSERT_EQUAL(19, first_leaf->next->next->keys[4]);

    ld = (struct uint64_t_ll_bpt_leading_data *)
            first_leaf->next->next->leading;
    TEST_ASSERT_EQUAL(4, ld->B->len);
    TEST_ASSERT_EQUAL(1, uint64_t_ll_contains(ld->B, 8));
    TEST_ASSERT_EQUAL(1, uint64_t_ll_contains(ld->B, 6));
    TEST_ASSERT_EQUAL(1, uint64_t_ll_contains(ld->B, 1));
    TEST_ASSERT_EQUAL(1, uint64_t_ll_contains(ld->B, 7));
    //}}}
#endif

    /*
     *  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19 20
     *  |0-------------------------|  |1----|     |9------| 
     *     |2----|  |3-------|
     *        |4-------|  |5----|     |7-------------------|
     *                    |6------------|
     *  |8----------------------------------------------|
     *
     */

    ///struct uint64_t_ll *R = (struct uint64_t_ll *)giggle_search(root, 2, 5);
    //uint32_t R_id = 
    struct uint64_t_ll *R = giggle_search(domain, root_id, 2, 5);
    
    TEST_ASSERT_TRUE(R != NULL);
    TEST_ASSERT_EQUAL(5, R->len);
    TEST_ASSERT_EQUAL(1, uint64_t_ll_contains(R, 0));
    TEST_ASSERT_EQUAL(1, uint64_t_ll_contains(R, 2));
    TEST_ASSERT_EQUAL(1, uint64_t_ll_contains(R, 4));
    TEST_ASSERT_EQUAL(1, uint64_t_ll_contains(R, 3));
    TEST_ASSERT_EQUAL(1, uint64_t_ll_contains(R, 8));

    uint64_t_ll_free((void **)&R);
    R = NULL;

    R = giggle_search(domain, root_id, 5, 15);
    TEST_ASSERT_EQUAL(9, R->len);
    TEST_ASSERT_EQUAL(1, uint64_t_ll_contains(R, 0));
    TEST_ASSERT_EQUAL(1, uint64_t_ll_contains(R, 1));
    TEST_ASSERT_EQUAL(1, uint64_t_ll_contains(R, 3));
    TEST_ASSERT_EQUAL(1, uint64_t_ll_contains(R, 4));
    TEST_ASSERT_EQUAL(1, uint64_t_ll_contains(R, 5));
    TEST_ASSERT_EQUAL(1, uint64_t_ll_contains(R, 6));
    TEST_ASSERT_EQUAL(1, uint64_t_ll_contains(R, 7));
    TEST_ASSERT_EQUAL(1, uint64_t_ll_contains(R, 8));
    TEST_ASSERT_EQUAL(1, uint64_t_ll_contains(R, 9));

    uint64_t_ll_free((void **)&R);
    R = NULL;

    R = giggle_search(domain, root_id, 19, 20);
    TEST_ASSERT_EQUAL(NULL, R);

    uint64_t_ll_free((void **)&R);
    R = NULL;

    R = giggle_search(domain, root_id, 18, 20);
    TEST_ASSERT_EQUAL(2, R->len);
    TEST_ASSERT_EQUAL(1, uint64_t_ll_contains(R, 9));
    TEST_ASSERT_EQUAL(1, uint64_t_ll_contains(R, 7));

    uint64_t_ll_free((void **)&R);
    R = NULL;

    cache.destroy();
}
//}}}

//{{{ void test_wah_giggle_search(void)
void test_wah_giggle_search(void)
{
    ORDER = 7;
    struct simple_cache *sc = simple_cache_init(5, 1, NULL);
    uint32_t domain = 0;
    wah_giggle_set_data_handler();

    /*
     *  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19 20
     *  |0-------------------------|  |1----|     |9------| 
     *     |2----|  |3-------|
     *        |4-------|  |5----|     |7-------------------|
     *                    |6------------|
     *  |8----------------------------------------------|
     *
     */
    uint32_t root_id = 0;

    uint32_t r = giggle_insert(domain, &root_id, 1, 10, 0);
    r = giggle_insert(domain, &root_id, 11, 13, 1);
    r = giggle_insert(domain, &root_id, 2, 4, 2);
    r = giggle_insert(domain, &root_id, 5, 8, 3);
    r = giggle_insert(domain, &root_id, 3, 6, 4);
    r = giggle_insert(domain, &root_id, 7, 9, 5);
    r = giggle_insert(domain, &root_id, 7, 12, 6);

    r = giggle_insert(domain, &root_id, 11, 18, 7);
    r = giggle_insert(domain, &root_id, 1, 17, 8);
    r = giggle_insert(domain, &root_id, 15, 18, 9);

    /*
    7 13
      (NULL)    1(SA:0,8) 2(SA:2) 3(SA:4) 5(SA:3 SE:2) 
      (0 3 4 8) 7(SA:5,6 SE:4) 9(SE: 3) 10(SE:5) 11(SA:1,7 SE:0) 
      (8 6 1 7) 13(SE:6) 14(SE:1) 15(SA:9) 18(SE:8) 19(SE:7,9)
    */

    /*
     *  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19 20
     *  |0-------------------------|  |1----|     |9------| 
     *     |2----|  |3-------|
     *        |4-------|  |5----|     |7-------------------|
     *                    |6------------|
     *  |8----------------------------------------------|
     *
     */

    ///struct uint64_t_ll *R = (struct uint64_t_ll *)giggle_search(root, 2, 5);
    //uint32_t R_id = 
    uint8_t *w = giggle_search(domain, root_id, 2, 5);
    uint32_t *R = NULL, R_size;
    R_size = wah_get_ints(w, &R);
    TEST_ASSERT_EQUAL(5, R_size);
    TEST_ASSERT_EQUAL(0, R[0] - 1);
    TEST_ASSERT_EQUAL(2, R[1] - 1);
    TEST_ASSERT_EQUAL(3, R[2] - 1);
    TEST_ASSERT_EQUAL(4, R[3] - 1);
    TEST_ASSERT_EQUAL(8, R[4] - 1);

    free(w);
    w = NULL;
    free(R);
    R = NULL;

    w = giggle_search(domain, root_id, 5, 15);
    R_size = wah_get_ints(w, &R);
    TEST_ASSERT_EQUAL(9, R_size);
    TEST_ASSERT_EQUAL(0, R[0] - 1);
    TEST_ASSERT_EQUAL(1, R[1] - 1);
    TEST_ASSERT_EQUAL(3, R[2] - 1);
    TEST_ASSERT_EQUAL(4, R[3] - 1);
    TEST_ASSERT_EQUAL(5, R[4] - 1);
    TEST_ASSERT_EQUAL(6, R[5] - 1);
    TEST_ASSERT_EQUAL(7, R[6] - 1);
    TEST_ASSERT_EQUAL(8, R[7] - 1);
    TEST_ASSERT_EQUAL(9, R[8] - 1);

    free(w);
    w = NULL;
    free(R);
    R = NULL;

    w = giggle_search(domain, root_id, 19, 20);
    if (w != NULL) {
        R_size = wah_get_ints(w, &R);
        TEST_ASSERT_EQUAL(0, R_size);
        free(R);
        R = NULL;
        free(w);
        w = NULL;
    } else {
        TEST_ASSERT_EQUAL(NULL, w);
    }

    w = giggle_search(domain, root_id, 18, 20);
    R_size = wah_get_ints(w, &R);
    TEST_ASSERT_EQUAL(2, R_size);
    TEST_ASSERT_EQUAL(7, R[0] - 1);
    TEST_ASSERT_EQUAL(9, R[1] - 1);

    free(w);
    w = NULL;
    free(R);
    R = NULL;

    cache.destroy();
}
//}}}

//{{{void test_giggle_get_chrm_id(void)
void test_giggle_get_chrm_id(void)
{
    struct giggle_index *gi = giggle_init_index(3, "tmp_offset_index.idx");
    TEST_ASSERT_EQUAL(0, giggle_get_chrm_id(gi, "chr2"));
    TEST_ASSERT_EQUAL(1, giggle_get_chrm_id(gi, "chr1"));
    TEST_ASSERT_EQUAL(2, giggle_get_chrm_id(gi, "chr3"));
    TEST_ASSERT_EQUAL(0, giggle_get_chrm_id(gi, "chr2"));

    TEST_ASSERT_EQUAL(3, giggle_get_chrm_id(gi, "chr9"));
    TEST_ASSERT_EQUAL(4, giggle_get_chrm_id(gi, "chr8"));
    TEST_ASSERT_EQUAL(5, giggle_get_chrm_id(gi, "chr7"));
    TEST_ASSERT_EQUAL(3, giggle_get_chrm_id(gi, "chr9"));
    giggle_index_destroy(&gi);
    remove("tmp_offset_index.idx");
}
//}}}

//{{{void test_giggle_index_file(void)
void test_giggle_index_file(void)
{
    ORDER = 10;
    struct simple_cache *sc = simple_cache_init(1000, 30, NULL);
    uint64_t_ll_giggle_set_data_handler();
    struct giggle_index *gi = giggle_init_index(30, "tmp_offset_index.idx");
    char *file_name = "../data/1k.unsort.bed.gz";
    uint32_t r = giggle_index_file(gi, file_name);

    valid_giggle_index(gi);

    TEST_ASSERT_EQUAL(1000, r);
    //TEST_ASSERT_EQUAL(23, gi->chrm_index->num); 
    TEST_ASSERT_EQUAL(23, gi->chrm_idx->index->num); 
    TEST_ASSERT_EQUAL(23, gi->num); 

    uint32_t sizes[23] = {152,
                          44,
                          66,
                          41,
                          9,
                          7,
                          30,
                          43,
                          38,
                          23,
                          77,
                          64,
                          21,
                          18,
                          27,
                          50,
                          20,
                          40,
                          73,
                          80,
                          33,
                          32,
                          12};

    char *chrms[23] = {"1",
                       "10",
                       "11",
                       "12",
                       "13",
                       "14",
                       "15",
                       "16",
                       "17",
                       "18",
                       "19",
                       "2",
                       "20",
                       "21",
                       "22",
                       "3",
                       "4",
                       "5",
                       "6",
                       "7",
                       "8",
                       "9",
                       "X"};
    /*
      for c in `gunzip -c 1k.unsort.bed.gz | cut -f1 | sort | uniq`
      do 
        V=`gunzip -c 1k.unsort.bed.gz | grep -n -e "$c\t" | cut -d":" -f1 | tr '\n' ',' | sed -e "s/,$//"`
        echo "{$V}"
      done
     */
    uint32_t ids[23][152] = {
        //{{{
        {18,22,41,42,59,62,67,69,71,85,86,87,88,98,103,104,106,114,131,133,141,148,151,156,168,170,177,180,189,190,193,197,204,207,211,212,218,221,241,247,255,257,261,267,276,282,311,322,327,333,338,358,359,365,372,386,408,414,418,458,466,467,468,473,480,483,485,496,499,500,504,512,525,549,552,558,568,570,574,578,581,587,601,607,612,613,625,626,628,638,641,650,657,666,676,677,678,679,684,686,687,688,693,697,715,718,724,729,730,735,740,745,749,753,755,763,764,766,776,781,783,785,800,801,817,821,822,828,844,846,865,869,884,885,889,924,928,930,934,941,943,952,954,956,962,964,977,978,988,994,995,1000},
        {3,27,34,76,79,92,113,130,164,215,216,222,271,272,332,390,409,425,463,486,519,569,585,593,606,611,622,645,665,669,689,696,719,748,775,780,790,853,855,871,886,898,908,998},
        {1,2,55,90,91,117,118,128,129,144,147,217,244,290,317,323,326,341,342,347,383,394,413,435,436,514,515,521,528,540,544,551,557,573,590,619,633,639,660,682,702,709,714,731,738,746,760,765,774,804,816,833,834,837,856,858,863,868,873,874,904,912,913,950,957,973},
        {47,70,80,107,124,145,159,205,210,228,236,253,258,288,393,398,401,410,474,481,507,531,537,554,579,584,588,635,640,670,704,712,713,727,770,805,810,887,896,925,974},
        {56,82,158,194,396,426,809,849,909},
        {102,198,501,516,672,705,918},
        {5,37,100,105,196,279,382,387,400,442,487,518,526,562,572,597,621,651,708,710,761,769,778,806,824,826,857,953,965,999},
        {4,50,64,115,127,149,152,214,225,248,259,275,281,293,310,321,357,421,424,439,491,506,547,609,767,768,795,862,870,872,878,893,905,906,911,915,919,921,935,944,975,979,996},
        {29,43,96,123,171,191,203,226,264,302,331,335,348,349,420,438,471,510,610,627,636,658,692,758,779,789,814,851,861,892,894,916,917,923,931,971,980,989},
        {8,30,31,73,122,140,162,200,213,251,309,336,378,490,524,527,556,637,711,784,802,907,983},
        {6,7,25,39,48,49,61,77,78,84,108,110,116,126,136,137,160,165,166,185,232,260,270,277,294,295,303,307,313,315,325,345,350,356,368,369,370,375,381,385,419,422,443,454,457,465,475,495,548,560,577,598,618,620,630,685,690,706,716,717,737,747,771,794,808,811,812,845,859,891,926,932,947,958,960,961,968},
        {17,51,60,68,101,142,175,178,227,234,235,252,254,274,286,287,292,296,298,314,316,339,344,353,361,397,433,445,448,451,455,523,530,539,594,596,599,608,629,642,674,681,732,742,754,757,762,807,815,818,825,847,854,875,900,901,902,914,938,948,963,966,981,986},
        {23,35,65,157,280,389,395,403,444,494,563,567,605,673,721,725,803,850,883,903,936},
        {9,112,181,209,324,402,411,520,536,564,643,683,759,798,835,836,929,997},
        {28,66,81,94,150,169,238,245,273,284,304,377,392,456,502,604,698,791,793,797,820,838,839,882,895,991,992},
        {24,54,72,83,89,93,139,153,161,163,224,230,242,243,269,283,300,328,346,362,363,367,373,406,428,434,460,470,476,479,488,497,513,517,532,534,543,580,617,649,664,680,723,788,827,880,881,942,945,949},
        {53,58,97,167,188,352,376,431,446,522,555,589,614,632,662,703,726,751,939,951},
        {16,21,57,121,174,195,223,233,237,285,297,312,320,329,337,355,404,437,441,477,498,505,511,535,550,583,586,694,699,734,787,823,830,840,866,876,879,959,976,984},
        {14,15,33,36,44,74,99,135,138,143,154,172,173,182,183,184,199,202,220,265,266,268,278,306,308,330,371,374,391,412,417,427,447,452,469,472,478,509,538,541,542,545,553,566,576,582,592,603,623,644,653,668,675,700,720,728,733,736,744,752,756,772,777,799,843,848,860,877,888,940,969,970,982},
        {10,11,12,13,26,40,46,52,63,95,109,132,134,176,179,201,208,219,231,240,250,256,299,305,318,319,334,340,351,354,360,384,388,415,416,423,449,453,459,461,462,529,533,546,559,561,565,571,575,595,600,615,631,634,648,655,656,663,695,701,722,741,773,782,792,819,829,864,890,897,910,922,927,933,946,967,972,985,990,993},
        {19,20,32,45,75,111,120,146,186,229,262,263,364,399,430,432,450,464,482,489,493,508,591,602,616,624,647,661,667,671,739,813,867},
        {38,119,125,155,187,192,206,246,289,291,301,366,379,380,429,440,484,492,646,659,691,743,750,796,831,832,841,842,899,937,955,987},
        {239,249,343,405,407,503,652,654,707,786,852,920}};
    //}}}

    int j,k;
    for (j = 0; j < 23; ++j) {
        uint32_t i = giggle_get_chrm_id(gi, chrms[j]);
        struct uint64_t_ll *R = (struct uint64_t_ll *)
            giggle_search(i, gi->root_ids[i], 0, 3000000000);
        TEST_ASSERT_EQUAL(sizes[j], R->len); 
        for (k = 0; k < R->len; ++k) {
            TEST_ASSERT_EQUAL(1, uint64_t_ll_contains(R, ids[j][k]-1));
        }
        uint64_t_ll_free((void **)&R);
    }

    giggle_index_destroy(&gi);
    cache.destroy();
    remove("tmp_offset_index.idx");
}
//}}}

//{{{ void test_uint64_t_ll_giggle_query_region(void)
void test_uint64_t_ll_giggle_query_region(void)
{
    ORDER = 10;
    struct simple_cache *sc = simple_cache_init(1000, 30, NULL);
    uint64_t_ll_giggle_set_data_handler();
    struct giggle_index *gi = giggle_init_index(30, "tmp_offset_index.idx");
    char *file_name = "../data/1k.unsort.bed.gz";
    uint32_t ret = giggle_index_file(gi, file_name);

    struct uint64_t_ll *R = (struct uint64_t_ll *)giggle_query_region(gi,
                                                                      "11",
                                                                      1000,
                                                                      3000000);
     //tabix 1k.sort.bed.gz chr11:1000-3000000
     //chr11    575808  576604  .   1000    .   ...
     //chr11    2950239 2952321 .   1000    .   ...
    TEST_ASSERT_EQUAL(2, R->len);

    /*
    struct file_id_offset_pair *r = 
        (struct file_id_offset_pair *)
        unordered_list_get(gi->offset_index, R->head->val);
    */
    struct file_id_offset_pair r = 
            offset_index_get(gi->offset_idx, R->head->val);
        //gi->offset_index->vals[R->head->val];
    struct file_data *fd = file_index_get(gi->file_idx, r.file_id);
    struct input_file *i = input_file_init(fd->file_name);

    i->input_file_seek(i, r.offset);

    int chrm_len = 10;
    char *chrm = (char *)malloc(chrm_len*sizeof(char));
    uint32_t start, end;
    long offset;
    kstring_t line = {0, 0, NULL};
            
    int x = i->input_file_get_next_interval(i,
                                            &chrm,
                                            &chrm_len,
                                            &start,
                                            &end,
                                            &offset,
                                            &line);
    
    TEST_ASSERT_EQUAL(0, strcmp("11", chrm));
    TEST_ASSERT_EQUAL(575808, start);
    TEST_ASSERT_EQUAL(576604, end);

    
    //r = (struct file_id_offset_pair *)
            //unordered_list_get(gi->offset_index, R->head->next->val);
    //r = gi->offset_index->vals[R->head->next->val];
    r = offset_index_get(gi->offset_idx, R->head->next->val);
    i->input_file_seek(i, r.offset);
    x = i->input_file_get_next_interval(i,
                                        &chrm,
                                        &chrm_len,
                                        &start,
                                        &end,
                                        &offset,
                                        &line);

    TEST_ASSERT_EQUAL(0, strcmp("11", chrm));
    TEST_ASSERT_EQUAL(2950239, start);
    TEST_ASSERT_EQUAL(2952321, end);

    uint64_t_ll_free((void **)&R);

    free(chrm);
    if (line.s != NULL)
        free(line.s);
    input_file_destroy(&i);
    giggle_index_destroy(&gi);
    cache.destroy();
    remove("tmp_offset_index.idx");
}
//}}}

////{{{void test_giggle_query_region(void)
//void test_wah_giggle_query_region(void)
//{
//    ORDER = 10;
//    struct simple_cache *sc = simple_cache_init(1000, 30, NULL);
//    wah_giggle_set_data_handler();
//    struct giggle_index *gi = giggle_init_index(30);
//    char *file_name = "../data/1k.unsort.bed.gz";
//    uint32_t ret = giggle_index_file(gi, file_name);
//
//    uint8_t *R_bm = (uint8_t*)giggle_query_region(gi,
//                                                  "11",
//                                                  1000,
//                                                  3000000);
//    /*
//     * tabix 1k.sort.bed.gz chr11:1000-3000000
//     * chr11    575808  576604  .   1000    .   ...
//     * chr11    2950239 2952321 .   1000    .   ...
//     */
//
//    uint32_t *R = NULL;
//    uint32_t R_len = wah_get_ints(R_bm, &R);
//    TEST_ASSERT_EQUAL(2, R_len);
//
//    //struct file_id_offset_pair *r = 
//        //(struct file_id_offset_pair *)
//        //unordered_list_get(gi->offset_index, R[0] - 1);
//    struct file_id_offset_pair r = gi->offset_index->vals[R[0] - 1];
//
//    struct file_data *fd = unordered_list_get(gi->file_index,
//                                              r.file_id);
//
//    struct input_file *i = input_file_init(fd->file_name);
//
//    i->input_file_seek(i, r.offset);
//
//    int chrm_len = 10;
//    char *chrm = (char *)malloc(chrm_len*sizeof(char));
//    uint32_t start, end;
//    long offset;
//            
//    int x = i->input_file_get_next_interval(i,
//                                            &chrm,
//                                            &chrm_len,
//                                            &start,
//                                            &end,
//                                            &offset);
//
//    TEST_ASSERT_EQUAL(0, strcmp("11", chrm));
//    TEST_ASSERT_EQUAL(575808, start);
//    TEST_ASSERT_EQUAL(576604, end);
//
//    //r = (struct file_id_offset_pair *)
//            //unordered_list_get(gi->offset_index, R[1] - 1);
//    r = gi->offset_index->vals[R[1] - 1];
//    i->input_file_seek(i, r.offset);
//    x = i->input_file_get_next_interval(i,
//                                        &chrm,
//                                        &chrm_len,
//                                        &start,
//                                        &end,
//                                        &offset);
//
//    TEST_ASSERT_EQUAL(0, strcmp("11", chrm));
//    TEST_ASSERT_EQUAL(2950239, start);
//    TEST_ASSERT_EQUAL(2952321, end);
//
//
//    free(R);
//    free(R_bm);
//    free(chrm);
//    input_file_destroy(&i);
//    giggle_index_destroy(&gi);
//    cache.destroy();
//}
////}}}

//{{{void test_giggle_index_directory(void)
void test_giggle_index_directory(void)
{
    ORDER = 10;
    struct simple_cache *sc = simple_cache_init(1000, 30, NULL);
    uint64_t_ll_giggle_set_data_handler();
    struct giggle_index *gi = giggle_init_index(30, "tmp_offset_index.idx");
    char *path_name = "../data/many/*bed.gz";
    uint32_t r = giggle_index_directory(gi, path_name, 0);

    valid_giggle_index(gi);

    TEST_ASSERT_EQUAL(21024, r);
    TEST_ASSERT_EQUAL(22, gi->file_idx->index->num);
    TEST_ASSERT_EQUAL(21024, gi->offset_idx->index->num);

    struct uint64_t_ll *R = (struct uint64_t_ll *)
        giggle_query_region(gi, "11", 1000, 3000000);

    uint64_t_ll_free((void **)&R);

    R = (struct uint64_t_ll *)giggle_query_region(gi,
                                                  "1",
                                                  1000,
                                                  3000000);
    TEST_ASSERT_EQUAL(71, R->len);
    uint64_t_ll_free((void **)&R);
    
    R = (struct uint64_t_ll *)giggle_query_region(gi,
                                                  "1",
                                                  249250622,
                                                  249250625);
    TEST_ASSERT_EQUAL(NULL, R);

    R = (struct uint64_t_ll *)giggle_query_region(gi,
                                                  "12",
                                                  52463173,
                                                  52464215);

    TEST_ASSERT_EQUAL(4, R->len);
    uint64_t_ll_free((void **)&R);
    giggle_index_destroy(&gi);
    cache.destroy();
    remove("tmp_offset_index.idx");
}
//}}}

//{{{ void test_giggle_init_store_load(void)
void test_giggle_init_store_load(void)
{
    struct giggle_index *gi = giggle_init(
                24,
                "tmp",
                1,
                uint64_t_ll_giggle_set_data_handler);

    char *path_name = "../data/many/*bed.gz";
    uint32_t r = giggle_index_directory(gi, path_name, 0);

    TEST_ASSERT_EQUAL(21024, r);
    TEST_ASSERT_EQUAL(22, gi->file_idx->index->num);
    TEST_ASSERT_EQUAL(21024, gi->offset_idx->index->num);

    struct uint64_t_ll *R = (struct uint64_t_ll *)
        giggle_query_region(gi, "11", 1000, 3000000);

    uint64_t_ll_free((void **)&R);

    R = (struct uint64_t_ll *)giggle_query_region(gi,
                                                  "1",
                                                  1000,
                                                  3000000);
    /*
     * ls *gz | xargs -I{} tabix {} chr1:1000-3000000 | wc -l
     * 71
     */
    TEST_ASSERT_EQUAL(71, R->len);

    uint64_t_ll_free((void **)&R);
    TEST_ASSERT_EQUAL(0, giggle_store(gi));
    giggle_index_destroy(&gi);
    cache.destroy();

    gi = giggle_load("tmp",
                     uint64_t_ll_giggle_set_data_handler);

    R = (struct uint64_t_ll *)giggle_query_region(gi,
                                                   "1",
                                                   1000,
                                                   1000000);
    TEST_ASSERT_EQUAL(16, R->len);
    uint64_t_ll_free((void **)&R);

    R = (struct uint64_t_ll *)giggle_query_region(gi,
                                                  "chr9",
                                                  112989628,
                                                  112989630);

    uint64_t_ll_free((void **)&R);
    giggle_index_destroy(&gi);
    cache.destroy();
    rmrf("tmp");
}
//}}}

////{{{void test_giggle_index_store(void)
//void test_giggle_index_store(void)
//{
//    struct stat st = {0};
//    if (stat("tmp", &st) == -1) {
//            mkdir("tmp", 0700);
//    } else {
//        rmrf("tmp/");
//        mkdir("tmp", 0700);
//    }
//
//    char *cache_names[30];
//    uint32_t i, ret;
//    for (i = 0; i < 30; ++i) {
//        ret = asprintf(&(cache_names[i]), "tmp/cache.%u", i);
//    }
//
//    ORDER = 10;
//    struct simple_cache *sc = simple_cache_init(1000, 30, cache_names);
//    uint64_t_ll_giggle_set_data_handler();
//    struct giggle_index *gi = giggle_init_index(30);
//    char *path_name = "../data/many/*bed.gz";
//    uint32_t r = giggle_index_directory(gi, path_name, 0);
//
//    for (i = 0; i < 30; ++i) {
//        bpt_write_tree(i, gi->root_ids[i]);
//    }
//
//    for (i = 0; i < 30; ++i) {
//        free(cache_names[i]);
//    }
//
//
//    char *chrm_index_file_name = "tmp/chrm_index.dat";
//    char *file_index_file_name = "tmp/file_index.dat";
//    char *offset_index_file_name = "tmp/offset_index.dat";
//    char *root_ids_file_name = "tmp/root_ids.dat";
//
//    FILE *f = fopen(root_ids_file_name, "wb");
//
//    if (fwrite(&(gi->len), sizeof(uint32_t), 1, f) != 1)
//        err(EX_IOERR, "Error writing len for root_ids'%s'.",
//            root_ids_file_name);
//
//    if (fwrite(gi->root_ids, sizeof(uint32_t), gi->len, f) != gi->len)
//        err(EX_IOERR, "Error writing root_ids '%s'.",
//            root_ids_file_name);
//    fclose(f);
//
//    f = fopen(chrm_index_file_name, "wb");
//    ordered_set_store(gi->chrm_idx->index,
//                      f,
//                      chrm_index_file_name,
//                      str_uint_pair_store);
//    fclose(f);
//
//    gi->file_idx->file_name = strdup(file_index_file_name);
//    file_index_store(gi->file_idx);
//
//    /*
//    f = fopen(file_index_file_name, "wb");
//    unordered_list_store(gi->file_index,
//                         f,
//                         file_index_file_name,
//                         c_str_store);
//    fclose(f);
//    */
//
//    gi->offset_idx->file_name = strdup(offset_index_file_name);
//    offset_index_store(gi->offset_idx);
//    /*
//    f = fopen(offset_index_file_name, "wb");
//    if (fwrite(&(gi->offset_index->num),
//               sizeof(uint64_t),1, f) != 1)
//        err(EX_IOERR, "Error writing offset_index num to '%s'.",
//            gi->offset_index_file_name);
//    if (fwrite(gi->offset_index->vals, 
//               sizeof(struct file_id_offset_pair), 
//               gi->offset_index->num, f) != gi->offset_index->num)
//        err(EX_IOERR, "Error writing file_id offset pairs to '%s'.",
//            gi->offset_index_file_name);
//    fclose(f);
//    */
//
//    giggle_index_destroy(&gi);
//    cache.destroy();
//    rmrf("tmp");
//}
////}}}

//{{{void valid_giggle_index(struct giggle_index *gi)
void valid_giggle_index(struct giggle_index *gi)
{
    uint32_t i,j;
    for (i = 0; i < gi->num; ++i) {
        struct uint64_t_ll *current_ids = NULL;
        struct bpt_node *r = cache.get(i,
                                       gi->root_ids[i] - 1,
                                       &bpt_node_cache_handler);
        while (BPT_IS_LEAF(r) == 0) {
            uint32_t left_most_id = BPT_POINTERS(r)[0];
            r = cache.get(i,
                          left_most_id - 1,
                          &bpt_node_cache_handler);
        }

        while (1) {

            /*
            if ((BPT_LEADING(r) == 0) && (current_ids != NULL)) {
                struct uint64_t_ll_node *c = current_ids->head;
                while (c != NULL) {
                    fprintf(stderr, "%u\n", c->val);
                    c = c->next;
                }
                
            }
            */


            /*
            fprintf(stderr, "%u %u\n",
                    BPT_LEADING(r) == 0, current_ids==NULL);
            */

            TEST_ASSERT_TRUE(
                    ((BPT_LEADING(r) == 0) && (current_ids==NULL))
                    ||
                    ((BPT_LEADING(r) != 0) && (current_ids!=NULL))
                    );


            if (BPT_LEADING(r) != 0) {
                struct uint64_t_ll_bpt_leading_data *ld = 
                        cache.get(i,
                                  BPT_LEADING(r) - 1,
                                  &uint64_t_ll_leading_cache_handler);

                if (ld->B != NULL) {
                    struct uint64_t_ll_node *curr = ld->B->head;
                    while (curr != NULL) {
                        TEST_ASSERT_EQUAL(1,
                                uint64_t_ll_contains(current_ids, curr->val));
                        curr = curr->next;
                    }
                }
            }


            for (j = 0; j < BPT_NUM_KEYS(r); ++j) {
                struct uint64_t_ll_bpt_non_leading_data *nld = 
                        cache.get(i,
                                  BPT_POINTERS(r)[j] - 1,
                                  &uint64_t_ll_non_leading_cache_handler);

                if (nld->SA != NULL) {
                    struct uint64_t_ll_node *curr = nld->SA->head;
                    while (curr != NULL) {
                        uint64_t_ll_uniq_append(&current_ids, curr->val);
                        curr = curr->next;
                    }
                }

                if (nld->SE != NULL) {
                    struct uint64_t_ll_node *curr = nld->SE->head;
                    while (curr != NULL) {
                        uint64_t_ll_remove(&current_ids, curr->val);
                        curr = curr->next;
                    }
                }
            }

            if (BPT_NEXT(r) == 0)
                break;
            r = cache.get(i,
                          BPT_NEXT(r) - 1,
                          &bpt_node_cache_handler);
        }
        TEST_ASSERT_EQUAL(NULL, current_ids);
    }
}
//}}}

//{{{ void test_valid_giggle_index_many(void)
void test_valid_giggle_index_many(void)
{
    ORDER=100;
    struct giggle_index *gi = giggle_init(24,
                                          "tmp",
                                          1,
                                          uint64_t_ll_giggle_set_data_handler);

    uint32_t total = 0;

    char *files[21] = { "../data/many/0.1.bed.gz",
                        "../data/many/0.2.bed.gz",
                        "../data/many/0.bed.gz",
                        "../data/many/1.1.bed.gz",
                        "../data/many/1.2.bed.gz",
                        "../data/many/1.bed.gz",
                        "../data/many/10.bed.gz",
                        "../data/many/2.1.bed.gz",
                        "../data/many/2.2.bed.gz",
                        "../data/many/2.bed.gz",
                        "../data/many/3.1.bed.gz",
                        "../data/many/3.2.bed.gz",
                        "../data/many/3.bed.gz",
                        "../data/many/4.1.bed.gz",
                        "../data/many/4.bed.gz",
                        "../data/many/5.1.bed.gz",
                        "../data/many/5.bed.gz",
                        "../data/many/6.bed.gz",
                        "../data/many/7.bed.gz",
                        "../data/many/8.bed.gz",
                        "../data/many/9.bed.gz"};
            
    uint32_t i;
    for (i = 0; i < 21; ++i) {
        {
            char *file_name = files[i];
            //fprintf(stderr, "%s\n", file_name);
            struct input_file *i = input_file_init(file_name);
            int chrm_len = 10;
            char *chrm = (char *)malloc(chrm_len*sizeof(char));
            uint32_t start, end;
            long offset;
            kstring_t line = {0, 0, NULL};

        
            /*
            struct file_data *fd = (struct file_data *)
                calloc(1, sizeof(struct file_data));
            fd->file_name = strdup(file_name);
       
            uint32_t file_id = unordered_list_add(gi->file_index, fd);
            */
            uint32_t file_id = file_index_add(gi->file_idx, file_name);
            struct file_data *fd = file_index_get(gi->file_idx, file_id);

            uint32_t j = 0;

            struct file_id_offset_pair *p;
            uint32_t intrv_id;

            while (i->input_file_get_next_interval(i,
                                                   &chrm,
                                                   &chrm_len,
                                                   &start,
                                                   &end,
                                                   &offset,
                                                   &line) >= 0) {
                intrv_id = offset_index_add(gi->offset_idx,
                                            offset,
                                            &line,
                                            file_id);

                uint32_t chrm_id = giggle_get_chrm_id(gi, chrm);
                uint32_t r = giggle_insert(chrm_id,
                                           &(gi->root_ids[chrm_id]),
                                           start,
                                           end,
                                           intrv_id);

                //fprintf(stderr, "%s %u %u\n", chrm, start, end);
                valid_giggle_index(gi);

                fd->mean_interval_size += end-start;
                fd->num_intervals += 1;
                j += 1;
            }
            fd->mean_interval_size = fd->mean_interval_size/fd->num_intervals;

            input_file_destroy(&i);
            if (line.s != NULL)
                free(line.s);
            free(chrm);
        }
    }
    giggle_index_destroy(&gi);
    cache.destroy();

    rmrf("tmp");
}
//}}}

//{{{ void test_giggle_init_store_load(void)
void test_giggle_index_search_store_search(void)
{
    struct giggle_index *gi = giggle_init(
                24,
                "tmp",
                1,
                uint64_t_ll_giggle_set_data_handler);

    char *path_name = "../data/many/*bed.gz";
    uint32_t r = giggle_index_directory(gi, path_name, 0);


    struct uint64_t_ll *R = (struct uint64_t_ll *)
                giggle_query_region(gi,
                                    "chr9",
                                    112989628,
                                    112989630);
   
    //ls *gz | xargs -I{} tabix {} chr1:1000-3000000 | wc -l
    //39
    //TEST_ASSERT_EQUAL(39, R->len);
    
    uint64_t_ll_free((void **)&R);
    TEST_ASSERT_EQUAL(0, giggle_store(gi));
    giggle_index_destroy(&gi);
    cache.destroy();

    gi = giggle_load("tmp",
                     uint64_t_ll_giggle_set_data_handler);

    R = (struct uint64_t_ll *)giggle_query_region(gi,
                                                  "chr9",
                                                  112989628,
                                                  112989630);
    uint64_t_ll_free((void **)&R);
    giggle_index_destroy(&gi);
    cache.destroy();
    rmrf("tmp");
}
//}}}

//{{{ void test_giggle_init_store_load_block(void)
void test_giggle_index_search_store_search_block(void)
{
    struct giggle_index *gi = giggle_init(
                24,
                "tmp",
                1,
                uint64_t_ll_giggle_set_data_handler);

    giggle_data_handler.write_tree = giggle_write_tree_leaf_data;

    char *path_name = "../data/many/*bed.gz";
    uint32_t r = giggle_index_directory(gi, path_name, 0);


    uint32_t num_tests = 10000, *chrs, *starts, *ends, *hits;
    chrs = (uint32_t *)malloc(num_tests * sizeof(uint32_t));
    starts = (uint32_t *)malloc(num_tests * sizeof(uint32_t));
    ends = (uint32_t *)malloc(num_tests * sizeof(uint32_t));
    hits = (uint32_t *)malloc(num_tests * sizeof(uint32_t));

    struct giggle_query_result *gqr;
    uint32_t count, i;
    for(i = 0; i < num_tests; i++) {
        chrs[i] = (rand() % 20) + 1;
        starts[i] = (rand()%10000000)+1;
        ends[i] = starts[i] + (rand()%10000000)+1;

        char *c;
        asprintf(&c, "%u", chrs[i]);
            

        gqr = giggle_query(gi,
                           c, 
                           starts[i],
                           ends[i],
                           NULL);

        free(c);

        count = 0;
        uint32_t j;
        for(j = 0; j < gqr->num_files; j++) 
            count += giggle_get_query_len(gqr, j);
        hits[i] = count;
        giggle_query_result_destroy(&gqr);
    }
   
    TEST_ASSERT_EQUAL(0, giggle_store(gi));
    giggle_index_destroy(&gi);
    cache.destroy();

    gi = giggle_load("tmp",
                     uint64_t_ll_giggle_set_data_handler);

    giggle_data_handler.giggle_collect_intersection =
            giggle_collect_intersection_data_in_block;

    giggle_data_handler.map_intersection_to_offset_list =
            leaf_data_map_intersection_to_offset_list;

    for(i = 0; i < num_tests; i++) {
        char *c;
        asprintf(&c, "%u", chrs[i]);
            
        gqr = giggle_query(gi,
                           c, 
                           starts[i],
                           ends[i],
                           NULL);


        count = 0;
        uint32_t j;
        for(j = 0; j < gqr->num_files; j++) 
            count += giggle_get_query_len(gqr, j);
        TEST_ASSERT_EQUAL(hits[i], count);
        giggle_query_result_destroy(&gqr);
        free(c);
    }

    free(chrs);
    free(starts);
    free(ends);
    free(hits);
 
    giggle_index_destroy(&gi);
    cache.destroy();
    rmrf("tmp");
}
//}}}

//{{{ void test_giggle_query_bug_0(void)
void test_giggle_query_bug_0(void)
{
    struct giggle_index *gi = giggle_init(
                24,
                "tmp",
                1,
                uint64_t_ll_giggle_set_data_handler);

    giggle_data_handler.write_tree = giggle_write_tree_leaf_data;

    char *path_name = "../data/many/*bed.gz";
    uint32_t r = giggle_index_directory(gi, path_name, 0);

    valid_giggle_index(gi);

    struct giggle_query_result *gqr;

    gqr = giggle_query(gi,
                       "18",
                       23669300,
                       23671590,
                       NULL);

    uint32_t o_count = 0;
    uint32_t j;
    for(j = 0; j < gqr->num_files; j++) 
        o_count += giggle_get_query_len(gqr, j);


    giggle_query_result_destroy(&gqr);

    TEST_ASSERT_EQUAL(0, giggle_store(gi));
    giggle_index_destroy(&gi);
    cache.destroy();

    gi = giggle_load("tmp",
                     uint64_t_ll_giggle_set_data_handler);

    giggle_data_handler.giggle_collect_intersection =
            giggle_collect_intersection_data_in_block;

    giggle_data_handler.map_intersection_to_offset_list =
            leaf_data_map_intersection_to_offset_list;

    gqr = giggle_query(gi,
                       "18",
                       23669300,
                       23671590,
                       NULL);

    uint32_t u_count = 0;
    for(j = 0; j < gqr->num_files; j++) 
        u_count += giggle_get_query_len(gqr, j);

    TEST_ASSERT_EQUAL(o_count, u_count);

    giggle_query_result_destroy(&gqr);
    giggle_index_destroy(&gi);
    cache.destroy();

    rmrf("tmp");
}
//}}}

//{{{ void test_giggle_query_bug_1(void)
void test_giggle_query_bug_1(void)
{
    struct giggle_index *gi = giggle_init(
                24,
                "tmp",
                1,
                uint64_t_ll_giggle_set_data_handler);

    giggle_data_handler.write_tree = giggle_write_tree_leaf_data;

    char *path_name = "../data/many/*bed.gz";
    uint32_t r = giggle_index_directory(gi, path_name, 0);

    valid_giggle_index(gi);

    struct giggle_query_result *gqr;

    gqr = giggle_query(gi,
                       "1",
                       19228888,
                       19229960,
                       NULL);

    uint32_t o_count = 0;
    uint32_t j;
    for(j = 0; j < gqr->num_files; j++) 
        o_count += giggle_get_query_len(gqr, j);


    giggle_query_result_destroy(&gqr);

    TEST_ASSERT_EQUAL(0, giggle_store(gi));
    giggle_index_destroy(&gi);
    cache.destroy();

    gi = giggle_load("tmp",
                     uint64_t_ll_giggle_set_data_handler);

    giggle_data_handler.giggle_collect_intersection =
            giggle_collect_intersection_data_in_block;

    giggle_data_handler.map_intersection_to_offset_list =
            leaf_data_map_intersection_to_offset_list;

    gqr = giggle_query(gi,
                       "1",
                       19228888,
                       19229960,
                       NULL);

    uint32_t u_count = 0;
    for(j = 0; j < gqr->num_files; j++) 
        u_count += giggle_get_query_len(gqr, j);

    TEST_ASSERT_EQUAL(o_count, u_count);

    giggle_query_result_destroy(&gqr);
    giggle_index_destroy(&gi);
    cache.destroy();
    rmrf("tmp");
}
//}}}

//{{{void test_giggle_bulk_insert(void)
void test_giggle_bulk_insert(void)
{
    char *input_path_name = "../data/many/*gz";
    char *output_path_name = "tmp_test_giggle_bulk_insert";

    uint64_t indexed_intervals = giggle_bulk_insert(input_path_name,
                                                    output_path_name,
                                                    1);

    char **names = NULL;
    uint32_t *num_intervals = NULL;
    double *mean_interval_sizes = NULL;
    uint32_t num_files = giggle_get_indexed_files(output_path_name,
                                                  &names,
                                                  &num_intervals,
                                                  &mean_interval_sizes);
    /*
     * ls ../data/many/ | wc -l
     * ls ../data/many/ 
     * ls ../data/many/  | xargs -I {} bash -c "gunzip -c {} | wc -l"
     * ls ../data/many/  | xargs -I {} bash -c "gunzip -c {} | awk '{s += \$3-\$2} END {print s/NR;}'"
     */
    TEST_ASSERT_EQUAL(22, num_intervals);
    
    char *A_names[22] = {"../data/many/0.1.bed.gz",
                         "../data/many/0.2.bed.gz",
                         "../data/many/0.bed.gz",
                         "../data/many/1.1.bed.gz",
                         "../data/many/1.2.bed.gz",
                         "../data/many/1.bed.gz",
                         "../data/many/10.bed.gz",
                         "../data/many/2.1.bed.gz",
                         "../data/many/2.2.bed.gz",
                         "../data/many/2.bed.gz",
                         "../data/many/3.1.bed.gz",
                         "../data/many/3.2.bed.gz",
                         "../data/many/3.bed.gz",
                         "../data/many/4.1.bed.gz",
                         "../data/many/4.bed.gz",
                         "../data/many/5.1.bed.gz",
                         "../data/many/5.bed.gz",
                         "../data/many/6.bed.gz",
                         "../data/many/7.bed.gz",
                         "../data/many/8.bed.gz",
                         "../data/many/9.bed.gz",
                         "../data/many/full_span.bed.gz"};

    uint32_t A_num_intervals[22] = { 1000,
                                     1000,
                                     1000,
                                     1000,
                                     1000,
                                     1000,
                                     1000,
                                     1000,
                                     1000,
                                     1000,
                                     1000,
                                     1000,
                                     1000,
                                     1000,
                                     1000,
                                     1000,
                                     1000,
                                     1000,
                                     1000,
                                     1000,
                                     1000,
                                     24};
    double A_mean_interval_sizes[22] = { 346.222,
                                         346.222,
                                         346.222,
                                         417.531,
                                         417.531,
                                         417.531,
                                         40.38,
                                         416.819,
                                         416.819,
                                         416.819,
                                         393.17,
                                         393.17,
                                         393.17,
                                         318.024,
                                         318.024,
                                         368.418,
                                         368.418,
                                         204.66,
                                         524.595,
                                         283.845,
                                         302.331,
                                         128986557.833333};
    uint32_t i;
    for (i = 0; i < num_files; ++i) {
        TEST_ASSERT_EQUAL(0,strcmp(A_names[i], names[i]));
        TEST_ASSERT_EQUAL(A_num_intervals[i],num_intervals[i]);
        TEST_ASSERT_EQUAL(A_mean_interval_sizes[i],mean_interval_sizes[i]);
        free(names[i]);
    }
    free(names);
    free(num_intervals);
    free(mean_interval_sizes);

    rmrf(output_path_name);
}
//}}}
