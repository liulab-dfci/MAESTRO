#include <stdio.h>
#include <stdlib.h>
#include <time.h> 
#include <math.h>
#include <stdio.h>
#include <inttypes.h>
#include <string.h>
#include <stdbool.h>
#include <unistd.h>

#include "bpt.h"
#include "lists.h"
#include "disk_store.h"

#include "unity.h"


void setUp(void)
{
}

void tearDown(void)
{
}

//{{{ int uint32_t_cmp(const void *a, const void *b)
static int uint32_t_cmp(const void *a, const void *b)
{
    uint32_t _a = *((uint32_t *)a), _b = *((uint32_t *)b);

    if (_a < _b)
        return -1;
    else if (_a > _b)
        return 1;
    else

    return 0;
}
//}}}

//{{{ void test_b_search(void)
void test_b_search(void)
{
    uint32_t D[10] = {2, 3, 4, 6, 8, 10, 12, 14, 16, 18};
    uint32_t A[20] = {0, // 0
                      0, // 1
                      0, // 2
                      1, // 3
                      2, // 4
                      3, // 5
                      3, // 6
                      4, // 7
                      4, // 8
                      5, // 9
                      5, // 10
                      6, // 11
                      6, // 12
                      7, // 13
                      7, // 14
                      8, // 15
                      8, // 16
                      9, // 17
                      9, // 18
                      10}; // 19

    int i;
    for (i = 0; i < 20; ++i)
        TEST_ASSERT_EQUAL(A[i], b_search(i, D, 10));
}
//}}}

//{{{void test_bpt_node_macros(void)
void test_bpt_node_macros(void)
{
    ORDER = 4;
    struct bpt_node *node = (struct bpt_node *)malloc(sizeof(struct bpt_node));

    uint32_t data[18] = {1,  // id
                         0,  // parent
                         0,  // is_leaf
                         0,  // leading
                         0,  // next
                         2,  // num_keys
                         2, 3, 0, 0, 0,  // keys (ORDER = 4)
                         0, // pointer head
                         5, 6, 7, 0, 0, 0}; //pointers (ORDER = 4)
    node->data = data;

    TEST_ASSERT_EQUAL(node->data[0], BPT_ID(node));
    TEST_ASSERT_EQUAL(node->data[1], BPT_PARENT(node));
    TEST_ASSERT_EQUAL(node->data[2], BPT_IS_LEAF(node));
    TEST_ASSERT_EQUAL(node->data[3], BPT_LEADING(node));
    TEST_ASSERT_EQUAL(node->data[4], BPT_NEXT(node));
    TEST_ASSERT_EQUAL(node->data[5], BPT_NUM_KEYS(node));
    TEST_ASSERT_EQUAL(node->data + 6, BPT_KEYS(node));
    TEST_ASSERT_EQUAL(node->data + 6 + ORDER + 1 + 1, BPT_POINTERS(node));

    TEST_ASSERT_EQUAL(2, BPT_KEYS(node)[0]);
    TEST_ASSERT_EQUAL(3, BPT_KEYS(node)[1]);
    TEST_ASSERT_EQUAL(0, BPT_KEYS(node)[2]);
    TEST_ASSERT_EQUAL(0, BPT_KEYS(node)[3]);
    TEST_ASSERT_EQUAL(0, BPT_KEYS(node)[4]);

    TEST_ASSERT_EQUAL(5, BPT_POINTERS(node)[0]);
    TEST_ASSERT_EQUAL(6, BPT_POINTERS(node)[1]);
    TEST_ASSERT_EQUAL(7, BPT_POINTERS(node)[2]);
    TEST_ASSERT_EQUAL(0, BPT_POINTERS(node)[3]);
    TEST_ASSERT_EQUAL(0, BPT_POINTERS(node)[4]);
    TEST_ASSERT_EQUAL(0, BPT_POINTERS(node)[5]);

    BPT_KEYS(node)[3] = 10;
    TEST_ASSERT_EQUAL(10, BPT_KEYS(node)[3]);

    BPT_ID(node) = 2;
    TEST_ASSERT_EQUAL(2, BPT_ID(node));

    ORDER=5;

    uint32_t data2[20] = {1,
                          0,  // parent
                          0,  // is_leaf
                          0,  // leading
                          0,  // next
                          2,  // num_keys
                          2, 3, 0, 0, 0, 0,  // keys (ORDER = 4)
                          0, 
                          5, 6, 7, 0, 0, 0, 0}; //pointers (ORDER = 4)

    node->data = data2;

    TEST_ASSERT_EQUAL(node->data[0], BPT_ID(node));
    TEST_ASSERT_EQUAL(node->data[1], BPT_PARENT(node));
    TEST_ASSERT_EQUAL(node->data[2], BPT_IS_LEAF(node));
    TEST_ASSERT_EQUAL(node->data[3], BPT_LEADING(node));
    TEST_ASSERT_EQUAL(node->data[4], BPT_NEXT(node));
    TEST_ASSERT_EQUAL(node->data[5], BPT_NUM_KEYS(node));
    TEST_ASSERT_EQUAL(node->data + 6, BPT_KEYS(node));
    TEST_ASSERT_EQUAL(node->data + 6 + ORDER + 2, BPT_POINTERS(node));

    TEST_ASSERT_EQUAL(2, BPT_KEYS(node)[0]);
    TEST_ASSERT_EQUAL(3, BPT_KEYS(node)[1]);
    TEST_ASSERT_EQUAL(0, BPT_KEYS(node)[2]);
    TEST_ASSERT_EQUAL(0, BPT_KEYS(node)[3]);
    TEST_ASSERT_EQUAL(0, BPT_KEYS(node)[4]);

    TEST_ASSERT_EQUAL(5, BPT_POINTERS(node)[0]);
    TEST_ASSERT_EQUAL(6, BPT_POINTERS(node)[1]);
    TEST_ASSERT_EQUAL(7, BPT_POINTERS(node)[2]);
    TEST_ASSERT_EQUAL(0, BPT_POINTERS(node)[3]);
    TEST_ASSERT_EQUAL(0, BPT_POINTERS(node)[4]);
    TEST_ASSERT_EQUAL(0, BPT_POINTERS(node)[5]);

    free(node);
}
//}}}

//{{{ void test_bpt_new_node(void)
void test_bpt_new_node(void)
{
    struct simple_cache *sc = simple_cache_init(5, 1, NULL);

    uint32_t domain = 0;

    struct bpt_node *n = bpt_new_node(domain);
    TEST_ASSERT_EQUAL(1, BPT_ID(n));
    TEST_ASSERT_EQUAL(1, cache.seen(domain));

    struct bpt_node *r = cache.get(domain, 1 - 1, &bpt_node_cache_handler);
    TEST_ASSERT_EQUAL(n, r);

    r = cache.get(domain, 2 - 1, &bpt_node_cache_handler);
    TEST_ASSERT_EQUAL(NULL, r);

    struct bpt_node *l = bpt_new_node(domain);
    TEST_ASSERT_EQUAL(2, BPT_ID(l));
    TEST_ASSERT_EQUAL(2, cache.seen(0));
    
    r = cache.get(domain, 2 - 1, &bpt_node_cache_handler);
    TEST_ASSERT_EQUAL(l, r);

    r = cache.get(domain, 3 - 1, &bpt_node_cache_handler);
    TEST_ASSERT_EQUAL(NULL, r);

    cache.destroy();
}
//}}}

//{{{ void test_bpt_find_leaf(void)
void test_bpt_find_leaf(void)
{
    struct simple_cache *sc = simple_cache_init(5, 1, NULL);

    uint32_t domain = 0;

    struct bpt_node *root = bpt_new_node(domain);
    BPT_NUM_KEYS(root) = 1;
    BPT_KEYS(root)[0] = 9;

    struct bpt_node *n1 = bpt_new_node(domain);
    BPT_NUM_KEYS(n1) = 1;
    BPT_KEYS(n1)[0] = 5;

    struct bpt_node *l1 = bpt_new_node(domain);
    BPT_IS_LEAF(l1) = 1;
    BPT_NUM_KEYS(l1) = 4;
    BPT_KEYS(l1)[0] = 1;
    BPT_KEYS(l1)[1] = 2;
    BPT_KEYS(l1)[2] = 3;
    BPT_KEYS(l1)[3] = 4;

    struct bpt_node *l2 = bpt_new_node(domain);
    BPT_IS_LEAF(l2) = 1;
    BPT_NUM_KEYS(l2) = 4;
    BPT_KEYS(l2)[0] = 5;
    BPT_KEYS(l2)[1] = 6;
    BPT_KEYS(l2)[2] = 7;
    BPT_KEYS(l2)[3] = 8;

    struct bpt_node *l3 = bpt_new_node(domain);
    BPT_IS_LEAF(l3) = 1;
    BPT_NUM_KEYS(l3) = 4;
    BPT_KEYS(l3)[0] = 9;
    BPT_KEYS(l3)[1] = 10;
    BPT_KEYS(l3)[2] = 11;
    BPT_KEYS(l3)[3] = 12;

    BPT_NEXT(l1) = BPT_ID(l2);
    BPT_NEXT(l2) = BPT_ID(l3);


    BPT_POINTERS(n1)[0] = BPT_ID(l1);
    BPT_POINTERS(n1)[1] = BPT_ID(l2);
    BPT_PARENT(l1) = BPT_ID(n1);
    BPT_PARENT(l2) = BPT_ID(n1);

    BPT_POINTERS(root)[0] = BPT_ID(n1);
    BPT_POINTERS(root)[1] = BPT_ID(l3);

    BPT_PARENT(n1) = BPT_ID(root);
    BPT_PARENT(l3) = BPT_ID(root);

    int i;
    for (i = 1; i <= 4; ++i)
        TEST_ASSERT_EQUAL(BPT_ID(l1),
                          bpt_find_leaf(domain, BPT_ID(root), i));

    for (i = 5; i <= 8; ++i) 
        TEST_ASSERT_EQUAL(BPT_ID(l2),
                          bpt_find_leaf(domain, BPT_ID(root), i));

    for (i = 9; i <= 12; ++i) 
        TEST_ASSERT_EQUAL(BPT_ID(l3),
                          bpt_find_leaf(domain, BPT_ID(root), i));

    cache.destroy();
}
//}}}

//{{{void test_bpt_split_node_non_root_parent_has_room(void)
void test_bpt_split_node_non_root_parent_has_room(void)
{
    struct simple_cache *sc = simple_cache_init(5, 1, NULL);
    ORDER = 4;

    uint32_t domain = 0;
    
    /*
     * 9
     * 6
     * 1,2,3,4,5
     *
     * Then split
     *
     * 9
     * 3,6
     * 1,2 3,4,5
     */
    struct bpt_node *root = bpt_new_node(domain);
    BPT_NUM_KEYS(root) = 1;
    BPT_KEYS(root)[0] = 9;

    struct bpt_node *n1 = bpt_new_node(domain);
    BPT_NUM_KEYS(n1) = 1;
    BPT_KEYS(n1)[0] = 6;

    BPT_POINTERS(root)[0] = BPT_ID(n1);
    BPT_PARENT(n1) = BPT_ID(root);

    struct bpt_node *l1 = bpt_new_node(domain);
    BPT_IS_LEAF(l1) = 1;
    BPT_NUM_KEYS(l1) = 5;
    BPT_KEYS(l1)[0] = 1;
    BPT_KEYS(l1)[1] = 2;
    BPT_KEYS(l1)[2] = 3;
    BPT_KEYS(l1)[3] = 4;
    BPT_KEYS(l1)[4] = 5;

    BPT_POINTERS(n1)[0] = BPT_ID(l1);
    BPT_PARENT(l1) = BPT_ID(n1);

    uint32_t *V_0 = (uint32_t *) malloc(sizeof(uint32_t));
    *V_0 = 2;
    uint32_t *V_1 = (uint32_t *) malloc(sizeof(uint32_t));
    *V_1 = 3;
    uint32_t *V_2 = (uint32_t *) malloc(sizeof(uint32_t));
    *V_2 = 4;
    uint32_t *V_3 = (uint32_t *) malloc(sizeof(uint32_t));
    *V_3 = 5;
    uint32_t *V_4 = (uint32_t *) malloc(sizeof(uint32_t));
    *V_4 = 6;

    // Put the data in the cache
    TEST_ASSERT_EQUAL(4, cache.seen(domain) + 1);
    BPT_POINTERS(l1)[0] = cache.seen(domain) + 1;
    cache.add(domain,
              cache.seen(domain),
              V_0,
              sizeof(uint32_t),
              &uint32_t_cache_handler);

    BPT_POINTERS(l1)[1] = cache.seen(domain) + 1;
    cache.add(domain,
              cache.seen(domain),
              V_1, 
              sizeof(uint32_t),
              &uint32_t_cache_handler);

    BPT_POINTERS(l1)[2] = cache.seen(domain) + 1;
    cache.add(domain,
              cache.seen(domain),
              V_2,
              sizeof(uint32_t),
              &uint32_t_cache_handler);

    BPT_POINTERS(l1)[3] = cache.seen(domain) + 1;
    cache.add(domain,
              cache.seen(domain),
              V_3,
              sizeof(uint32_t),
              &uint32_t_cache_handler);


    BPT_POINTERS(l1)[4] = cache.seen(domain) + 1;
    cache.add(domain,
              cache.seen(domain),
              V_4,
              sizeof(uint32_t),
              &uint32_t_cache_handler);



    uint32_t lo_id, hi_id;
    int split;
    uint32_t root_id = bpt_split_node(domain,
                                      BPT_ID(root),
                                      BPT_ID(l1),
                                      &lo_id,
                                      &hi_id,
                                      &split,
                                      NULL);

    TEST_ASSERT_EQUAL(BPT_ID(l1), lo_id);
    TEST_ASSERT_EQUAL(BPT_NEXT(l1), hi_id);
    TEST_ASSERT_EQUAL(2, split);


    TEST_ASSERT_EQUAL(2, BPT_NUM_KEYS(l1));
    struct bpt_node *hi_node = cache.get(domain,
                                         hi_id - 1,
                                         &bpt_node_cache_handler);
    TEST_ASSERT_EQUAL(3, BPT_NUM_KEYS(hi_node));
    TEST_ASSERT_EQUAL(2, BPT_NUM_KEYS(n1));
    TEST_ASSERT_EQUAL(3, BPT_KEYS(n1)[0]);
    TEST_ASSERT_EQUAL(6, BPT_KEYS(n1)[1]);
    TEST_ASSERT_EQUAL(lo_id, BPT_POINTERS(n1)[0]);
    TEST_ASSERT_EQUAL(hi_id, BPT_POINTERS(n1)[1]);

    TEST_ASSERT_EQUAL(lo_id, bpt_find_leaf(domain, root_id, 1));
    TEST_ASSERT_EQUAL(lo_id, bpt_find_leaf(domain, root_id, 2));
    TEST_ASSERT_EQUAL(hi_id, bpt_find_leaf(domain, root_id, 3));
    TEST_ASSERT_EQUAL(hi_id, bpt_find_leaf(domain, root_id, 4));
    TEST_ASSERT_EQUAL(hi_id, bpt_find_leaf(domain, root_id, 5));

    cache.destroy();
}
//}}}

//{{{void test_bpt_split_node_root(void)
void test_bpt_split_node_root(void)
{
    struct simple_cache *sc = simple_cache_init(5, 1, NULL);
    uint32_t domain = 0;
    ORDER = 4;

    uint32_t *V_0 = (uint32_t *) malloc(sizeof(int));
    *V_0 = 2;
    uint32_t *V_1 = (uint32_t *) malloc(sizeof(int));
    *V_1 = 3;
    uint32_t *V_2 = (uint32_t *) malloc(sizeof(int));
    *V_2 = 4;
    uint32_t *V_3 = (uint32_t *) malloc(sizeof(int));
    *V_3 = 5;
    uint32_t *V_4 = (uint32_t *) malloc(sizeof(int));
    *V_4 = 6;

    /*
     * 1,2,3,4,5
     *
     * SPLIT
     *
     * 3
     * 1,2 3,4,5
     */

    struct bpt_node *root = bpt_new_node(domain);
    BPT_IS_LEAF(root) = 1;
    BPT_NUM_KEYS(root) = 5;
    BPT_KEYS(root)[0] = 1;
    BPT_KEYS(root)[1] = 2;
    BPT_KEYS(root)[2] = 3;
    BPT_KEYS(root)[3] = 4;
    BPT_KEYS(root)[4] = 5;

    BPT_POINTERS(root)[0] = cache.seen(domain) + 1;
    cache.add(domain,
            cache.seen(domain),
            V_0,
            sizeof(uint32_t),
            &uint32_t_cache_handler);

    BPT_POINTERS(root)[1] = cache.seen(domain) + 1;
    cache.add(domain,
              cache.seen(domain),
              V_1,
              sizeof(uint32_t),
              &uint32_t_cache_handler);

    BPT_POINTERS(root)[2] = cache.seen(domain) + 1;
    cache.add(domain,
              cache.seen(domain),
              V_2,
              sizeof(uint32_t),
              &uint32_t_cache_handler);

    BPT_POINTERS(root)[3] = cache.seen(domain) + 1;
    cache.add(domain,
              cache.seen(domain),
              V_3,
              sizeof(uint32_t),
              &uint32_t_cache_handler);

    BPT_POINTERS(root)[4] = cache.seen(domain) + 1;
    cache.add(domain,
              cache.seen(domain),
              V_4,
              sizeof(uint32_t),
              &uint32_t_cache_handler);

    uint32_t lo_id, hi_id;
    int split;
    uint32_t root_id = bpt_split_node(domain,
                                      BPT_ID(root),
                                      BPT_ID(root),
                                      &lo_id,
                                      &hi_id,
                                      &split,
                                      NULL);

    TEST_ASSERT_EQUAL(BPT_ID(root), lo_id);
    TEST_ASSERT_EQUAL(BPT_NEXT(root), hi_id);
    TEST_ASSERT_EQUAL(2, split);


    TEST_ASSERT_EQUAL(2, BPT_NUM_KEYS(root));
    TEST_ASSERT_EQUAL(1, BPT_IS_LEAF(root));
    TEST_ASSERT_EQUAL(root_id, BPT_PARENT(root));

    struct bpt_node *next_node = cache.get(domain,
                                           BPT_NEXT(root) - 1,
                                           &bpt_node_cache_handler);
    TEST_ASSERT_EQUAL(3, BPT_NUM_KEYS(next_node));
    TEST_ASSERT_EQUAL(root_id, BPT_PARENT(next_node));

    struct bpt_node *new_root = cache.get(domain,
                                          root_id - 1,
                                          &bpt_node_cache_handler);
    TEST_ASSERT_EQUAL(1, BPT_NUM_KEYS(new_root));
    TEST_ASSERT_EQUAL(3, BPT_KEYS(new_root)[0]);
    TEST_ASSERT_EQUAL(lo_id, BPT_POINTERS(new_root)[0]);
    TEST_ASSERT_EQUAL(hi_id, BPT_POINTERS(new_root)[1]);

    cache.destroy();
}
//}}}

//{{{void test_bpt_insert(void)
void test_bpt_insert(void)
{
    struct simple_cache *sc = simple_cache_init(5, 1, NULL);
    uint32_t domain = 0;

    /*
     * 2,3,4,5
     *
     * 4
     *  2,3 4,5,6
     *
     * 4,6
     *  2,3 4,5 6,7,8
     *
     * 4,6,8
     *  2,3 4,5 6,7 8,9,10
     *
     * 4,6,8,10
     *  2,3 4,5 6,7 8,9 10,11,12
     *
     * 8 
     *  4,6 8,10,12
     *   2,3 4,5 6,7 8,9 10,11 12,13,14
     */
    ORDER = 4;
 
    uint32_t *V_0 = (uint32_t *) malloc(sizeof(uint32_t));
    *V_0 = 1;
    uint32_t *V_1 = (uint32_t *) malloc(sizeof(uint32_t));
    *V_1 = 2;
    uint32_t *V_2 = (uint32_t *) malloc(sizeof(uint32_t));
    *V_2 = 3;
    uint32_t *V_3 = (uint32_t *) malloc(sizeof(uint32_t));
    *V_3 = 4;
    uint32_t *V_4 = (uint32_t *) malloc(sizeof(uint32_t));
    *V_4 = 5;
    uint32_t *V_5 = (uint32_t *) malloc(sizeof(uint32_t));
    *V_5 = 6;
    uint32_t *V_6 = (uint32_t *) malloc(sizeof(uint32_t));
    *V_6 = 7;
    uint32_t *V_7 = (uint32_t *) malloc(sizeof(uint32_t));
    *V_7 = 8;
    uint32_t *V_8 = (uint32_t *) malloc(sizeof(uint32_t));
    *V_8 = 9;
    uint32_t *V_9 = (uint32_t *) malloc(sizeof(uint32_t));
    *V_9 = 10;
    uint32_t *V_10 = (uint32_t *) malloc(sizeof(uint32_t));
    *V_10 = 11;
    uint32_t *V_11 = (uint32_t *) malloc(sizeof(uint32_t));
    *V_11 = 12;
    uint32_t *V_12 = (uint32_t *) malloc(sizeof(uint32_t));
    *V_12 = 13;

    uint32_t leaf_id;
    int pos;
    uint32_t V_id = cache.seen(domain) + 1;
    cache.add(domain,
              V_id - 1,
              V_0,
              sizeof(uint32_t),
              &uint32_t_cache_handler);

    //2
    uint32_t root_id = bpt_insert(domain,
                                  0,
                                  2,
                                  V_id,
                                  &uint32_t_cache_handler,
                                  &leaf_id,
                                  &pos);

    struct bpt_node *root = cache.get(domain,
                                      root_id - 1,
                                      &bpt_node_cache_handler);
    TEST_ASSERT_EQUAL(1, BPT_NUM_KEYS(root));
    TEST_ASSERT_EQUAL(2, BPT_KEYS(root)[0]);
    TEST_ASSERT_EQUAL(1, BPT_IS_LEAF(root));
    TEST_ASSERT_EQUAL(root_id, leaf_id);
    TEST_ASSERT_EQUAL(0, pos);

    // 4
    //  2,3 4,5,6
    V_id = cache.seen(domain) + 1;
    cache.add(domain,
              V_id - 1,
              V_1,
              sizeof(uint32_t),
              &uint32_t_cache_handler);
    root_id = bpt_insert(domain,
                         root_id,
                         3,
                         V_id,
                         &uint32_t_cache_handler,
                         &leaf_id,
                         &pos);
    TEST_ASSERT_EQUAL(root_id, leaf_id);
    TEST_ASSERT_EQUAL(1, pos);

    V_id = cache.seen(domain) + 1;
    cache.add(domain,
              V_id - 1,
              V_2,
              sizeof(uint32_t),
              &uint32_t_cache_handler);
    root_id = bpt_insert(domain,
                         root_id,
                         4,
                         V_id,
                         &uint32_t_cache_handler,
                         &leaf_id,
                         &pos);
    TEST_ASSERT_EQUAL(root_id, leaf_id);
    TEST_ASSERT_EQUAL(2, pos);

    V_id = cache.seen(domain) + 1;
    cache.add(domain,
              V_id - 1,
              V_3,
              sizeof(uint32_t),
              &uint32_t_cache_handler);
    root_id = bpt_insert(domain,
                         root_id,
                         5,
                         V_id,
                         &uint32_t_cache_handler,
                         &leaf_id,
                         &pos);
    TEST_ASSERT_EQUAL(root_id, leaf_id);
    TEST_ASSERT_EQUAL(3, pos);

    V_id = cache.seen(domain) + 1;
    cache.add(domain,
              V_id - 1,
              V_4,
              sizeof(uint32_t),
              &uint32_t_cache_handler);
    root_id = bpt_insert(domain,
                         root_id,
                         6,
                         V_id,
                         &uint32_t_cache_handler,
                         &leaf_id,
                         &pos);

    TEST_ASSERT_EQUAL(2, pos);

    TEST_ASSERT_FALSE(root_id == leaf_id);
    root = cache.get(domain, root_id - 1, &bpt_node_cache_handler);
    TEST_ASSERT_EQUAL(1, BPT_NUM_KEYS(root));
    TEST_ASSERT_EQUAL(4, BPT_KEYS(root)[0]);
    TEST_ASSERT_EQUAL(0, BPT_IS_LEAF(root));
    TEST_ASSERT_EQUAL(leaf_id, BPT_POINTERS(root)[1]);

    struct bpt_node *node = cache.get(domain,
                                      BPT_POINTERS(root)[0] - 1,
                                      &bpt_node_cache_handler);
    TEST_ASSERT_EQUAL(2, BPT_NUM_KEYS(node));
    TEST_ASSERT_EQUAL(2, BPT_KEYS(node)[0]);
    TEST_ASSERT_EQUAL(3, BPT_KEYS(node)[1]);
    TEST_ASSERT_EQUAL(1, BPT_IS_LEAF(node));

    node = cache.get(domain, BPT_NEXT(node) - 1, &bpt_node_cache_handler);
    TEST_ASSERT_EQUAL(3, BPT_NUM_KEYS(node));
    TEST_ASSERT_EQUAL(4, BPT_KEYS(node)[0]);
    TEST_ASSERT_EQUAL(5, BPT_KEYS(node)[1]);
    TEST_ASSERT_EQUAL(6, BPT_KEYS(node)[2]);
    TEST_ASSERT_EQUAL(1, BPT_IS_LEAF(node));

    // 4,6
    //  2,3 4,5 6,7,8
    V_id = cache.seen(domain) + 1;
    cache.add(domain,
              V_id - 1,
              V_5,
              sizeof(uint32_t),
              &uint32_t_cache_handler);
    root_id = bpt_insert(domain,
                         root_id,
                         7,
                         V_id,
                         &uint32_t_cache_handler,
                         &leaf_id,
                         &pos);

    V_id = cache.seen(domain) + 1;
    cache.add(domain,
              V_id - 1,
              V_6,
              sizeof(uint32_t),
              &uint32_t_cache_handler);
    root_id = bpt_insert(domain,
                         root_id,
                         8,
                         V_id,
                         &uint32_t_cache_handler,
                         &leaf_id,
                         &pos);

    root = cache.get(domain, root_id - 1, &bpt_node_cache_handler);
    TEST_ASSERT_EQUAL(2, BPT_NUM_KEYS(root));
    TEST_ASSERT_EQUAL(4, BPT_KEYS(root)[0]);
    TEST_ASSERT_EQUAL(6, BPT_KEYS(root)[1]);
    TEST_ASSERT_EQUAL(0, BPT_IS_LEAF(root));

    node = cache.get(domain,
                     BPT_POINTERS(root)[0] - 1,
                     &bpt_node_cache_handler);
    TEST_ASSERT_EQUAL(2, BPT_NUM_KEYS(node));
    TEST_ASSERT_EQUAL(2, BPT_KEYS(node)[0]);
    TEST_ASSERT_EQUAL(3, BPT_KEYS(node)[1]);
    TEST_ASSERT_EQUAL(1, BPT_IS_LEAF(node));

    TEST_ASSERT_EQUAL(BPT_POINTERS(root)[1], BPT_NEXT(node));

    node = cache.get(domain,
                     BPT_POINTERS(root)[1] - 1,
                     &bpt_node_cache_handler);
    TEST_ASSERT_EQUAL(2, BPT_NUM_KEYS(node));
    TEST_ASSERT_EQUAL(4, BPT_KEYS(node)[0]);
    TEST_ASSERT_EQUAL(5, BPT_KEYS(node)[1]);
    TEST_ASSERT_EQUAL(1, BPT_IS_LEAF(node));

    TEST_ASSERT_EQUAL(BPT_POINTERS(root)[2], BPT_NEXT(node));

    node = cache.get(domain,
                     BPT_POINTERS(root)[2] - 1,
                     &bpt_node_cache_handler);
    TEST_ASSERT_EQUAL(3, BPT_NUM_KEYS(node));
    TEST_ASSERT_EQUAL(6, BPT_KEYS(node)[0]);
    TEST_ASSERT_EQUAL(7, BPT_KEYS(node)[1]);
    TEST_ASSERT_EQUAL(8, BPT_KEYS(node)[2]);
    TEST_ASSERT_EQUAL(1, BPT_IS_LEAF(node));

    // 8
    //   4,6 8,10,12
    //     2,3 4,5 6,7 8,9 10,11 12,13,14
    V_id = cache.seen(domain) + 1;
    cache.add(domain,
              V_id - 1,
              V_7,
              sizeof(uint32_t),
              &uint32_t_cache_handler);
    root_id = bpt_insert(domain,
                         root_id,
                         9,
                         V_id,
                         &uint32_t_cache_handler,
                         &leaf_id,
                         &pos);

    V_id = cache.seen(domain) + 1;
    cache.add(domain,
              V_id - 1,
              V_8,
              sizeof(uint32_t),
              &uint32_t_cache_handler);
    root_id = bpt_insert(domain,
                         root_id,
                         10,
                         V_id,
                         &uint32_t_cache_handler,
                         &leaf_id,
                         &pos);

    V_id = cache.seen(domain) + 1;
    cache.add(domain,
              V_id - 1,
              V_9,
              sizeof(uint32_t),
              &uint32_t_cache_handler);
    root_id = bpt_insert(domain,
                         root_id,
                         11,
                         V_id,
                         &uint32_t_cache_handler,
                         &leaf_id,
                         &pos);

    V_id = cache.seen(domain) + 1;
    cache.add(domain,
              V_id - 1,
              V_10,
              sizeof(uint32_t),
              &uint32_t_cache_handler);
    root_id = bpt_insert(domain,
                         root_id,
                         12,
                         V_id,
                         &uint32_t_cache_handler,
                         &leaf_id,
                         &pos);

    V_id = cache.seen(domain) + 1;
    cache.add(domain,
              V_id - 1,
              V_11,
              sizeof(uint32_t),
              &uint32_t_cache_handler);
    root_id = bpt_insert(domain,
                         root_id,
                         13,
                         V_id,
                         &uint32_t_cache_handler,
                         &leaf_id,
                         &pos);

    V_id = cache.seen(domain) + 1;
    cache.add(domain,
              V_id - 1,
              V_12,
              sizeof(uint32_t),
              &uint32_t_cache_handler);
    root_id = bpt_insert(domain,
                         root_id,
                         14,
                         V_id,
                         &uint32_t_cache_handler,
                         &leaf_id,
                         &pos);

    // 8
    //   4,6 8,10,12
    //     2,3 4,5 6,7 8,9 10,11 12,13,14
    root = cache.get(domain,
                     root_id - 1,
                     &bpt_node_cache_handler);
    TEST_ASSERT_EQUAL(1, BPT_NUM_KEYS(root));
    TEST_ASSERT_EQUAL(8, BPT_KEYS(root)[0]);
    TEST_ASSERT_EQUAL(0, BPT_IS_LEAF(root));

    node = cache.get(domain,
                     BPT_POINTERS(root)[0] - 1,
                     &bpt_node_cache_handler);
    TEST_ASSERT_EQUAL(2, BPT_NUM_KEYS(node));
    TEST_ASSERT_EQUAL(4, BPT_KEYS(node)[0]);
    TEST_ASSERT_EQUAL(6, BPT_KEYS(node)[1]);
    TEST_ASSERT_EQUAL(0, BPT_IS_LEAF(node));

    node = cache.get(domain,
                     BPT_POINTERS(root)[1] - 1,
                     &bpt_node_cache_handler);
    TEST_ASSERT_EQUAL(3, BPT_NUM_KEYS(node));
    TEST_ASSERT_EQUAL(8, BPT_KEYS(node)[0]);
    TEST_ASSERT_EQUAL(10, BPT_KEYS(node)[1]);
    TEST_ASSERT_EQUAL(12, BPT_KEYS(node)[2]);
    TEST_ASSERT_EQUAL(0, BPT_IS_LEAF(node));

    node = cache.get(domain,
                     BPT_POINTERS(root)[0] - 1,
                     &bpt_node_cache_handler);
    struct bpt_node *leaf_1 = cache.get(domain,
                                        BPT_POINTERS(node)[0] - 1,
                                        &bpt_node_cache_handler);
    TEST_ASSERT_EQUAL(2, BPT_NUM_KEYS(leaf_1));
    TEST_ASSERT_EQUAL(2, BPT_KEYS(leaf_1)[0]);
    TEST_ASSERT_EQUAL(3, BPT_KEYS(leaf_1)[1]);
    TEST_ASSERT_EQUAL(1, BPT_IS_LEAF(leaf_1));

    node = cache.get(domain,
                     BPT_POINTERS(root)[0] - 1,
                     &bpt_node_cache_handler);
    struct bpt_node *leaf_2 = cache.get(domain,
                                        BPT_POINTERS(node)[1] - 1,
                                        &bpt_node_cache_handler);
    TEST_ASSERT_EQUAL(2, BPT_NUM_KEYS(leaf_2));
    TEST_ASSERT_EQUAL(4, BPT_KEYS(leaf_2)[0]);
    TEST_ASSERT_EQUAL(5, BPT_KEYS(leaf_2)[1]);
    TEST_ASSERT_EQUAL(1, BPT_IS_LEAF(leaf_2));

    node = cache.get(domain,
                     BPT_POINTERS(root)[0] - 1,
                     &bpt_node_cache_handler);
    struct bpt_node *leaf_3 = cache.get(domain,
                                        BPT_POINTERS(node)[2] - 1,
                                        &bpt_node_cache_handler);
    TEST_ASSERT_EQUAL(2, BPT_NUM_KEYS(leaf_3));
    TEST_ASSERT_EQUAL(6, BPT_KEYS(leaf_3)[0]);
    TEST_ASSERT_EQUAL(7, BPT_KEYS(leaf_3)[1]);
    TEST_ASSERT_EQUAL(1, BPT_IS_LEAF(leaf_3));

    node = cache.get(domain,
                     BPT_POINTERS(root)[1] - 1,
                     &bpt_node_cache_handler);
    struct bpt_node *leaf_4 = cache.get(domain,
                                        BPT_POINTERS(node)[1] - 1,
                                        &bpt_node_cache_handler);
    TEST_ASSERT_EQUAL(2, BPT_NUM_KEYS(leaf_4));
    TEST_ASSERT_EQUAL(8, BPT_KEYS(leaf_4)[0]);
    TEST_ASSERT_EQUAL(9, BPT_KEYS(leaf_4)[1]);
    TEST_ASSERT_EQUAL(1, BPT_IS_LEAF(leaf_4));

    node = cache.get(domain,
                     BPT_POINTERS(root)[1] - 1,
                     &bpt_node_cache_handler);
    struct bpt_node *leaf_5 = cache.get(domain,
                                        BPT_POINTERS(node)[2] - 1,
                                        &bpt_node_cache_handler);
    TEST_ASSERT_EQUAL(2, BPT_NUM_KEYS(leaf_5));
    TEST_ASSERT_EQUAL(10, BPT_KEYS(leaf_5)[0]);
    TEST_ASSERT_EQUAL(11, BPT_KEYS(leaf_5)[1]);
    TEST_ASSERT_EQUAL(1, BPT_IS_LEAF(leaf_5));

    node = cache.get(domain,
                     BPT_POINTERS(root)[1] - 1,
                     &bpt_node_cache_handler);
    struct bpt_node *leaf_6 = cache.get(domain,
                                        BPT_POINTERS(node)[3] - 1,
                                        &bpt_node_cache_handler);
    TEST_ASSERT_EQUAL(3, BPT_NUM_KEYS(leaf_6));
    TEST_ASSERT_EQUAL(12, BPT_KEYS(leaf_6)[0]);
    TEST_ASSERT_EQUAL(13, BPT_KEYS(leaf_6)[1]);
    TEST_ASSERT_EQUAL(14, BPT_KEYS(leaf_6)[2]);
    TEST_ASSERT_EQUAL(1, BPT_IS_LEAF(leaf_6));

    TEST_ASSERT_EQUAL(BPT_ID(leaf_2), BPT_NEXT(leaf_1));
    TEST_ASSERT_EQUAL(BPT_ID(leaf_3), BPT_NEXT(leaf_2));
    TEST_ASSERT_EQUAL(BPT_ID(leaf_4), BPT_NEXT(leaf_3));
    TEST_ASSERT_EQUAL(BPT_ID(leaf_5), BPT_NEXT(leaf_4));
    TEST_ASSERT_EQUAL(BPT_ID(leaf_6), BPT_NEXT(leaf_5));

    cache.destroy();
}
//}}}

//{{{ void test_bpt_insert_new_value(void)
void test_bpt_insert_new_value(void)
{
    /*
     * 8 
     *  4,6 8,10,12
     *   2,3 4,5 6,7 8,9 10,11 12,13,14
     */
    ORDER = 4;
 
    uint32_t *V_0 = (uint32_t *) malloc(sizeof(uint32_t));
    *V_0 = 1;
    uint32_t *V_1 = (uint32_t *) malloc(sizeof(uint32_t));
    *V_1 = 2;
    uint32_t *V_2 = (uint32_t *) malloc(sizeof(uint32_t));
    *V_2 = 3;
    uint32_t *V_3 = (uint32_t *) malloc(sizeof(uint32_t));
    *V_3 = 4;
    uint32_t *V_4 = (uint32_t *) malloc(sizeof(uint32_t));
    *V_4 = 5;
    uint32_t *V_5 = (uint32_t *) malloc(sizeof(uint32_t));
    *V_5 = 6;
    uint32_t *V_6 = (uint32_t *) malloc(sizeof(uint32_t));
    *V_6 = 7;
    uint32_t *V_7 = (uint32_t *) malloc(sizeof(uint32_t));
    *V_7 = 8;
    uint32_t *V_8 = (uint32_t *) malloc(sizeof(uint32_t));
    *V_8 = 9;
    uint32_t *V_9 = (uint32_t *) malloc(sizeof(uint32_t));
    *V_9 = 10;
    uint32_t *V_10 = (uint32_t *) malloc(sizeof(uint32_t));
    *V_10 = 11;
    uint32_t *V_11 = (uint32_t *) malloc(sizeof(uint32_t));
    *V_11 = 12;
    uint32_t *V_12 = (uint32_t *) malloc(sizeof(uint32_t));
    *V_12 = 13;

    struct simple_cache *sc = simple_cache_init(5, 1, NULL);
    uint32_t domain = 0;
   
    uint32_t leaf_id;
    int pos;
    uint32_t V_id, root_id;

    root_id = bpt_insert_new_value(domain, 
                                   0,
                                   2,
                                   V_0,
                                   &uint32_t_cache_handler,
                                   &V_id,
                                   &leaf_id,
                                   &pos);

    root_id = bpt_insert_new_value(domain,
                                   root_id,
                                   3,
                                   V_1,
                                   &uint32_t_cache_handler,
                                   &V_id,
                                   &leaf_id,
                                   &pos);

    root_id = bpt_insert_new_value(domain,
                                   root_id,
                                   4,
                                   V_2,
                                   &uint32_t_cache_handler,
                                   &V_id,
                                   &leaf_id,
                                   &pos);

    root_id = bpt_insert_new_value(domain,
                                   root_id,
                                   5,
                                   V_3,
                                   &uint32_t_cache_handler,
                                   &V_id,
                                   &leaf_id,
                                   &pos);

    root_id = bpt_insert_new_value(domain,
                                   root_id,
                                   6,
                                   V_4,
                                   &uint32_t_cache_handler,
                                   &V_id,
                                   &leaf_id,
                                   &pos);

    root_id = bpt_insert_new_value(domain,
                                   root_id,
                                   7,
                                   V_5,
                                   &uint32_t_cache_handler,
                                   &V_id,
                                   &leaf_id,
                                   &pos);

    root_id = bpt_insert_new_value(domain,
                                   root_id,
                                   8,
                                   V_6,
                                   &uint32_t_cache_handler,
                                   &V_id,
                                   &leaf_id,
                                   &pos);

    root_id = bpt_insert_new_value(domain,
                                   root_id,
                                   9,
                                   V_7,
                                   &uint32_t_cache_handler,
                                   &V_id,
                                   &leaf_id,
                                   &pos);

    root_id = bpt_insert_new_value(domain,
                                   root_id,
                                   10,
                                   V_8,
                                   &uint32_t_cache_handler,
                                   &V_id,
                                   &leaf_id,
                                   &pos);


    root_id = bpt_insert_new_value(domain,
                                   root_id,
                                   11,
                                   V_9,
                                   &uint32_t_cache_handler,
                                   &V_id,
                                   &leaf_id,
                                   &pos);


    root_id = bpt_insert_new_value(domain,
                                   root_id,
                                   12,
                                   V_10,
                                   &uint32_t_cache_handler,
                                   &V_id,
                                   &leaf_id,
                                   &pos);


    root_id = bpt_insert_new_value(domain,
                                   root_id,
                                   13,
                                   V_11,
                                   &uint32_t_cache_handler,
                                   &V_id,
                                   &leaf_id,
                                   &pos);


    root_id = bpt_insert_new_value(domain,
                                   root_id,
                                   14,
                                   V_12,
                                   &uint32_t_cache_handler,
                                   &V_id,
                                   &leaf_id,
                                   &pos);

    // 8
    //   4,6 8,10,12
    //     2,3 4,5 6,7 8,9 10,11 12,13,14
    struct bpt_node *root = cache.get(domain,
                                      root_id - 1,
                                      &bpt_node_cache_handler);
    TEST_ASSERT_EQUAL(1, BPT_NUM_KEYS(root));
    TEST_ASSERT_EQUAL(8, BPT_KEYS(root)[0]);
    TEST_ASSERT_EQUAL(0, BPT_IS_LEAF(root));

    struct bpt_node *node = cache.get(domain,
                                      BPT_POINTERS(root)[0] - 1,
                                      &bpt_node_cache_handler);
    TEST_ASSERT_EQUAL(2, BPT_NUM_KEYS(node));
    TEST_ASSERT_EQUAL(4, BPT_KEYS(node)[0]);
    TEST_ASSERT_EQUAL(6, BPT_KEYS(node)[1]);
    TEST_ASSERT_EQUAL(0, BPT_IS_LEAF(node));

    node = cache.get(domain,
                     BPT_POINTERS(root)[1] - 1,
                     &bpt_node_cache_handler);
    TEST_ASSERT_EQUAL(3, BPT_NUM_KEYS(node));
    TEST_ASSERT_EQUAL(8, BPT_KEYS(node)[0]);
    TEST_ASSERT_EQUAL(10, BPT_KEYS(node)[1]);
    TEST_ASSERT_EQUAL(12, BPT_KEYS(node)[2]);
    TEST_ASSERT_EQUAL(0, BPT_IS_LEAF(node));

    node = cache.get(domain,
                     BPT_POINTERS(root)[0] - 1,
                     &bpt_node_cache_handler);
    struct bpt_node *leaf_1 = cache.get(domain,
                                        BPT_POINTERS(node)[0] - 1,
                                        &bpt_node_cache_handler);
    TEST_ASSERT_EQUAL(2, BPT_NUM_KEYS(leaf_1));
    TEST_ASSERT_EQUAL(2, BPT_KEYS(leaf_1)[0]);
    TEST_ASSERT_EQUAL(3, BPT_KEYS(leaf_1)[1]);
    TEST_ASSERT_EQUAL(1, BPT_IS_LEAF(leaf_1));

    node = cache.get(domain,
                     BPT_POINTERS(root)[0] - 1,
                     &bpt_node_cache_handler);
    struct bpt_node *leaf_2 = cache.get(domain,
                                        BPT_POINTERS(node)[1] - 1,
                                        &bpt_node_cache_handler);
    TEST_ASSERT_EQUAL(2, BPT_NUM_KEYS(leaf_2));
    TEST_ASSERT_EQUAL(4, BPT_KEYS(leaf_2)[0]);
    TEST_ASSERT_EQUAL(5, BPT_KEYS(leaf_2)[1]);
    TEST_ASSERT_EQUAL(1, BPT_IS_LEAF(leaf_2));

    node = cache.get(domain,
                     BPT_POINTERS(root)[0] - 1,
                     &bpt_node_cache_handler);
    struct bpt_node *leaf_3 = cache.get(domain,
                                        BPT_POINTERS(node)[2] - 1,
                                        &bpt_node_cache_handler);
    TEST_ASSERT_EQUAL(2, BPT_NUM_KEYS(leaf_3));
    TEST_ASSERT_EQUAL(6, BPT_KEYS(leaf_3)[0]);
    TEST_ASSERT_EQUAL(7, BPT_KEYS(leaf_3)[1]);
    TEST_ASSERT_EQUAL(1, BPT_IS_LEAF(leaf_3));

    node = cache.get(domain,
                     BPT_POINTERS(root)[1] - 1,
                     &bpt_node_cache_handler);
    struct bpt_node *leaf_4 = cache.get(domain,
                                        BPT_POINTERS(node)[1] - 1,
                                        &bpt_node_cache_handler);
    TEST_ASSERT_EQUAL(2, BPT_NUM_KEYS(leaf_4));
    TEST_ASSERT_EQUAL(8, BPT_KEYS(leaf_4)[0]);
    TEST_ASSERT_EQUAL(9, BPT_KEYS(leaf_4)[1]);
    TEST_ASSERT_EQUAL(1, BPT_IS_LEAF(leaf_4));

    node = cache.get(domain,
                     BPT_POINTERS(root)[1] - 1,
                     &bpt_node_cache_handler);
    struct bpt_node *leaf_5 = cache.get(domain,
                                        BPT_POINTERS(node)[2] - 1,
                                        &bpt_node_cache_handler);
    TEST_ASSERT_EQUAL(2, BPT_NUM_KEYS(leaf_5));
    TEST_ASSERT_EQUAL(10, BPT_KEYS(leaf_5)[0]);
    TEST_ASSERT_EQUAL(11, BPT_KEYS(leaf_5)[1]);
    TEST_ASSERT_EQUAL(1, BPT_IS_LEAF(leaf_5));

    node = cache.get(domain,
                     BPT_POINTERS(root)[1] - 1,
                     &bpt_node_cache_handler);
    struct bpt_node *leaf_6 = cache.get(domain,
                                        BPT_POINTERS(node)[3] - 1,
                                        &bpt_node_cache_handler);
    TEST_ASSERT_EQUAL(3, BPT_NUM_KEYS(leaf_6));
    TEST_ASSERT_EQUAL(12, BPT_KEYS(leaf_6)[0]);
    TEST_ASSERT_EQUAL(13, BPT_KEYS(leaf_6)[1]);
    TEST_ASSERT_EQUAL(14, BPT_KEYS(leaf_6)[2]);
    TEST_ASSERT_EQUAL(1, BPT_IS_LEAF(leaf_6));

    TEST_ASSERT_EQUAL(BPT_ID(leaf_2), BPT_NEXT(leaf_1));
    TEST_ASSERT_EQUAL(BPT_ID(leaf_3), BPT_NEXT(leaf_2));
    TEST_ASSERT_EQUAL(BPT_ID(leaf_4), BPT_NEXT(leaf_3));
    TEST_ASSERT_EQUAL(BPT_ID(leaf_5), BPT_NEXT(leaf_4));
    TEST_ASSERT_EQUAL(BPT_ID(leaf_6), BPT_NEXT(leaf_5));

    cache.destroy();
}
//}}}

//{{{ void test_bpt_insert_new_value(void)
void test_bpt_rand_order_insert_new_value(void)
{
    /*
     * 8 
     *  4,6 8,10,12
     *   2,3 4,5 6,7 8,9 10,11 12,13,14
     */
    ORDER = 4;
 
    uint32_t *V_0 = (uint32_t *) malloc(sizeof(uint32_t));
    *V_0 = 1;
    uint32_t *V_1 = (uint32_t *) malloc(sizeof(uint32_t));
    *V_1 = 2;
    uint32_t *V_2 = (uint32_t *) malloc(sizeof(uint32_t));
    *V_2 = 3;
    uint32_t *V_3 = (uint32_t *) malloc(sizeof(uint32_t));
    *V_3 = 4;
    uint32_t *V_4 = (uint32_t *) malloc(sizeof(uint32_t));
    *V_4 = 5;
    uint32_t *V_5 = (uint32_t *) malloc(sizeof(uint32_t));
    *V_5 = 6;
    uint32_t *V_6 = (uint32_t *) malloc(sizeof(uint32_t));
    *V_6 = 7;
    uint32_t *V_7 = (uint32_t *) malloc(sizeof(uint32_t));
    *V_7 = 8;
    uint32_t *V_8 = (uint32_t *) malloc(sizeof(uint32_t));
    *V_8 = 9;
    uint32_t *V_9 = (uint32_t *) malloc(sizeof(uint32_t));
    *V_9 = 10;
    uint32_t *V_10 = (uint32_t *) malloc(sizeof(uint32_t));
    *V_10 = 11;
    uint32_t *V_11 = (uint32_t *) malloc(sizeof(uint32_t));
    *V_11 = 12;
    uint32_t *V_12 = (uint32_t *) malloc(sizeof(uint32_t));
    *V_12 = 13;

    struct simple_cache *sc = simple_cache_init(5, 1, NULL);
    uint32_t domain = 0;
   
    uint32_t l;
    int p;
    uint32_t v, r;

    r = bpt_insert_new_value(domain,
                             0,
                             4,
                             V_2,
                             &uint32_t_cache_handler,
                             &v,
                             &l,
                             &p);

    r = bpt_insert_new_value(domain,
                             r,
                             3,
                             V_1,
                             &uint32_t_cache_handler,
                             &v,
                             &l,
                             &p);

    r = bpt_insert_new_value(domain,
                             r,
                             2,
                             V_0,
                             &uint32_t_cache_handler,
                             &v,
                             &l,
                             &p);

    r = bpt_insert_new_value(domain,
                             r,
                             14,
                             V_12,
                             &uint32_t_cache_handler,
                             &v,
                             &l,
                             &p);

    r = bpt_insert_new_value(domain,
                             r,
                             5,
                             V_3,
                             &uint32_t_cache_handler,
                             &v,
                             &l,
                             &p);

    r = bpt_insert_new_value(domain,
                             r,
                             12,
                             V_10,
                             &uint32_t_cache_handler,
                             &v,
                             &l,
                             &p);

    r = bpt_insert_new_value(domain,
                             r,
                             7,
                             V_5,
                             &uint32_t_cache_handler,
                             &v,
                             &l,
                             &p);

    r = bpt_insert_new_value(domain,
                             r,
                             9,
                             V_7,
                             &uint32_t_cache_handler,
                             &v,
                             &l,
                             &p);

    r = bpt_insert_new_value(domain,
                             r,
                             8,
                             V_6,
                             &uint32_t_cache_handler,
                             &v,
                             &l,
                             &p);

    r = bpt_insert_new_value(domain,
                             r,
                             13,
                             V_11,
                             &uint32_t_cache_handler,
                             &v,
                             &l,
                             &p);

    r = bpt_insert_new_value(domain,
                             r,
                             11,
                             V_9,
                             &uint32_t_cache_handler,
                             &v,
                             &l,
                             &p);

    r = bpt_insert_new_value(domain,
                             r,
                             6,
                             V_4,
                             &uint32_t_cache_handler,
                             &v,
                             &l,
                             &p);

    r = bpt_insert_new_value(domain,
                             r,
                             10,
                             V_8,
                             &uint32_t_cache_handler,
                             &v,
                             &l,
                             &p);

    /*
     * 4,7,9,12
     * 2,3, 4,5,6 7,8, 9,10,11 12,13,14
     */

    struct bpt_node *root = cache.get(domain, r - 1, &bpt_node_cache_handler);
    TEST_ASSERT_EQUAL(4, BPT_NUM_KEYS(root));
    TEST_ASSERT_EQUAL(4, BPT_KEYS(root)[0]);
    TEST_ASSERT_EQUAL(7, BPT_KEYS(root)[1]);
    TEST_ASSERT_EQUAL(9, BPT_KEYS(root)[2]);
    TEST_ASSERT_EQUAL(12, BPT_KEYS(root)[3]);
    TEST_ASSERT_EQUAL(0, BPT_IS_LEAF(root));

    struct bpt_node *leaf_1 = cache.get(domain,
                                        BPT_POINTERS(root)[0] - 1,
                                        &bpt_node_cache_handler);
    TEST_ASSERT_EQUAL(2, BPT_NUM_KEYS(leaf_1));
    TEST_ASSERT_EQUAL(2, BPT_KEYS(leaf_1)[0]);
    TEST_ASSERT_EQUAL(3, BPT_KEYS(leaf_1)[1]);
    TEST_ASSERT_EQUAL(1, BPT_IS_LEAF(leaf_1));

    struct bpt_node *leaf_2 = cache.get(domain,
                                        BPT_POINTERS(root)[1] - 1,
                                        &bpt_node_cache_handler);
    TEST_ASSERT_EQUAL(3, BPT_NUM_KEYS(leaf_2));
    TEST_ASSERT_EQUAL(4, BPT_KEYS(leaf_2)[0]);
    TEST_ASSERT_EQUAL(5, BPT_KEYS(leaf_2)[1]);
    TEST_ASSERT_EQUAL(6, BPT_KEYS(leaf_2)[2]);
    TEST_ASSERT_EQUAL(1, BPT_IS_LEAF(leaf_2));

    struct bpt_node *leaf_3 = cache.get(domain,
                                        BPT_POINTERS(root)[2] - 1,
                                        &bpt_node_cache_handler);
    TEST_ASSERT_EQUAL(2, BPT_NUM_KEYS(leaf_3));
    TEST_ASSERT_EQUAL(7, BPT_KEYS(leaf_3)[0]);
    TEST_ASSERT_EQUAL(8, BPT_KEYS(leaf_3)[1]);
    TEST_ASSERT_EQUAL(1, BPT_IS_LEAF(leaf_3));

    struct bpt_node *leaf_4 = cache.get(domain,
                                        BPT_POINTERS(root)[3] - 1,
                                        &bpt_node_cache_handler);
    TEST_ASSERT_EQUAL(3, BPT_NUM_KEYS(leaf_4));
    TEST_ASSERT_EQUAL(9, BPT_KEYS(leaf_4)[0]);
    TEST_ASSERT_EQUAL(10, BPT_KEYS(leaf_4)[1]);
    TEST_ASSERT_EQUAL(11, BPT_KEYS(leaf_4)[2]);
    TEST_ASSERT_EQUAL(1, BPT_IS_LEAF(leaf_4));

    struct bpt_node *leaf_5 = cache.get(domain,
                                        BPT_POINTERS(root)[4] - 1,
                                        &bpt_node_cache_handler);
    TEST_ASSERT_EQUAL(3, BPT_NUM_KEYS(leaf_5));
    TEST_ASSERT_EQUAL(12, BPT_KEYS(leaf_5)[0]);
    TEST_ASSERT_EQUAL(13, BPT_KEYS(leaf_5)[1]);
    TEST_ASSERT_EQUAL(14, BPT_KEYS(leaf_5)[2]);
    TEST_ASSERT_EQUAL(1, BPT_IS_LEAF(leaf_5));

    TEST_ASSERT_EQUAL(BPT_ID(leaf_2), BPT_NEXT(leaf_1));
    TEST_ASSERT_EQUAL(BPT_ID(leaf_3), BPT_NEXT(leaf_2));
    TEST_ASSERT_EQUAL(BPT_ID(leaf_4), BPT_NEXT(leaf_3));
    TEST_ASSERT_EQUAL(BPT_ID(leaf_5), BPT_NEXT(leaf_4));

    cache.destroy();
}
//}}}

//{{{void test_bpt_insert_id_updated_bpt_node(void)
void test_bpt_insert_id_updated_bpt_node(void)
{
    int V[14] = {1,2,3,4,5,6,7,8,9,10,11,12,13,14};
    /*
    uint32_t *V[14];
    uint32_t i;
    for (i = 0; i < 14; ++i) {
        V[i] = (uint32_t *)calloc(1,sizeof(uint32_t));
        *(V[i]) = i;
    }
    */

    struct simple_cache *sc = simple_cache_init(5, 1, NULL);
    uint32_t domain = 0;
 
    uint32_t r = 0, l, v = 0, v_id;
    int p;

    //r = bpt_insert_new_value(0, 4, V_2, free_wrapper, &v, &l, &p);

    r = bpt_insert_new_value(domain,
                             r,
                             2,
                             (V + v++),
                             NULL,
                             &v_id,
                             &l,
                             &p);
    r = bpt_insert_new_value(domain,
                             r,
                             5,
                             (V + v++),
                             NULL,
                             &v_id,
                             &l,
                             &p);
    r = bpt_insert_new_value(domain,
                             r,
                             10,
                             (V + v++),
                             NULL,
                             &v_id,
                             &l,
                             &p);
    r = bpt_insert_new_value(domain,
                             r,
                             15,
                             (V + v++),
                             NULL,
                             &v_id,
                             &l,
                             &p);

    r = bpt_insert_new_value(domain,
                             r,
                             1,
                             (V + v++),
                             NULL,
                             &v_id,
                             &l,
                             &p);

    struct bpt_node *root = cache.get(domain, r - 1, &bpt_node_cache_handler);
    TEST_ASSERT_EQUAL(BPT_POINTERS(root)[0], l);

    cache.destroy();
}
//}}}

//{{{void test_find(void)
void test_find(void)
{
    /*
     *  -> 2
     *  -> 4
     *  -> 6
     *  -> 8
     *  2,4,6,8
     *
     *  -> 10 
     *  6
     *  2,4, 6,8,10
     *
     *  -> 12 
     *  6
     *  2,4, 6,8,10,12
     *
     *  -> 14 
     *  6,10
     *  2,4, 6,8  10,12,14
     *
     *  -> 16 
     *  6,10
     *  2,4, 6,8  10,12,14,16
     *
     *  -> 18 
     *  6,10,14
     *  2,4, 6,8  10,12 14,16,18
     *
     *  -> 20 
     *  6,10,14
     *  2,4, 6,8  10,12 14,16,18,20
     *
     *  -> 22 
     *  6,10,14,18
     *  2,4, 6,8  10,12 14,16 18,20,22
     *
     *  -> 24 
     *  6,10,14,18
     *  2,4, 6,8  10,12 14,16 18,20,22,24
     *
     *  -> 26 
     *  14
     *  6,10 14,18,22
     *  2,4, 6,8  10,12 14,16 18,20 22,24,26
     *
     */

    struct simple_cache *sc = simple_cache_init(5, 1, NULL);
    uint32_t domain = 0;

    int V[13] = {1,2,3,4,5,6,7,8,9,10,11,12,13};
    int v=0;

    uint32_t r_id = 0, l_id, v_id;
    int pos;

    int i;
    for (i = 0; i < 13; ++i)
        r_id = bpt_insert_new_value(domain,
                                    r_id,
                                    (i+1)*2,
                                    (V + v++),
                                    NULL,
                                    &v_id,
                                    &l_id,
                                    &pos);

    uint32_t r_i;
    int *r, pos_r, pos_r_i = 0;
    v=0;
    int POS_R[13] = {0,1,0,1,0,1,0,1,0,1,0,1,2};
    for (i = 1; i <= 13; ++i) {
        r_i = bpt_find(domain, r_id, &l_id, &pos_r, i);
        if (i % 2 != 0)
            TEST_ASSERT_EQUAL(0, r_i);
        else  {
            r = cache.get(domain, r_i - 1, NULL);
            TEST_ASSERT_EQUAL(i/2, *r);
            TEST_ASSERT_EQUAL(POS_R[pos_r_i], pos_r);
            pos_r_i += 1;
        }
    }

    cache.destroy();
}
//}}}

//{{{void test_split_repair(void)
void decrement_repair(uint32_t domain, struct bpt_node *a, struct bpt_node *b)
{
    uint32_t i;
    for (i = 0; i < BPT_NUM_KEYS(a); ++i)
        BPT_KEYS(a)[i] = BPT_KEYS(a)[i] - 1;

    for (i = 0; i < BPT_NUM_KEYS(b); ++i)
        BPT_KEYS(b)[i] = BPT_KEYS(b)[i] + 1;

}

void test_split_repair(void)
{
    /*
     * 2,3,4,5
     *
     * 4
     *  2,3 4,5,6
     *
     * 4,6
     *  2,3 4,5 6,7,8
     *
     * 4,6,8
     *  2,3 4,5 6,7 8,9,10
     *
     * 4,6,8,10
     *  2,3 4,5 6,7 8,9 10,11,12
     *
     * 8 
     *  4,6 8,10,12
     *   2,3 4,5 6,7 8,9 10,11 12,13,14
     */
    int V[14] = {1,2,3,4,5,6,7,8,9,10,11,12,13,14};
    int v=0;

    bpt_node_repair = decrement_repair;

    struct simple_cache *sc = simple_cache_init(5, 1, NULL);
    uint32_t domain = 0;

    uint32_t r_id = 0, l_id, v_id;
    int pos;

    int i;
    r_id = bpt_insert_new_value(domain,
                                r_id,
                                2,
                                (V + v++),
                                NULL,
                                &v_id,
                                &l_id,
                                &pos);


    struct bpt_node *root = cache.get(domain,
                                      r_id - 1,
                                      &bpt_node_cache_handler);
    TEST_ASSERT_EQUAL(1, BPT_NUM_KEYS(root));
    TEST_ASSERT_EQUAL(2, BPT_KEYS(root)[0]);
    TEST_ASSERT_EQUAL(1, BPT_IS_LEAF(root));

    //  2,3,4,5
    r_id = bpt_insert_new_value(domain,
                                r_id,
                                3,
                                (V + v++),
                                NULL,
                                &v_id,
                                &l_id,
                                &pos);
    r_id = bpt_insert_new_value(domain,
                                r_id,
                                4,
                                (V + v++),
                                NULL,
                                &v_id,
                                &l_id,
                                &pos);
    r_id = bpt_insert_new_value(domain,
                                r_id,
                                5,
                                (V + v++),
                                NULL,
                                &v_id,
                                &l_id,
                                &pos);

    root = cache.get(domain, r_id - 1, &bpt_node_cache_handler);
    TEST_ASSERT_EQUAL(4, BPT_NUM_KEYS(root));
    TEST_ASSERT_EQUAL(2, BPT_KEYS(root)[0]);
    TEST_ASSERT_EQUAL(3, BPT_KEYS(root)[1]);
    TEST_ASSERT_EQUAL(4, BPT_KEYS(root)[2]);
    TEST_ASSERT_EQUAL(5, BPT_KEYS(root)[3]);
    TEST_ASSERT_EQUAL(1, BPT_IS_LEAF(root));

    // 4
    //  2,3 4,5,6
    r_id = bpt_insert_new_value(domain,
                                r_id,
                                6,
                                (V + v++),
                                NULL,
                                &v_id,
                                &l_id,
                                &pos);

    // Using a silly split repair function just to test, 
    // -1 the keys in the left bpt_node and + those in the right

    // 5
    //  1,2 5,6,7
    root = cache.get(domain, r_id - 1, &bpt_node_cache_handler);
    TEST_ASSERT_EQUAL(1, BPT_NUM_KEYS(root));
    TEST_ASSERT_EQUAL(5, BPT_KEYS(root)[0]);
    TEST_ASSERT_EQUAL(0, BPT_IS_LEAF(root));

    struct bpt_node *node_l = cache.get(domain,
                                        BPT_POINTERS(root)[0] - 1,
                                        &bpt_node_cache_handler);
    struct bpt_node *node_r = cache.get(domain,
                                        BPT_POINTERS(root)[1] - 1,
                                        &bpt_node_cache_handler);
    TEST_ASSERT_EQUAL(2, BPT_NUM_KEYS(node_l));
    TEST_ASSERT_EQUAL(1, BPT_KEYS(node_l)[0]);
    TEST_ASSERT_EQUAL(2, BPT_KEYS(node_l)[1]);
    TEST_ASSERT_EQUAL(1, BPT_IS_LEAF(node_l));
    TEST_ASSERT_EQUAL(BPT_POINTERS(root)[1], BPT_NEXT(node_l));

    TEST_ASSERT_EQUAL(3, BPT_NUM_KEYS(node_r));
    TEST_ASSERT_EQUAL(5, BPT_KEYS(node_r)[0]);
    TEST_ASSERT_EQUAL(6, BPT_KEYS(node_r)[1]);
    TEST_ASSERT_EQUAL(7, BPT_KEYS(node_r)[2]);
    TEST_ASSERT_EQUAL(1, BPT_IS_LEAF(node_r));
    TEST_ASSERT_EQUAL(0, BPT_NEXT(node_r));

    bpt_node_repair = NULL;

    cache.destroy();
}
//}}}

//{{{ void test_rand_test(void)
void test_rand_test(void)
{
    bpt_node_repair = NULL;
    ORDER = 10;
    struct simple_cache *sc = simple_cache_init(5, 1, NULL);
    uint32_t domain = 0;

    uint32_t r_id = 0, l_id, v_id;
    int pos;

    uint32_t size = 1000;

    uint32_t *d = (uint32_t *)malloc(size * sizeof(uint32_t));
    uint32_t *v = (uint32_t *)malloc(size * sizeof(uint32_t));

    time_t t = time(NULL);
    srand(t);

    uint32_t i, j;
    for (i = 0; i < size; ++i) {
        d[i] = rand() % 100000;
        v[i] = d[i] /2;
        r_id = bpt_insert_new_value(domain,
                                    r_id,
                                    d[i],
                                    (v + i),
                                    NULL,
                                    &v_id,
                                    &l_id,
                                    &pos);
    }

    int pos_i;
    for (i = 0; i < size; ++i) {
        uint32_t r_i = bpt_find(domain, r_id, &l_id, &pos_i, d[i]);
        int *r = cache.get(domain, r_i - 1, NULL);
        TEST_ASSERT_EQUAL(d[i]/2, *r);
    }

    free(d);
    free(v);
    cache.destroy();
}
//}}}

//{{{ void test_bpt_insert_repeat(void)
void test_bpt_insert_repeat(void)
{
    bpt_node_repair = NULL;
    ORDER = 4;
    struct simple_cache *sc = simple_cache_init(5, 1, NULL);
    uint32_t domain = 0;

    int V[14] = {1,2,3,4,5,6,7,8,9,10,11,12,13,14};

    uint32_t r_id = 0, v_id, l_id;
    int pos;

    r_id = bpt_insert_new_value(domain,
                                r_id,
                                3,
                                V,
                                NULL,
                                &v_id,
                                &l_id,
                                &pos);
    r_id = bpt_insert_new_value(domain,
                                r_id,
                                5,
                                V + 1,
                                NULL,
                                &v_id,
                                &l_id,
                                &pos);
    r_id = bpt_insert_new_value(domain,
                                r_id,
                                4,
                                V + 2,
                                NULL,
                                &v_id,
                                &l_id,
                                &pos);
    r_id = bpt_insert_new_value(domain,
                                r_id,
                                6,
                                V + 3,
                                NULL,
                                &v_id,
                                &l_id,
                                &pos);
    r_id = bpt_insert_new_value(domain,
                                r_id,
                                2,
                                V + 4,
                                NULL,
                                &v_id,
                                &l_id,
                                &pos);

    r_id = bpt_insert_new_value(domain,
                                r_id,
                                2,
                                V + 5,
                                NULL,
                                &v_id,
                                &l_id,
                                &pos);
    r_id = bpt_insert_new_value(domain,
                                r_id,
                                3,
                                V + 6,
                                NULL,
                                &v_id,
                                &l_id,
                                &pos);
    r_id = bpt_insert_new_value(domain,
                                r_id,
                                4,
                                V + 7,
                                NULL,
                                &v_id,
                                &l_id,
                                &pos);
    r_id = bpt_insert_new_value(domain,
                                r_id,
                                5,
                                V + 8,
                                NULL,
                                &v_id,
                                &l_id,
                                &pos);
    r_id = bpt_insert_new_value(domain,
                                r_id,
                                6,
                                V + 9,
                                NULL,
                                &v_id,
                                &l_id,
                                &pos);


    uint32_t r_i = bpt_find(domain, r_id, &l_id, &pos, 2);
    int *r = cache.get(domain, r_i - 1, NULL);
    TEST_ASSERT_EQUAL(V[5], *r);

    r_i = bpt_find(domain, r_id, &l_id, &pos, 3);
    r = cache.get(domain, r_i - 1, NULL);
    TEST_ASSERT_EQUAL(V[6], *r);

    r_i = bpt_find(domain, r_id, &l_id, &pos, 4);
    r = cache.get(domain, r_i - 1, NULL);
    TEST_ASSERT_EQUAL(V[7], *r);

    r_i = bpt_find(domain, r_id, &l_id, &pos, 5);
    r = cache.get(domain, r_i - 1, NULL);
    TEST_ASSERT_EQUAL(V[8], *r);

    r_i = bpt_find(domain, r_id, &l_id, &pos, 6);
    r = cache.get(domain, r_i - 1, NULL);
    TEST_ASSERT_EQUAL(V[9], *r);

    cache.destroy();
}
//}}}

//{{{ void test_bpt_insert_repeat_append(void)

void append_sum(uint32_t domain,
                uint32_t new_value_id,
                uint32_t curr_value_id,
                struct cache_handler *handler)
{

    uint32_t *new_value = cache.get(domain, new_value_id - 1, handler);
    uint32_t *curr_value = cache.get(domain, curr_value_id - 1, handler);
    *curr_value = *curr_value + *new_value;
}

void test_bpt_insert_repeat_append(void)
{
    append = append_sum;
    bpt_node_repair = NULL;
    ORDER = 4;
    struct simple_cache *sc = simple_cache_init(5, 1, NULL);
    uint32_t domain = 0;

    int R[14] = {1,2,3,4,5,6,7,8,9,10,11,12,13,14};
    int V[14] = {1,2,3,4,5,6,7,8,9,10,11,12,13,14};

    uint32_t r_id = 0, v_id, l_id;
    int pos;

    r_id = bpt_insert_new_value(domain,
                                r_id,
                                3,
                                V,
                                NULL,
                                &v_id,
                                &l_id,
                                &pos);
    r_id = bpt_insert_new_value(domain,
                                r_id,
                                5,
                                V + 1,
                                NULL,
                                &v_id,
                                &l_id,
                                &pos);
    r_id = bpt_insert_new_value(domain,
                                r_id,
                                4,
                                V + 2,
                                NULL,
                                &v_id,
                                &l_id,
                                &pos);
    r_id = bpt_insert_new_value(domain,
                                r_id,
                                6,
                                V + 3,
                                NULL,
                                &v_id,
                                &l_id,
                                &pos);
    r_id = bpt_insert_new_value(domain,
                                r_id,
                                2,
                                V + 4,
                                NULL,
                                &v_id,
                                &l_id,
                                &pos);

    r_id = bpt_insert_new_value(domain,
                                r_id,
                                2,
                                V + 5,
                                NULL,
                                &v_id,
                                &l_id,
                                &pos);
    r_id = bpt_insert_new_value(domain,
                                r_id,
                                3,
                                V + 6,
                                NULL,
                                &v_id,
                                &l_id,
                                &pos);
    r_id = bpt_insert_new_value(domain,
                                r_id,
                                4,
                                V + 7,
                                NULL,
                                &v_id,
                                &l_id,
                                &pos);
    r_id = bpt_insert_new_value(domain,
                                r_id,
                                5,
                                V + 8,
                                NULL,
                                &v_id,
                                &l_id,
                                &pos);
    r_id = bpt_insert_new_value(domain,
                                r_id,
                                6,
                                V + 9,
                                NULL,
                                &v_id,
                                &l_id,
                                &pos);


    uint32_t r_i = bpt_find(domain, r_id, &l_id, &pos, 2);
    int *r = cache.get(domain, r_i - 1, NULL);
    TEST_ASSERT_EQUAL(R[5] + R[4], *r);

    r_i = bpt_find(domain, r_id, &l_id, &pos, 3);
    r = cache.get(domain, r_i - 1, NULL);
    TEST_ASSERT_EQUAL(R[6] + R[0], *r);

    r_i = bpt_find(domain, r_id, &l_id, &pos, 4);
    r = cache.get(domain, r_i - 1, NULL);
    TEST_ASSERT_EQUAL(R[7] + R[2], *r);

    r_i = bpt_find(domain, r_id, &l_id, &pos, 5);
    r = cache.get(domain, r_i - 1, NULL);
    TEST_ASSERT_EQUAL(R[8] + R[1], *r);

    r_i = bpt_find(domain, r_id, &l_id, &pos, 6);
    r = cache.get(domain, r_i - 1, NULL);
    TEST_ASSERT_EQUAL(R[9] + R[3], *r);

    cache.destroy();

    append = NULL;
}
//}}}

//{{{ void test_rand_test_high_order(void)
void test_rand_test_high_order(void)
{
    struct simple_cache *sc = simple_cache_init(5, 1, NULL);
    uint32_t domain = 0;
    ORDER=50;
    bpt_node_repair = NULL;
    uint32_t size = 100000;

    uint32_t *d = (uint32_t *)malloc(size * sizeof(uint32_t));
    uint32_t *v = (uint32_t *)malloc(size * sizeof(uint32_t));

    time_t t = time(NULL);
    srand(t);

    uint32_t r_id = 0, v_id, l_id;
    int pos;

 
    uint32_t i, j;
    for (i = 0; i < size; ++i) {
        d[i] = rand();
        v[i] = d[i] /2;
        r_id = bpt_insert_new_value(domain,
                                    r_id,
                                    d[i],
                                    v + i,
                                    NULL,
                                    &v_id,
                                    &l_id,
                                    &pos);
    }

    uint32_t *r;
    int pos_r;
    for (i = 0; i < size; ++i) {
        uint32_t r_i = bpt_find(domain, r_id, &l_id, &pos, d[i]);
        int *r = cache.get(domain, r_i - 1, NULL);
        TEST_ASSERT_EQUAL(v[i], *r);
    }

    free(d);
    free(v);
    cache.destroy();
}
//}}}

//{{{void test_bpt_write_tree(void)
void test_bpt_write_tree(void)
{
    char *bpt_file_out[1] = {"test_bpt_write_tree.out"};

    if (access("test_bpt_write_tree.out.idx", F_OK) != -1)
        remove("test_bpt_write_tree.out.idx");

    if (access("test_bpt_write_tree.out.dat", F_OK) != -1)
        remove("test_bpt_write_tree.out.dat");

    struct simple_cache *sc = simple_cache_init(5, 1, bpt_file_out);
    uint32_t domain = 0;
    ORDER=5;
    bpt_node_repair = NULL;
    uint32_t size = 20;

    uint32_t *keys = (uint32_t *)malloc(size * sizeof(uint32_t));
    //uint32_t *v = (uint32_t *)malloc(size * sizeof(uint32_t));

    time_t t = time(NULL);
    srand(t);

    uint32_t r_id = 0, v_id, l_id;
    int pos;

    uint32_t i, j;
    for (i = 0; i < size; ++i) {
        keys[i] = rand();
        uint32_t *value = (uint32_t *)calloc(1,sizeof(uint32_t));
        *value = keys[i] /2;
        r_id = bpt_insert_new_value(domain,
                                    r_id,
                                    keys[i],
                                    value,
                                    &uint32_t_cache_handler,
                                    &v_id,
                                    &l_id,
                                    &pos);
    }


    bool ret = bpt_write_tree(domain, r_id);

    cache.destroy();

    sc = simple_cache_init(5, 1, bpt_file_out);

    uint32_t *r;
    int pos_r;
    for (i = 0; i < size; ++i) {
        uint32_t r_i = bpt_find(domain, r_id, &l_id, &pos, keys[i]);
        uint32_t *r = cache.get(domain, r_i - 1, &uint32_t_cache_handler);
        TEST_ASSERT_EQUAL(keys[i] / 2, *r);
    }

    cache.destroy();
    remove("test_bpt_write_tree.out.idx");
    remove("test_bpt_write_tree.out.dat");

    /*
    struct disk_store *ds = disk_store_load(&f, bpt_file_out);

    struct bpt_node *read_root = (struct bpt_node *)
            malloc(sizeof(struct bpt_node));
    uint64_t node_size;
    read_root->data = disk_store_get(ds, 0, &node_size);

    struct bpt_node *mem_root = cache.get(cache.cache, r_id);

    TEST_ASSERT_EQUAL(1, BPT_ID(read_root));
    TEST_ASSERT_EQUAL(BPT_NUM_KEYS(mem_root), BPT_NUM_KEYS(read_root));

    for (i = 0; i < BPT_NUM_KEYS(mem_root); ++i) 
        TEST_ASSERT_EQUAL(BPT_KEYS(mem_root)[i], BPT_KEYS(read_root)[i]);

    for (i = 0; i <= BPT_NUM_KEYS(mem_root); ++i) {
        struct bpt_node *mem_node = cache.get(cache.cache,
                                              BPT_POINTERS(mem_root)[i]);
        struct bpt_node *read_node = (struct bpt_node *)
                malloc(sizeof(struct bpt_node));
        read_node->data = disk_store_get(ds,
                                         BPT_POINTERS(read_root)[i] - 1,
                                         &node_size);
        TEST_ASSERT_EQUAL(BPT_NUM_KEYS(mem_node), BPT_NUM_KEYS(read_node));

        for (j = 0; j < BPT_NUM_KEYS(mem_node); ++j) {
            TEST_ASSERT_EQUAL(BPT_KEYS(mem_node)[j], BPT_KEYS(read_node)[j]);

            uint32_t *mem_val = cache.get(cache.cache,
                                          BPT_POINTERS(mem_node)[j]);
            uint32_t *read_val = disk_store_get(ds,
                                                BPT_POINTERS(read_node)[j] - 1,
                                                &node_size);
            TEST_ASSERT_EQUAL(*mem_val, *read_val);
            free(read_val);
        }

        free(read_node->data);
        free(read_node);
    }

    free(read_root->data);
    free(read_root);

    free(d);
    free(v);
    cache.destroy(&(cache.cache));
    disk_store_destroy(&ds);
    */
}
//}}}

#if 0
////{{{ void test_bpt_destroy_tree(void)
//void test_bpt_destroy_tree(void)
//{
//    TEST_IGNORE()
//
//    /*
//     * 2,3,4,5
//     *
//     * 4
//     *  2,3 4,5,6
//     *
//     * 4,6
//     *  2,3 4,5 6,7,8
//     *
//     * 4,6,8
//     *  2,3 4,5 6,7 8,9,10
//     *
//     * 4,6,8,10
//     *  2,3 4,5 6,7 8,9 10,11,12
//     *
//     * 8 
//     *  4,6 8,10,12
//     *   2,3 4,5 6,7 8,9 10,11 12,13,14
//     */
//    int V[14] = {1,2,3,4,5,6,7,8,9,10,11,12,13,14};
//    int v=0;
//
//    struct bpt_node *root = NULL;
//    struct bpt_node *leaf = NULL;
//    int pos;
//
//    root = bpt_insert(root, 2, (void *)(V + v++), &leaf, &pos);
//    root = bpt_insert(root, 3, (void *)(V + v++), &leaf, &pos);
//    root = bpt_insert(root, 4, (void *)(V + v++), &leaf, &pos);
//    root = bpt_insert(root, 5, (void *)(V + v++), &leaf, &pos);
//    root = bpt_insert(root, 6, (void *)(V + v++), &leaf, &pos);
//    root = bpt_insert(root, 7, (void *)(V + v++), &leaf, &pos);
//    root = bpt_insert(root, 8, (void *)(V + v++), &leaf, &pos);
//    root = bpt_insert(root, 9, (void *)(V + v++), &leaf, &pos);
//    root = bpt_insert(root, 10, (void *)(V + v++), &leaf, &pos);
//    root = bpt_insert(root, 11, (void *)(V + v++), &leaf, &pos);
//    root = bpt_insert(root, 12, (void *)(V + v++), &leaf, &pos);
//    root = bpt_insert(root, 13, (void *)(V + v++), &leaf, &pos);
//    root = bpt_insert(root, 14, (void *)(V + v++), &leaf, &pos);
//
//    bpt_destroy_tree(&root);
//}
////}}}
////{{{ void test_bpt_find_null(void)
//void test_bpt_find_null(void)
//{
//    TEST_IGNORE()
//    int V[14] = {1,2,3,4,5,6,7,8,9,10,11,12,13,14};
//    int v=0;
//
//    struct bpt_node *root = NULL;
//    struct bpt_node *leaf = NULL;
//    int pos_r;
//
//    int *r = (int *)bpt_find(root, &leaf, &pos_r, 0);
//
//    TEST_ASSERT_EQUAL(NULL, r);
//
//    root = bpt_insert(root, 2, (void *)(V + v++), &leaf, &pos_r);
//
//    r = (int *)bpt_find(root, &leaf, &pos_r, 2);
//    TEST_ASSERT_EQUAL(V, r);
//
//    r = (int *)bpt_find(root, &leaf, &pos_r, 1);
//    TEST_ASSERT_EQUAL(NULL, r);
//
//    bpt_destroy_tree(&root);
//}
////}}}
#endif
