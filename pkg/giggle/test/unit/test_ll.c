#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <inttypes.h>
#include <string.h>

#include "unity.h"
#include "util.h"
#include "bpt.h"
#include "ll.h"

void setUp(void) { }
void tearDown(void) { }

//{{{void test_uint64_t_ll_append(void)
void test_uint64_t_ll_append(void)
{
    struct uint64_t_ll *ll = NULL;

    uint64_t_ll_append(&ll, 1);
    uint64_t_ll_append(&ll, 2);
    uint64_t_ll_append(&ll, 5);

    TEST_ASSERT_EQUAL(1, ll->head->val);
    TEST_ASSERT_EQUAL(2, ll->head->next->val);
    TEST_ASSERT_EQUAL(5, ll->head->next->next->val);

    uint64_t_ll_free((void **)&ll);

    TEST_ASSERT_EQUAL(NULL, ll);
}
//}}}

//{{{void test_uint64_t_ll_remove(void)
void test_uint64_t_ll_remove(void)
{
    struct uint64_t_ll *ll = NULL;

    uint64_t_ll_append(&ll, 1);
    uint64_t_ll_remove(&ll, 1);

    TEST_ASSERT_EQUAL(NULL, ll);


    uint64_t_ll_append(&ll, 1);
    uint64_t_ll_append(&ll, 2);

    uint64_t_ll_remove(&ll, 1);
    TEST_ASSERT_EQUAL(2, ll->head->val);
    uint64_t_ll_remove(&ll, 1);
    TEST_ASSERT_EQUAL(2, ll->head->val);
    uint64_t_ll_remove(&ll, 2);
    TEST_ASSERT_EQUAL(NULL, ll);


    uint64_t_ll_append(&ll, 1);
    uint64_t_ll_append(&ll, 2);
    uint64_t_ll_append(&ll, 5);

    uint64_t_ll_remove(&ll, 2);
    TEST_ASSERT_EQUAL(1, ll->head->val);
    TEST_ASSERT_EQUAL(5, ll->head->next->val);
    uint64_t_ll_remove(&ll, 5);
    uint64_t_ll_remove(&ll, 1);
    TEST_ASSERT_EQUAL(NULL, ll);

    uint64_t_ll_free((void **)&ll);

    TEST_ASSERT_EQUAL(NULL, ll);
}
//}}}

//{{{void test_uint64_t_ll_contains(void)
void test_uint64_t_ll_contains(void)
{
    struct uint64_t_ll *ll = NULL;

    TEST_ASSERT_EQUAL(0, uint64_t_ll_contains(ll, 1));
    
    uint64_t_ll_append(&ll, 1);
    TEST_ASSERT_EQUAL(1, uint64_t_ll_contains(ll, 1));

    uint64_t_ll_append(&ll, 1);
    TEST_ASSERT_EQUAL(2, uint64_t_ll_contains(ll, 1));

    uint64_t_ll_remove(&ll, 1);
    TEST_ASSERT_EQUAL(0, uint64_t_ll_contains(ll, 1));


    uint64_t_ll_append(&ll, 1);
    uint64_t_ll_append(&ll, 2);
    uint64_t_ll_append(&ll, 5);

    TEST_ASSERT_EQUAL(1, uint64_t_ll_contains(ll, 5));
    TEST_ASSERT_EQUAL(1, uint64_t_ll_contains(ll, 1));
    TEST_ASSERT_EQUAL(1, uint64_t_ll_contains(ll, 2));

    TEST_ASSERT_EQUAL(0, uint64_t_ll_contains(ll, 7));

    uint64_t_ll_free((void **)&ll);

}
//}}}
