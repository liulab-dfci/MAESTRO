#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <inttypes.h>
#include <string.h>

#include "util.h"
#include "unity.h"
#include "lists.h"
#include "cache.h"
#include "giggle_index.h"

void setUp(void) { }
void tearDown(void) { }

//{{{void test_hash_list_init(void)
void test_hash_list_init(void)
{
    struct hash_list *hashl = hash_list_init(100000);

    uint64_t *ret = (uint64_t *)hash_list_get(hashl, 1);
    TEST_ASSERT_EQUAL(NULL, ret);

    uint64_t d = 2;
    uint32_t r = hash_list_add(hashl, 1, &d, sizeof(uint64_t));
    ret = (uint64_t *)hash_list_get(hashl, 1);
    TEST_ASSERT_EQUAL(2, *ret);

    ret = (uint64_t *)hash_list_get(hashl, 100);
    TEST_ASSERT_EQUAL(NULL, ret);
    d = 4;
    r = hash_list_add(hashl, 100, &d, sizeof(uint64_t));
    ret = (uint64_t *)hash_list_get(hashl, 100);
    TEST_ASSERT_EQUAL(4, *ret);

    ret = (uint64_t *)hash_list_get(hashl, 1000);
    TEST_ASSERT_EQUAL(NULL, ret);
    d = 6;
    r = hash_list_add(hashl, 1000, &d, sizeof(uint64_t));
    ret = (uint64_t *)hash_list_get(hashl, 1000);
    TEST_ASSERT_EQUAL(6, *ret);

    ret = (uint64_t *)hash_list_remove(hashl, 1);
    free(ret);
    ret = (uint64_t *)hash_list_remove(hashl, 100);
    free(ret);
    ret = (uint64_t *)hash_list_remove(hashl, 1000);
    free(ret);

    ret = (uint64_t *)hash_list_get(hashl, 1);
    TEST_ASSERT_EQUAL(NULL, ret);
    ret = (uint64_t *)hash_list_get(hashl, 100);
    TEST_ASSERT_EQUAL(NULL, ret);
    ret = (uint64_t *)hash_list_get(hashl, 1000);
    TEST_ASSERT_EQUAL(NULL, ret);

    uint32_t i;
    for (i = 1001; i < 100000; ++i) {
        d = i+1;
        r = hash_list_add(hashl, i, &d, sizeof(uint64_t));
    }

    hash_list_destroy(&hashl);
}
//}}}
