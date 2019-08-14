#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <inttypes.h>
#include <string.h>

#include "unity.h"
#include "lists.h"
#include "cache.h"
#include "disk_store.h"

void setUp(void) { }
void tearDown(void) { }

//{{{ tripple
struct tripple {
    uint32_t a;
    char b;
    double c;
};

uint64_t tripple_serialize(void *deserialized, void **serialized)
{
    struct tripple *de = (struct tripple *)deserialized;

    uint8_t *data = (uint8_t *)calloc(1, sizeof(uint32_t) + 
                                         sizeof(char) + 
                                         sizeof(double));

    memcpy(data, &(de->a), sizeof(uint32_t));
    memcpy(data + sizeof(uint32_t), &(de->b), sizeof(char));
    memcpy(data + sizeof(uint32_t) + sizeof(char), &(de->c), sizeof(double));

    *serialized = (void *)data;
    return sizeof(uint32_t) + sizeof(char) + sizeof(double);
}

uint64_t tripple_deserialize(void *serialized,
                             uint64_t serialized_size,
                             void **deserialized)
{
    uint8_t *data = (uint8_t *)serialized;

    struct tripple *t = (struct tripple *)calloc(1, sizeof(struct tripple));

    memcpy(&(t->a), data, sizeof(uint32_t));
    memcpy(&(t->b), data + sizeof(uint32_t),sizeof(char));
    memcpy(&(t->c), data + sizeof(uint32_t) + sizeof(char), sizeof(double));

    *deserialized = (void *)t;

    return sizeof(struct tripple);
}

void tripple_free_mem(void **deserialized)
{
    struct tripple **de = (struct tripple **)deserialized;
    free(*de);
    *de = NULL;
}

struct cache_handler tripple_cache_handler = {tripple_serialize, 
                                              tripple_deserialize,
                                              tripple_free_mem};
//}}}

//{{{void test_simple_cache(void)
void test_simple_cache(void)
{
    struct simple_cache *sc = simple_cache_init(5, 3, NULL);

    struct tripple *t = (struct tripple *)calloc(1, sizeof(struct tripple));
    t->a = 10;
    t->b = 'x';
    t->c = 10.1;

    TEST_ASSERT_EQUAL(0, simple_cache_seen(0));
    TEST_ASSERT_EQUAL(0, simple_cache_seen(1));
    TEST_ASSERT_EQUAL(0, simple_cache_seen(2));

    simple_cache_add(1, 0, t, sizeof(struct tripple), &tripple_cache_handler);

    TEST_ASSERT_EQUAL(0, simple_cache_seen(0));
    TEST_ASSERT_EQUAL(1, simple_cache_seen(1));
    TEST_ASSERT_EQUAL(0, simple_cache_seen(2));

    t = (struct tripple *)calloc(1, sizeof(struct tripple));
    t->a = 11;
    t->b = 'a';
    t->c = 11.1;

    simple_cache_add(0, 0, t, sizeof(struct tripple), &tripple_cache_handler);

    t = (struct tripple *)calloc(1, sizeof(struct tripple));
    t->a = 12;
    t->b = 'b';
    t->c = 12.1;

    simple_cache_add(0, 1, t, sizeof(struct tripple), &tripple_cache_handler);

    t = (struct tripple *)calloc(1, sizeof(struct tripple));
    t->a = 13;
    t->b = 'c';
    t->c = 13.1;

    simple_cache_add(0, 2, t, sizeof(struct tripple), &tripple_cache_handler);


    TEST_ASSERT_EQUAL(3, simple_cache_seen(0));
    TEST_ASSERT_EQUAL(1, simple_cache_seen(1));
    TEST_ASSERT_EQUAL(0, simple_cache_seen(2));

    t = simple_cache_get(0, 0, &tripple_cache_handler);

    TEST_ASSERT_EQUAL(11, t->a);
    TEST_ASSERT_EQUAL('a', t->b);
    TEST_ASSERT_EQUAL_FLOAT(11.1 , t->c);

    t = simple_cache_get(0, 1, &tripple_cache_handler);
    TEST_ASSERT_EQUAL(12, t->a);
    TEST_ASSERT_EQUAL('b', t->b);
    TEST_ASSERT_EQUAL_FLOAT(12.1 , t->c);

    t = simple_cache_get(0, 2, &tripple_cache_handler);
    TEST_ASSERT_EQUAL(13, t->a);
    TEST_ASSERT_EQUAL('c', t->b);
    TEST_ASSERT_EQUAL_FLOAT(13.1 , t->c);

    t = simple_cache_get(0, 15, &tripple_cache_handler);
    TEST_ASSERT_EQUAL(NULL, t);

    simple_cache_destroy();

}
///}}}
//
////{{{ void test_simple_cache_with_disk(void)
//void test_simple_cache_with_disk(void)
//{
//
//    char *file_names[3] = {"tmp.one", "tmp.two", "tmp.three"};
//    struct simple_cache *sc = simple_cache_init(5, 3, file_names);
//
//    //{{{
//    struct tripple *t = (struct tripple *)calloc(1, sizeof(struct tripple));
//    t->a = 10;
//    t->b = 'x';
//    t->c = 10.1;
//
//    TEST_ASSERT_EQUAL(0, simple_cache_seen(0));
//    TEST_ASSERT_EQUAL(0, simple_cache_seen(1));
//    TEST_ASSERT_EQUAL(0, simple_cache_seen(2));
//
//    simple_cache_add(1, 0, t, sizeof(struct tripple), &tripple_cache_handler);
//
//    TEST_ASSERT_EQUAL(0, simple_cache_seen(0));
//    TEST_ASSERT_EQUAL(1, simple_cache_seen(1));
//    TEST_ASSERT_EQUAL(0, simple_cache_seen(2));
//
//    t = (struct tripple *)calloc(1, sizeof(struct tripple));
//    t->a = 11;
//    t->b = 'a';
//    t->c = 11.1;
//
//    simple_cache_add(0, 0, t, sizeof(struct tripple), &tripple_cache_handler);
//
//    t = (struct tripple *)calloc(1, sizeof(struct tripple));
//    t->a = 12;
//    t->b = 'b';
//    t->c = 12.1;
//
//    simple_cache_add(0, 1, t, sizeof(struct tripple),&tripple_cache_handler);
//
//    t = (struct tripple *)calloc(1, sizeof(struct tripple));
//    t->a = 13;
//    t->b = 'c';
//    t->c = 13.1;
//
//    simple_cache_add(0, 2, t, sizeof(struct tripple), &tripple_cache_handler);
//    //}}}
//
//    simple_cache_store(0, NULL);
//    simple_cache_store(1, NULL);
//    simple_cache_store(2, NULL);
//
//    simple_cache_destroy();
//
//    sc = simple_cache_init(5, 3, file_names);
//
//    TEST_ASSERT_EQUAL(3, simple_cache_seen(0));
//    TEST_ASSERT_EQUAL(1, simple_cache_seen(1));
//    TEST_ASSERT_EQUAL(0, simple_cache_seen(2));
//
//    t = simple_cache_get(0, 0, &tripple_cache_handler);
//    TEST_ASSERT_EQUAL(11, t->a);
//    TEST_ASSERT_EQUAL('a', t->b);
//    TEST_ASSERT_EQUAL_FLOAT(11.1 , t->c);
//
//    t = simple_cache_get(0, 1, &tripple_cache_handler);
//    TEST_ASSERT_EQUAL(12, t->a);
//    TEST_ASSERT_EQUAL('b', t->b);
//    TEST_ASSERT_EQUAL_FLOAT(12.1 , t->c);
//
//    t = simple_cache_get(0, 2, &tripple_cache_handler);
//    TEST_ASSERT_EQUAL(13, t->a);
//    TEST_ASSERT_EQUAL('c', t->b);
//    TEST_ASSERT_EQUAL_FLOAT(13.1 , t->c);
//    TEST_ASSERT_FALSE(13.1 == 9.9);
//
//    t = simple_cache_get(0, 15, &tripple_cache_handler);
//    TEST_ASSERT_EQUAL(NULL, t);
//
//    simple_cache_destroy();
//
//    remove("tmp.one.idx");
//    remove("tmp.two.idx");
//    remove("tmp.three.idx");
//
//    remove("tmp.one.dat");
//    remove("tmp.two.dat");
//    remove("tmp.three.dat");
//}
/////}}}
//
////{{{ void test_simple_cache_with_disk(void)
//void test_simple_cache_with_new_disk_order(void)
//{
//    char *file_names[3] = {"tmp.one", "tmp.two", "tmp.three"};
//    struct simple_cache *sc = simple_cache_init(5, 3, file_names);
//
//    //{{{
//    struct tripple *t = (struct tripple *)calloc(1, sizeof(struct tripple));
//    t->a = 10;
//    t->b = 'x';
//    t->c = 10.1;
//
//    TEST_ASSERT_EQUAL(0, simple_cache_seen(0));
//    TEST_ASSERT_EQUAL(0, simple_cache_seen(1));
//    TEST_ASSERT_EQUAL(0, simple_cache_seen(2));
//
//    simple_cache_add(1, 0, t, sizeof(struct tripple), &tripple_cache_handler);
//
//    TEST_ASSERT_EQUAL(0, simple_cache_seen(0));
//    TEST_ASSERT_EQUAL(1, simple_cache_seen(1));
//    TEST_ASSERT_EQUAL(0, simple_cache_seen(2));
//
//    t = (struct tripple *)calloc(1, sizeof(struct tripple));
//    t->a = 11;
//    t->b = 'a';
//    t->c = 11.1;
//
//    simple_cache_add(0, 0, t, sizeof(struct tripple), &tripple_cache_handler);
//
//    t = (struct tripple *)calloc(1, sizeof(struct tripple));
//    t->a = 12;
//    t->b = 'b';
//    t->c = 12.1;
//
//    simple_cache_add(0, 1, t, sizeof(struct tripple), &tripple_cache_handler);
//
//    t = (struct tripple *)calloc(1, sizeof(struct tripple));
//    t->a = 13;
//    t->b = 'c';
//    t->c = 13.1;
//
//    simple_cache_add(0, 2, t, sizeof(struct tripple), &tripple_cache_handler);
//    //}}}
//
//    uint32_t new_order[3] = {1,2,0};
//    simple_cache_store(0, new_order);
//    simple_cache_store(1, NULL);
//    simple_cache_store(2, NULL);
//
//    simple_cache_destroy();
//
//    sc = simple_cache_init(5, 3, file_names);
//
//    TEST_ASSERT_EQUAL(3, simple_cache_seen(0));
//    TEST_ASSERT_EQUAL(1, simple_cache_seen(1));
//    TEST_ASSERT_EQUAL(0, simple_cache_seen(2));
//
//    //t = simple_cache_get(sc, 0, 0, &tripple_cache_handler);
//    t = simple_cache_get(0, 2, &tripple_cache_handler);
//    TEST_ASSERT_EQUAL(11, t->a);
//    TEST_ASSERT_EQUAL('a', t->b);
//    TEST_ASSERT_EQUAL_FLOAT(11.1 , t->c);
//
//    //t = simple_cache_get(sc, 0, 1, &tripple_cache_handler);
//    t = simple_cache_get(0, 0, &tripple_cache_handler);
//    TEST_ASSERT_EQUAL(12, t->a);
//    TEST_ASSERT_EQUAL('b', t->b);
//    TEST_ASSERT_EQUAL_FLOAT(12.1 , t->c);
//
//    //t = simple_cache_get(sc, 0, 2, &tripple_cache_handler);
//    t = simple_cache_get(0, 1, &tripple_cache_handler);
//    TEST_ASSERT_EQUAL(13, t->a);
//    TEST_ASSERT_EQUAL('c', t->b);
//    TEST_ASSERT_EQUAL_FLOAT(13.1 , t->c);
//    TEST_ASSERT_FALSE(13.1 == 9.9);
//
//    t = simple_cache_get(0, 15, &tripple_cache_handler);
//    TEST_ASSERT_EQUAL(NULL, t);
//
//    simple_cache_destroy();
//
//    remove("tmp.one.idx");
//    remove("tmp.two.idx");
//    remove("tmp.three.idx");
//
//    remove("tmp.one.dat");
//    remove("tmp.two.dat");
//    remove("tmp.three.dat");
//}
/////}}}
//
////{{{ void test_simple_cache_with_disk(void)
//void test_simple_cache_as_cache_with_new_disk_order(void)
//{
//    char *file_names[3] = {"tmp.one", "tmp.two", "tmp.three"};
//    struct simple_cache *sc = simple_cache_init(5, 3, file_names);
//
//    //{{{
//    struct tripple *t = (struct tripple *)calloc(1, sizeof(struct tripple));
//    t->a = 10;
//    t->b = 'x';
//    t->c = 10.1;
//
//    TEST_ASSERT_EQUAL(0, cache.seen(0));
//    TEST_ASSERT_EQUAL(0, cache.seen(1));
//    TEST_ASSERT_EQUAL(0, cache.seen(2));
//
//    cache.add(1, 0, t, sizeof(struct tripple), &tripple_cache_handler);
//
//    TEST_ASSERT_EQUAL(0, cache.seen(0));
//    TEST_ASSERT_EQUAL(1, cache.seen(1));
//    TEST_ASSERT_EQUAL(0, cache.seen(2));
//
//    t = (struct tripple *)calloc(1, sizeof(struct tripple));
//    t->a = 11;
//    t->b = 'a';
//    t->c = 11.1;
//
//    cache.add(0, 0, t, sizeof(struct tripple), &tripple_cache_handler);
//
//    t = (struct tripple *)calloc(1, sizeof(struct tripple));
//    t->a = 12;
//    t->b = 'b';
//    t->c = 12.1;
//
//    cache.add(0, 1, t, sizeof(struct tripple),  &tripple_cache_handler);
//
//    t = (struct tripple *)calloc(1, sizeof(struct tripple));
//    t->a = 13;
//    t->b = 'c';
//    t->c = 13.1;
//
//    cache.add(0, 2, t, sizeof(struct tripple),  &tripple_cache_handler);
//    //}}}
//
//    uint32_t new_order[3] = {1,2,0};
//    cache.store(0, new_order);
//    cache.store(1, NULL);
//    cache.store(2, NULL);
//
//    cache.destroy();
//
//    sc = simple_cache_init(5, 3, file_names);
//
//    TEST_ASSERT_EQUAL(3, cache.seen(0));
//    TEST_ASSERT_EQUAL(1, cache.seen(1));
//    TEST_ASSERT_EQUAL(0, cache.seen(2));
//
//    //t = simple_cache_get(sc, 0, 0, &tripple_cache_handler);
//    t = cache.get(0, 2, &tripple_cache_handler);
//    TEST_ASSERT_EQUAL(11, t->a);
//    TEST_ASSERT_EQUAL('a', t->b);
//    TEST_ASSERT_EQUAL_FLOAT(11.1 , t->c);
//
//    //t = simple_cache_get(sc, 0, 1, &tripple_cache_handler);
//    t = cache.get(0, 0, &tripple_cache_handler);
//    TEST_ASSERT_EQUAL(12, t->a);
//    TEST_ASSERT_EQUAL('b', t->b);
//    TEST_ASSERT_EQUAL_FLOAT(12.1 , t->c);
//
//    //t = simple_cache_get(sc, 0, 2, &tripple_cache_handler);
//    t = cache.get(0, 1, &tripple_cache_handler);
//    TEST_ASSERT_EQUAL(13, t->a);
//    TEST_ASSERT_EQUAL('c', t->b);
//    TEST_ASSERT_EQUAL_FLOAT(13.1 , t->c);
//    TEST_ASSERT_FALSE(13.1 == 9.9);
//
//    t = cache.get(0, 15, &tripple_cache_handler);
//    TEST_ASSERT_EQUAL(NULL, t);
//
//    cache.destroy();
//
//    remove("tmp.one.idx");
//    remove("tmp.two.idx");
//    remove("tmp.three.idx");
//
//    remove("tmp.one.dat");
//    remove("tmp.two.dat");
//    remove("tmp.three.dat");
//}
/////}}}
