#define _GNU_SOURCE
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


//{{{void test_unordered_list_init(void)
void test_unordered_list_init(void)
{
    struct unordered_list *ul = unordered_list_init(10);

    TEST_ASSERT_EQUAL(ul->num, 0);
    TEST_ASSERT_EQUAL(ul->size, 10);

    unordered_list_destroy(&ul, NULL);

    TEST_ASSERT_EQUAL(ul, NULL);
}
//}}}

//{{{void test_unordered_list_add(void)
void test_unordered_list_add(void)
{
    struct unordered_list *ul = unordered_list_init(10);

    int i, V[20];

    for (i = 0; i < 20; ++i)
        V[i] = (i+1)*2; 

    for (i = 0; i < 9; ++i) {
        TEST_ASSERT_EQUAL(i, unordered_list_add(ul, (void *)(V + i)));
        TEST_ASSERT_EQUAL(V + i, ul->data[i]);
        TEST_ASSERT_EQUAL(i + 1, ul->num);
        TEST_ASSERT_EQUAL(10, ul->size);
    }
    for (; i < 19; ++i) {
        TEST_ASSERT_EQUAL(i, unordered_list_add(ul, (void *)(V + i)));
        TEST_ASSERT_EQUAL(V + i, ul->data[i]);
        TEST_ASSERT_EQUAL(i + 1, ul->num);
        TEST_ASSERT_EQUAL(20, ul->size);
    }

    unordered_list_destroy(&ul, NULL);
}
//}}}

//{{{void test_unordered_list_get(void)
void test_unordered_list_get(void)
{
    struct unordered_list *ul = unordered_list_init(10);

    int i, V[20];

    for (i = 0; i < 20; ++i)
        V[i] = (i+1)*2; 

    for (i = 0; i < 20; ++i) 
        TEST_ASSERT_EQUAL(i, unordered_list_add(ul, (void *)(V + i)));

    for (i = 0; i < 20; ++i)  {
        int *r = (int *)unordered_list_get(ul, i);
        TEST_ASSERT_EQUAL(V[i], (int)(*r));
    }

    void *r = unordered_list_get(ul, 5000);
    TEST_ASSERT_EQUAL(NULL, r);
    
    unordered_list_destroy(&ul, NULL);
}
//}}}

//{{{ void test_unordered_list_store_load_file_id_offset_pair(void)
void test_unordered_list_store_load_file_id_offset_pair(void)
{
    struct unordered_list *ul = unordered_list_init(10);


    uint32_t uints[20];
    long longs[20];

    uint32_t i;
    for (i = 0; i < 20; ++i) {
        uints[i] = rand();
        longs[i] = rand();

        struct file_id_offset_pair *p = 
                (struct file_id_offset_pair *)
                malloc(sizeof(struct file_id_offset_pair));
        p->file_id = uints[i];
        p->offset = longs[i];
        unordered_list_add(ul, p);
    }

    for (i = 0; i < 20; ++i) {
        struct file_id_offset_pair *p = unordered_list_get(ul, i);
        TEST_ASSERT_EQUAL(uints[i], p->file_id);
        TEST_ASSERT_EQUAL(longs[i], p->offset);
    }

    char *file_name = "test_unordered_list_store_load_file_id_offset_pair.dat";
    FILE *f = fopen(file_name, "wb");

    unordered_list_store(ul, f, file_name, file_id_offset_pair_store);

    fclose(f);

    f = fopen(file_name, "rb");

    struct unordered_list *ul_2 = unordered_list_load(f,
                                                    file_name,
                                                    file_id_offset_pair_load);

    for (i = 0; i < 20; ++i) {
        struct file_id_offset_pair *p = unordered_list_get(ul, i);
        struct file_id_offset_pair *p_2 = unordered_list_get(ul_2, i);
        TEST_ASSERT_EQUAL(p->file_id, p_2->file_id);
        TEST_ASSERT_EQUAL(p->offset, p_2->offset);
    }

    unordered_list_destroy(&ul, free_wrapper);
    unordered_list_destroy(&ul_2, free_wrapper);
    remove("test_unordered_list_store_load_file_id_offset_pair.dat");
}
//}}}

//{{{ void test_unordered_list_store_c_str(void)
void test_unordered_list_store_load_c_str(void)
{
    struct unordered_list *ul = unordered_list_init(5);

    char *A[10];
    asprintf(&(A[0]), "zero");
    asprintf(&(A[1]), "one");
    asprintf(&(A[2]), "two");
    asprintf(&(A[3]), "three");
    asprintf(&(A[4]), "four");
    asprintf(&(A[5]), "five");
    asprintf(&(A[6]), "six");
    asprintf(&(A[7]), "seven");
    asprintf(&(A[8]), "eight");
    asprintf(&(A[9]), "nine");

    uint32_t i;
    for (i = 0; i < 10; ++i)
        unordered_list_add(ul, A[i]);

    for (i = 0; i < 10; ++i) {
        char *s = unordered_list_get(ul, i);
        TEST_ASSERT_TRUE(strcmp(A[i], s) == 0)
    }

    char *file_name = "test_unordered_list_store_load_c_str.dat";
    FILE *f = fopen(file_name, "wb");

    unordered_list_store(ul, f, file_name, c_str_store);

    fclose(f);

    f = fopen(file_name, "rb");

    struct unordered_list *ul_2 = unordered_list_load(f,
                                                    file_name,
                                                    c_str_load);

    for (i = 0; i < 10; ++i) {
        char *s = unordered_list_get(ul, i);
        char *s_2 = unordered_list_get(ul_2, i);
        TEST_ASSERT_TRUE(strcmp(s, s_2) == 0)
    }

    unordered_list_destroy(&ul, free_wrapper);
    unordered_list_destroy(&ul_2, free_wrapper);
    remove(file_name);
}
//}}}

//{{{void test_ordered_set_add(void)
void test_ordered_set_add(void)
{
    struct ordered_set *os = ordered_set_init(3,
                                              str_uint_pair_sort_element_cmp,
                                              str_uint_pair_search_element_cmp,
                                              str_uint_pair_search_key_cmp);
    struct str_uint_pair *p = (struct str_uint_pair *)
            malloc(sizeof(struct str_uint_pair));
    p->uint = os->num;
    p->str = strdup("one");

    struct str_uint_pair *r = (struct str_uint_pair *)
            ordered_set_add(os, p);

    TEST_ASSERT_EQUAL(0, r->uint);


    p = (struct str_uint_pair *) malloc(sizeof(struct str_uint_pair));
    p->uint = os->num;
    p->str = strdup("two");

    r = (struct str_uint_pair *) ordered_set_add(os, p);
    TEST_ASSERT_EQUAL(1, r->uint);

    p = (struct str_uint_pair *) malloc(sizeof(struct str_uint_pair));
    p->uint = os->num;
    p->str = strdup("one");

    r = (struct str_uint_pair *) ordered_set_add(os, p);
    TEST_ASSERT_EQUAL(0, r->uint);
    str_uint_pair_free((void **)&p);

    p = (struct str_uint_pair *) malloc(sizeof(struct str_uint_pair));
    p->uint = os->num;
    p->str = strdup("three");

    r = (struct str_uint_pair *) ordered_set_add(os, p);
    TEST_ASSERT_EQUAL(2, r->uint);

    p = (struct str_uint_pair *) malloc(sizeof(struct str_uint_pair));
    p->uint = os->num;
    p->str = strdup("four");

    r = (struct str_uint_pair *) ordered_set_add(os, p);
    TEST_ASSERT_EQUAL(3, r->uint);

    ordered_set_destroy(&os, str_uint_pair_free);
}
//}}}

//{{{void test_ordered_set_get(void)
void test_ordered_set_get(void)
{
    struct ordered_set *os = ordered_set_init(3,
                                              str_uint_pair_sort_element_cmp,
                                              str_uint_pair_search_element_cmp,
                                              str_uint_pair_search_key_cmp);
    struct str_uint_pair *p = (struct str_uint_pair *)
            malloc(sizeof(struct str_uint_pair));
    p->uint = os->num;
    p->str = strdup("one");
    struct str_uint_pair *r = (struct str_uint_pair *)ordered_set_add(os, p);


    p = (struct str_uint_pair *) malloc(sizeof(struct str_uint_pair));
    p->uint = os->num;
    p->str = strdup("three");
    r = (struct str_uint_pair *)ordered_set_add(os, p);

    p = (struct str_uint_pair *) malloc(sizeof(struct str_uint_pair));
    p->uint = os->num;
    p->str = strdup("four");
    r = (struct str_uint_pair *)ordered_set_add(os, p);

    p = (struct str_uint_pair *) malloc(sizeof(struct str_uint_pair));
    p->uint = os->num;
    p->str = strdup("two");
    r = (struct str_uint_pair *)ordered_set_add(os, p);

    char *one = "one";
    r = (struct str_uint_pair *) ordered_set_get(os, one);
    TEST_ASSERT_EQUAL(0, r->uint);

    char *two = "two";
    r = (struct str_uint_pair *) ordered_set_get(os, two);
    TEST_ASSERT_EQUAL(3, r->uint);

    char *three = "three";
    r = (struct str_uint_pair *) ordered_set_get(os, three);
    TEST_ASSERT_EQUAL(1, r->uint);

    char *four = "four";
    r = (struct str_uint_pair *) ordered_set_get(os, four);
    TEST_ASSERT_EQUAL(2, r->uint);

    ordered_set_destroy(&os, str_uint_pair_free);
}
//}}}

//{{{void test_ordered_set_store_load(void)
void test_ordered_set_store_load(void)
{
    struct ordered_set *os = ordered_set_init(3,
                                              str_uint_pair_sort_element_cmp,
                                              str_uint_pair_search_element_cmp,
                                              str_uint_pair_search_key_cmp);
    struct str_uint_pair *p = (struct str_uint_pair *)
            malloc(sizeof(struct str_uint_pair));
    p->uint = os->num;
    p->str = strdup("one");
    struct str_uint_pair *r = (struct str_uint_pair *)ordered_set_add(os, p);


    p = (struct str_uint_pair *) malloc(sizeof(struct str_uint_pair));
    p->uint = os->num;
    p->str = strdup("three");
    r = (struct str_uint_pair *)ordered_set_add(os, p);

    p = (struct str_uint_pair *) malloc(sizeof(struct str_uint_pair));
    p->uint = os->num;
    p->str = strdup("four");
    r = (struct str_uint_pair *)ordered_set_add(os, p);

    p = (struct str_uint_pair *) malloc(sizeof(struct str_uint_pair));
    p->uint = os->num;
    p->str = strdup("two");
    r = (struct str_uint_pair *)ordered_set_add(os, p);

    char *file_name = "test_ordered_set_store_load.dat"; 
    FILE *f = fopen(file_name, "wb");
    ordered_set_store(os, f, file_name,  str_uint_pair_store);
    fclose(f);

    f = fopen(file_name, "rb");

    struct ordered_set *os_2 = ordered_set_load(
                f,
                file_name,
                str_uint_pair_load,
                str_uint_pair_sort_element_cmp,
                str_uint_pair_search_element_cmp,
                str_uint_pair_search_key_cmp);

    fclose(f);

    TEST_ASSERT_EQUAL(os->num, os_2->num);
    TEST_ASSERT_EQUAL(os->size, os_2->size);

    struct str_uint_pair *r_2;

    char *one = "one";
    r = (struct str_uint_pair *) ordered_set_get(os, one);
    r_2 = (struct str_uint_pair *) ordered_set_get(os_2, one);
    TEST_ASSERT_EQUAL(r->uint, r_2->uint);

    char *two = "two";
    r = (struct str_uint_pair *) ordered_set_get(os, two);
    r_2 = (struct str_uint_pair *) ordered_set_get(os_2, two);
    TEST_ASSERT_EQUAL(r->uint, r_2->uint);

    char *three = "three";
    r = (struct str_uint_pair *) ordered_set_get(os, three);
    r_2 = (struct str_uint_pair *) ordered_set_get(os_2, three);
    TEST_ASSERT_EQUAL(r->uint, r_2->uint);

    char *four = "four";
    r = (struct str_uint_pair *) ordered_set_get(os, four);
    r_2 = (struct str_uint_pair *) ordered_set_get(os_2, four);
    TEST_ASSERT_EQUAL(r->uint, r_2->uint);

    ordered_set_destroy(&os, str_uint_pair_free);
    ordered_set_destroy(&os_2, str_uint_pair_free);
    remove(file_name);
}
//}}}

//{{{void test_uint_pair(void)
void test_uint_pair(void)
{
    /*
    struct uint_pair
    {
        uint32_t first,second;
    };
    int uint_pair_sort_by_first_element_cmp(const void *a, const void *b);
    int uint_pair_search_by_first_element_cmp(const void *a, const void *b);
    int uint_pair_search_by_first_cmp(const void *a, const void *b);
    */


    struct ordered_set *os =
            ordered_set_init(3,
                             uint_pair_sort_by_first_element_cmp,
                             uint_pair_search_by_first_element_cmp,
                             uint_pair_search_by_first_key_cmp);

    struct uint_pair *p = (struct uint_pair *)
            malloc(sizeof(struct uint_pair));
    p->first = 10;
    p->second = os->num;
    struct uint_pair *r = ordered_set_add(os, p);

    p = (struct uint_pair *) malloc(sizeof(struct uint_pair));
    p->first = 50;
    p->second = os->num;
    r = ordered_set_add(os, p);

    p = (struct uint_pair *) malloc(sizeof(struct uint_pair));
    p->first = 1;
    p->second = os->num;
    r = ordered_set_add(os, p);

    p = (struct uint_pair *) malloc(sizeof(struct uint_pair));
    p->first = 100;
    p->second = os->num;
    r = ordered_set_add(os, p);


    uint32_t k = 100;
    r = ordered_set_get(os, &k);
    TEST_ASSERT_EQUAL(3, r->second);

    k = 1;
    r = ordered_set_get(os, &k);
    TEST_ASSERT_EQUAL(2, r->second);

    k = 50;
    r = ordered_set_get(os, &k);
    TEST_ASSERT_EQUAL(1, r->second);

    k = 10;
    r = ordered_set_get(os, &k);
    TEST_ASSERT_EQUAL(0, r->second);

    k = 70;
    r = ordered_set_get(os, &k);
    TEST_ASSERT_EQUAL(NULL, r);

    p = (struct uint_pair *) malloc(sizeof(struct uint_pair));
    p->first = 100;
    p->second = os->num;
    r = ordered_set_add(os, p);
    TEST_ASSERT_EQUAL(3, r->second);


    ordered_set_destroy(&os, free_wrapper);

    TEST_ASSERT_EQUAL(NULL, os);
    free(p);
}
//}}}

//{{{void test_fifo_q(void)
void test_fifo_q(void)
{
    int V[5] = {1,2,3,4,5};

    struct fifo_q *q = NULL;
    int *r = (int *)fifo_q_peek(q);
    TEST_ASSERT_EQUAL(NULL, r);

    fifo_q_push(&q, V);
    r = (int *)fifo_q_peek(q);
    TEST_ASSERT_EQUAL(V[0], *r);

    r = (int *)fifo_q_pop(&q);
    TEST_ASSERT_EQUAL(V[0], *r);
    TEST_ASSERT_EQUAL(NULL, q);
    r = (int *)fifo_q_peek(q);
    TEST_ASSERT_EQUAL(NULL, r);

    int v = 0;
    fifo_q_push(&q, V + v++);
    fifo_q_push(&q, V + v++);
    fifo_q_push(&q, V + v++);
    fifo_q_push(&q, V + v++);
    fifo_q_push(&q, V + v++);

    r = (int *)fifo_q_peek(q);
    TEST_ASSERT_EQUAL(V[0], *r);

    r = (int *)fifo_q_pop(&q);
    TEST_ASSERT_EQUAL(V[0], *r);
    r = (int *)fifo_q_pop(&q);
    TEST_ASSERT_EQUAL(V[1], *r);
    r = (int *)fifo_q_pop(&q);
    TEST_ASSERT_EQUAL(V[2], *r);
    r = (int *)fifo_q_pop(&q);
    TEST_ASSERT_EQUAL(V[3], *r);
    r = (int *)fifo_q_pop(&q);
    TEST_ASSERT_EQUAL(V[4], *r);
    TEST_ASSERT_EQUAL(NULL, q);
}
//}}}

//{{{ void test_pointer_uint_pair(void)
void test_pointer_uint_pair(void)
{

    int V[10] = {1,2,3,4};

    struct ordered_set *os =
            ordered_set_init(3,
                             pointer_uint_pair_sort_element_cmp,
                             pointer_uint_pair_search_element_cmp,
                             pointer_uint_pair_search_key_cmp);

    struct pointer_uint_pair *p = (struct pointer_uint_pair *)
            malloc(sizeof(struct pointer_uint_pair));
    p->pointer = V+3;
    p->uint = 3;

    struct pointer_uint_pair *r = (struct pointer_uint_pair *)
            ordered_set_add(os, p);

    p = (struct pointer_uint_pair *)
            malloc(sizeof(struct pointer_uint_pair));
    p->pointer = V+1;
    p->uint = 1;
    r = (struct pointer_uint_pair *) ordered_set_add(os, p);

    p = (struct pointer_uint_pair *)
            malloc(sizeof(struct pointer_uint_pair));
    p->pointer = V;
    p->uint = 0;
    r = (struct pointer_uint_pair *) ordered_set_add(os, p);

    r = (struct pointer_uint_pair *) ordered_set_get(os, V);
    TEST_ASSERT_EQUAL(0,r->uint);
    r = (struct pointer_uint_pair *) ordered_set_get(os, V+1);
    TEST_ASSERT_EQUAL(1,r->uint);
    r = (struct pointer_uint_pair *) ordered_set_get(os, V+3);
    TEST_ASSERT_EQUAL(3,r->uint);
    r = (struct pointer_uint_pair *) ordered_set_get(os, V+2);
    TEST_ASSERT_EQUAL(NULL,r);

    ordered_set_destroy(&os, free_wrapper);
}
//}}}

//{{{ void test_pointer_uint_pair(void)
void test_uint_offset_size_pair(void)
{
    struct ordered_set *os =
            ordered_set_init(3,
                             uint_offset_size_pair_sort_element_cmp,
                             uint_offset_size_pair_search_element_cmp,
                             uint_offset_size_pair_search_key_cmp);

    struct uint_offset_size_pair *p = (struct uint_offset_size_pair *)
            malloc(sizeof(struct uint_offset_size_pair));
    p->uint = 1;
    p->offset = 50;
    struct uint_offset_size_pair *r = (struct uint_offset_size_pair *)
            ordered_set_add(os, p);

    p = (struct uint_offset_size_pair *)
            malloc(sizeof(struct uint_offset_size_pair));
    p->uint = 5;
    p->offset = 200;
    r = (struct uint_offset_size_pair *) ordered_set_add(os, p);

    p = (struct uint_offset_size_pair *)
            malloc(sizeof(struct uint_offset_size_pair));
    p->uint = 2;
    p->offset = 100;
    r = (struct uint_offset_size_pair *) ordered_set_add(os, p);

    p = (struct uint_offset_size_pair *)
            malloc(sizeof(struct uint_offset_size_pair));
    p->uint = 4;
    p->offset = 175;
    r = (struct uint_offset_size_pair *) ordered_set_add(os, p);

    p = (struct uint_offset_size_pair *)
            malloc(sizeof(struct uint_offset_size_pair));
    p->uint = 3;
    p->offset = 150;
    r = (struct uint_offset_size_pair *) ordered_set_add(os, p);


    uint32_t s = 1;
    r = (struct uint_offset_size_pair *) ordered_set_get(os, &s);
    TEST_ASSERT_EQUAL(50,r->offset);

    s = 2;
    r = (struct uint_offset_size_pair *) ordered_set_get(os, &s);
    TEST_ASSERT_EQUAL(100,r->offset);

    s = 3;
    r = (struct uint_offset_size_pair *) ordered_set_get(os, &s);
    TEST_ASSERT_EQUAL(150,r->offset);

    s = 4;
    r = (struct uint_offset_size_pair *) ordered_set_get(os, &s);
    TEST_ASSERT_EQUAL(175,r->offset);

    s = 5;
    r = (struct uint_offset_size_pair *) ordered_set_get(os, &s);
    TEST_ASSERT_EQUAL(200,r->offset);

    s = 7;
    r = (struct uint_offset_size_pair *) ordered_set_get(os, &s);
    TEST_ASSERT_EQUAL(NULL, r);

    s = 0;
    r = (struct uint_offset_size_pair *) ordered_set_get(os, &s);
    TEST_ASSERT_EQUAL(NULL, r);

    ordered_set_destroy(&os, free_wrapper);
}
//}}}

//{{{void test_byte_array(void)
void test_byte_array(void)
{
    struct byte_array *ba = byte_array_init(3 * sizeof(uint32_t));
    TEST_ASSERT_EQUAL(0, ba->num);
    TEST_ASSERT_EQUAL(3 * sizeof(uint32_t), ba->size);

    uint64_t A[3] = {2,4,6};
    uint32_t B[3] = {1,2,3};

    byte_array_append(ba, A, 3*sizeof(uint64_t));
    TEST_ASSERT_EQUAL(3*sizeof(uint64_t), ba->num);

    byte_array_append(ba, B, 3*sizeof(uint32_t));
    TEST_ASSERT_EQUAL(3*sizeof(uint64_t) + 3*sizeof(uint32_t), ba->num);

    uint64_t *A_r = (uint64_t *)(ba->data);
    TEST_ASSERT_EQUAL(A[0], A_r[0]);
    TEST_ASSERT_EQUAL(A[1], A_r[1]);
    TEST_ASSERT_EQUAL(A[2], A_r[2]);

    uint32_t *B_r = (uint32_t *)(ba->data + 3*sizeof(uint64_t));
    TEST_ASSERT_EQUAL(B[0], B_r[0]);
    TEST_ASSERT_EQUAL(B[1], B_r[1]);
    TEST_ASSERT_EQUAL(B[2], B_r[2]);

    byte_array_append_zeros(ba, 5*sizeof(uint32_t));
    byte_array_append(ba, B, 3*sizeof(uint32_t));
    B_r = (uint32_t *)(ba->data + 3*sizeof(uint64_t));
    TEST_ASSERT_EQUAL(0,  B_r[3]);
    TEST_ASSERT_EQUAL(0,  B_r[4]);
    TEST_ASSERT_EQUAL(0,  B_r[5]);
    TEST_ASSERT_EQUAL(0,  B_r[6]);
    TEST_ASSERT_EQUAL(0,  B_r[7]);
    TEST_ASSERT_EQUAL(B[0], B_r[8]);
    TEST_ASSERT_EQUAL(B[1], B_r[9]);
    TEST_ASSERT_EQUAL(B[2], B_r[10]);

    byte_array_destory(&ba);
    TEST_ASSERT_EQUAL(NULL, ba);
}
//}}}

//{{{void test_indexed_list(void)
void test_indexed_list(void)
{
    struct offset_size_pair p, *r;

    uint64_t max_size = 50;

    struct indexed_list *il =
            indexed_list_init(5,
                              sizeof(struct offset_size_pair));
    uint64_t i;
    for (i = 0; i < 5; ++i) {
        p.offset = i;
        p.size = i*2;
        TEST_ASSERT_EQUAL(0, indexed_list_add(il, i, &p));
    }

    p.offset = i;
    p.size = i*2;
    TEST_ASSERT_EQUAL(1, indexed_list_add(il, i, &p));

    i += 1;

    for (; i < max_size; ++i) {
        p.offset = i;
        p.size = i*2;
        indexed_list_add(il, i, &p);
    }

   for (i = 0; i < max_size; ++i) {
        r = (struct offset_size_pair*)indexed_list_get(il, i);
        TEST_ASSERT_EQUAL(i, r->offset);
        TEST_ASSERT_EQUAL(i*2, r->size);
    }

    p.offset = 100;
    p.size = 200;
    TEST_ASSERT_EQUAL(1, indexed_list_add(il, 100, &p));

    r = (struct offset_size_pair*)indexed_list_get(il, 100);
    TEST_ASSERT_EQUAL(100, r->offset);
    TEST_ASSERT_EQUAL(200, r->size);

    for (i = 0; i < max_size; ++i) {
        r = (struct offset_size_pair*)indexed_list_get(il, i);
        TEST_ASSERT_EQUAL(i, r->offset);
        TEST_ASSERT_EQUAL(i*2, r->size);
    }

    FILE *f = fopen("test_indexed_list.tmp", "wb");
    indexed_list_write(il, f, "test_indexed_list.tmp");
    fclose(f);

    f = fopen("test_indexed_list.tmp", "rb");
    struct indexed_list *il_1 = indexed_list_load(f, "test_indexed_list.tmp");
    fclose(f);

    struct offset_size_pair *a, *b;
    for (i = 0; i < max_size; ++i) {
        a = (struct offset_size_pair*)indexed_list_get(il, i);
        b = (struct offset_size_pair*)indexed_list_get(il_1, i);
        TEST_ASSERT_EQUAL(a->offset, b->offset);
        TEST_ASSERT_EQUAL(a->size, b->size);
    }

    indexed_list_destroy(&il);
    indexed_list_destroy(&il_1);
    TEST_ASSERT_EQUAL(NULL, il);
    TEST_ASSERT_EQUAL(NULL, il_1);
    remove("test_indexed_list.tmp");
}
//}}}

//{{{ void test_indexed_list_nulls(void)
void test_indexed_list_nulls(void)
{
    struct offset_size_pair p, *r;

    uint32_t max_size = 50;

    struct indexed_list *il =
            indexed_list_init(max_size,
                              sizeof(struct offset_size_pair));

    uint32_t i;
    for (i = 0; i < max_size; ++i) {
        r = (struct offset_size_pair*)indexed_list_get(il, i);
        TEST_ASSERT_EQUAL(NULL, r);
    }

    indexed_list_destroy(&il);
    TEST_ASSERT_EQUAL(NULL, il);
}
//}}}

//////{{{void test_cc_hash(void)
////void test_cc_hash(void)
////{
////
////    struct cc_hash *hash = cc_hash_init(2500);
////
////    uint32_t size = 1000;
////    uint32_t K[size],V[size];
////
////    uint32_t i;
////
////    for (i = 0; i < size; ++i){
////        V[i] = i * 3;
////        K[i] = rand();
////        int r = cc_hash_add(hash, K[i], V+i);
////    }
////
////    for (i = 0; i < size; ++i){
////        uint32_t *r = (uint32_t *)cc_hash_get(hash, K[i]);
////        TEST_ASSERT_EQUAL(V[i], *r);
////    }
////
////    for (i = 0; i < size; ++i){
////        uint32_t *r = (uint32_t *)cc_hash_remove(hash, K[i]);
////        TEST_ASSERT_EQUAL(V[i], *r);
////        r = (uint32_t *)cc_hash_remove(hash, K[i]);
////        TEST_ASSERT_EQUAL(NULL, r);
////    }
////
////    cc_hash_destroy(&hash);
////}
//////}}}
//
////{{{void test_lru_cache(void)
////void test_lru_cache(void)
////{
////    struct lru_cache *lruc = lru_cache_init(5, NULL);
////
////    int V[10] = {2,4,6,8,10,12,14,16,18,20};
////
////    lru_cache_add(lruc, 1, V, NULL);
////    TEST_ASSERT_EQUAL(1, lruc->seen);
////
////    TEST_ASSERT_EQUAL(V, lruc->head->value);
////    TEST_ASSERT_EQUAL(V, lruc->tail->value);
////
////    int *r = (int *)lru_cache_get(lruc, 1);
////
////    TEST_ASSERT_EQUAL(V[0], *r);
////    TEST_ASSERT_EQUAL(NULL, lru_cache_get(lruc, 2));
////
////    lru_cache_add(lruc, 2, V+1, NULL);
////    TEST_ASSERT_EQUAL(2, lruc->seen);
////
////    TEST_ASSERT_EQUAL(V, lruc->head->value);
////    TEST_ASSERT_EQUAL(V+1, lruc->tail->value);
////
////    r = (int *)lru_cache_get(lruc, 1);
////
////    TEST_ASSERT_EQUAL(V+1, lruc->head->value);
////    TEST_ASSERT_EQUAL(V, lruc->tail->value);
////
////    lru_cache_add(lruc, 3, V+2, NULL);
////    TEST_ASSERT_EQUAL(3, lruc->seen);
////    lru_cache_add(lruc, 4, V+3, NULL);
////    TEST_ASSERT_EQUAL(4, lruc->seen);
////    lru_cache_add(lruc, 5, V+4, NULL);
////    TEST_ASSERT_EQUAL(5, lruc->seen);
////    lru_cache_add(lruc, 6, V+5, NULL);
////    TEST_ASSERT_EQUAL(6, lruc->seen);
////
////    TEST_ASSERT_EQUAL(NULL, lru_cache_get(lruc, 2));
////    r = (int *)lru_cache_get(lruc, 1);
////    TEST_ASSERT_EQUAL(V[0], *r);
////    r = (int *)lru_cache_get(lruc, 3);
////    TEST_ASSERT_EQUAL(V[2], *r);
////    r = (int *)lru_cache_get(lruc, 4);
////    TEST_ASSERT_EQUAL(V[3], *r);
////    r = (int *)lru_cache_get(lruc, 5);
////    TEST_ASSERT_EQUAL(V[4], *r);
////    r = (int *)lru_cache_get(lruc, 6);
////    TEST_ASSERT_EQUAL(V[5], *r);
////
////    TEST_ASSERT_EQUAL(6, lruc->seen);
////
////    lru_cache_destroy((void **)&lruc);
////}
/////}}}
//
////{{{void test_simple_cache(void)
//void test_simple_cache(void)
//{
//    struct simple_cache *sc = simple_cache_init(5, NULL);
//
//    int V[10] = {2,4,6,8,10,12,14,16,18,20};
//
//    simple_cache_add(sc, 1, V, NULL);
//    TEST_ASSERT_EQUAL(1, simple_cache_seen(sc));
//
//    int *r = (int *)simple_cache_get(sc, 1);
//    TEST_ASSERT_EQUAL(V[0], *r);
//
//    uint32_t i;
//    for (i = 1; i<10; ++i)
//        simple_cache_add(sc, i+1, V+i, NULL);
//    TEST_ASSERT_EQUAL(10, simple_cache_seen(sc));
//
//    for (i = 0; i<10; ++i) {
//        int *r = (int *)simple_cache_get(sc, i+1);
//        TEST_ASSERT_EQUAL(V[i], *r);
//    }
//
//    simple_cache_destroy((void **)&sc);
//}
/////}}}
//
////{{{ void test_simple_cache_with_disk(void)
//void test_simple_cache_with_disk(void)
//{
//    char *file_name = "test_simple_cache_with_disk.out";
//    FILE *f = NULL;
//    struct disk_store *ds = disk_store_init(10, &f, file_name);
//    uint32_t V[5] = {2,4,6,8,10};
//    uint32_t id = disk_store_append(ds, V, sizeof(uint32_t));
//    id = disk_store_append(ds, V+1, sizeof(uint32_t));
//    id = disk_store_append(ds, V+2, sizeof(uint32_t));
//    id = disk_store_append(ds, V+3, sizeof(uint32_t));
//    id = disk_store_append(ds, V+4, sizeof(uint32_t));
//
//    disk_store_destroy(&ds);
//
//    f = fopen(file_name, "r+");
//
//    struct simple_cache *sc = simple_cache_init(10, f);
//
//    TEST_ASSERT_EQUAL(5, sc->seen);
//    TEST_ASSERT_EQUAL(5, sc->num);
//    
//    uint32_t *r = (uint32_t *)simple_cache_get(sc, 0);
//    TEST_ASSERT_EQUAL(V[0], *r);
//    r = (uint32_t *)simple_cache_get(sc, 1);
//    TEST_ASSERT_EQUAL(V[1], *r);
//
//    r = (uint32_t *)simple_cache_get(sc, 1);
//    TEST_ASSERT_EQUAL(V[1], *r);
//
//    r = (uint32_t *)simple_cache_get(sc, 2);
//    TEST_ASSERT_EQUAL(V[2], *r);
//    r = (uint32_t *)simple_cache_get(sc, 3);
//    TEST_ASSERT_EQUAL(V[3], *r);
//    r = (uint32_t *)simple_cache_get(sc, 4);
//    TEST_ASSERT_EQUAL(V[4], *r);
//
//    simple_cache_destroy((void **)&sc);
//}
/////}}}

//{{{ void test_bit_map(void)
void test_bit_map(void)
{
    struct bit_map *bm = bit_map_init(100);

    TEST_ASSERT_EQUAL(0, bit_map_get(bm, 1000));

    bit_map_set(bm, 1000);

    TEST_ASSERT_EQUAL(1, bit_map_get(bm, 1000));


    uint32_t A[500];

    uint32_t i;

    A[0] = 1000;

    for (i = 1; i < 500; ++i) {
        A[i] = rand() % 10000;
        bit_map_set(bm, A[i]);
    }

    qsort(A, 500, sizeof(uint32_t), uint32_t_cmp);

    uint32_t j;
    for (i = 0; i < 499; ++i)  {
        TEST_ASSERT_EQUAL(1, bit_map_get(bm, A[i]));
        for (j = A[i]+1; j < A[i+1]; ++j)
            TEST_ASSERT_EQUAL(0, bit_map_get(bm, j));
    }

    char *file_name = "test_bit_map.out";

    FILE *f = fopen(file_name, "wb");

    bit_map_store(bm, f, file_name);

    fclose(f);

    f = fopen(file_name, "rb");
    struct bit_map *bm_r = bit_map_load(f, file_name);
    fclose(f);

    TEST_ASSERT_EQUAL(bm->num_bits, bm_r->num_bits);
    TEST_ASSERT_EQUAL(bm->num_ints, bm_r->num_ints);

    for (i = 0; i < 499; ++i)  {
        TEST_ASSERT_EQUAL(1, bit_map_get(bm_r, A[i]));
        for (j = A[i]+1; j < A[i+1]; ++j)
            TEST_ASSERT_EQUAL(0, bit_map_get(bm_r, j));
    }


    bit_map_destroy(&bm);
    bit_map_destroy(&bm_r);
    remove(file_name);
}
//}}}

//{{{void test_uint32_t_array(void)
void test_uint32_t_array(void)
{
    /*
struct uint32_t_array
{
    uint32_t num, size, *data;
};

struct uint32_t_array *uint32_t_array_init(uint32_t init_size);
void uint32_t_array_destroy(struct uint32_t_array **ua);
uint32_t uint32_t_array_add(struct uint32_t_array *ua, uint32_t val);
uint32_t *uint32_t_array_get(struct uint32_t_array *ua, uint32_t index);
*/

    struct uint32_t_array *ua = uint32_t_array_init(10);

    TEST_ASSERT_EQUAL(10, ua->size);
    TEST_ASSERT_EQUAL(0, ua->num);

    uint32_t ret = uint32_t_array_add(ua, 5);

    TEST_ASSERT_EQUAL(10, ua->size);
    TEST_ASSERT_EQUAL(0, ret);
    TEST_ASSERT_EQUAL(1, ua->num);
    TEST_ASSERT_EQUAL(5, ua->data[0]);

    uint32_t *ret_p = uint32_t_array_get(ua, 0);
    TEST_ASSERT_EQUAL(5, *ret_p);

    ret_p = uint32_t_array_get(ua, 100);
    TEST_ASSERT_EQUAL(NULL, ret_p);


    uint32_t i;
    for (i = 0; i < 10; ++i) {
        ret = uint32_t_array_add(ua, 5);
        TEST_ASSERT_EQUAL(i + 1, ret);
    }

    TEST_ASSERT_EQUAL(20, ua->size);
    TEST_ASSERT_EQUAL(11, ua->num);
   
    uint32_t_array_destroy(&ua);

    TEST_ASSERT_EQUAL(NULL, ua);
}
//}}}
