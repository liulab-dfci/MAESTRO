#ifndef __LISTS_H__
#define __LISTS_H__

#include <stdint.h>
#include <htslib/khash.h>
#include "disk_store.h"

// BITMAP
struct bit_map
{
    uint64_t num_bits;
    uint32_t num_ints;
    uint32_t *bm;
};

struct bit_map *bit_map_init(uint64_t bits);
struct bit_map *bit_map_load(FILE *f, char *file_name);
void bit_map_store(struct bit_map *b, FILE *f, char *file_name);
void bit_map_destroy(struct bit_map **b);
void bit_map_set(struct bit_map *b, uint64_t i);
uint32_t bit_map_get(struct bit_map *b, uint64_t q);

// INDEXED LIST
struct indexed_list
{
    char *data;
    uint64_t size;
    uint32_t element_size;
    struct bit_map *bm;
};
struct indexed_list *indexed_list_init(uint64_t init_size,
                                       uint32_t element_size);
void indexed_list_destroy(struct indexed_list **il);
uint32_t indexed_list_add(struct indexed_list *il,
                          uint64_t index,
                          void *data);
void *indexed_list_get(struct indexed_list *il, uint64_t index);
void indexed_list_write(struct indexed_list *il, FILE *f, char *file_name);
struct indexed_list *indexed_list_load(FILE *f, char *file_name);

struct offset_size_pair
{
    long offset;
    uint32_t size;
};

// HASH LIST
struct hash_list
{
    void *hash;
};

struct hash_list *hash_list_init();
uint32_t hash_list_add(struct hash_list *il,
                       uint32_t index,
                       void *data,
                       uint32_t data_size);
void *hash_list_get(struct hash_list *il,
                    uint32_t index);
void hash_list_write(struct indexed_list *il, FILE *f, char *file_name);
void *hash_list_remove(struct hash_list *il,
                       uint32_t index);
void hash_list_destroy(struct hash_list **hashl);
void hash_list_value_cache_handler_pair_destroy(struct hash_list **hashl);

// UNORDERD LIST
struct unordered_list
{
    void **data;
    uint32_t num, size;
};

struct unordered_list *unordered_list_init(uint32_t init_size);
void unordered_list_destroy(struct unordered_list **ul,
                            void (*free_data)(void **data));
uint32_t unordered_list_add(struct unordered_list *ul,
                            void *data);
void *unordered_list_get(struct unordered_list *ul, uint32_t i);
struct unordered_list *unordered_list_load(
                FILE *f,
                char *file_name,
                void *(*ul_load)(FILE *f, char *file_name));
void unordered_list_store(struct unordered_list *ul, 
                          FILE *f,
                          char *file_name,
                          void (*ul_store)(void *v, FILE *f, char *file_name));



// ORDERD SET
struct ordered_set
{
    void **data;
    uint32_t num, size;
    int (*sort_element_cmp)(const void *a, const void *b);
    int (*search_element_cmp)(const void *a, const void *b);
    int (*search_key_cmp)(const void *a, const void *b);
};

struct ordered_set 
        *ordered_set_init(
                uint32_t init_size,
                int (*sort_element_cmp)(const void *a, const void *b),
                int (*search_element_cmp)(const void *a, const void *b),
                int (*search_key_cmp)(const void *a, const void *b));

struct ordered_set 
        *ordered_set_load(
                FILE *f,
                char *file_name,
                void *(*os_load)(FILE *f, char *file_name),
                int (*sort_element_cmp)(const void *a, const void *b),
                int (*search_element_cmp)(const void *a, const void *b),
                int (*search_key_cmp)(const void *a, const void *b));

void ordered_set_destroy(struct ordered_set **os,
                         void (*free_data)(void **data));
void *ordered_set_add(struct ordered_set *os,
                       void *data);
void *ordered_set_get(struct ordered_set *os, void *key);

void ordered_set_store(struct ordered_set *os, 
                       FILE *f,
                       char *file_name,
                       void (*os_store)(void *v, FILE *f, char *file_name));
 

struct str_uint_pair
{
    char *str;
    uint32_t uint;
};

int str_uint_pair_sort_element_cmp(const void *a, const void *b);
int str_uint_pair_search_element_cmp(const void *a, const void *b);
int str_uint_pair_search_key_cmp(const void *a, const void *b);
void str_uint_pair_store(void *v, FILE *f, char *file_name);
void *str_uint_pair_load(FILE *f, char *file_name);
void str_uint_pair_free(void **v);

struct pointer_uint_pair
{
    void *pointer;
    uint32_t uint;
};
int pointer_uint_pair_sort_element_cmp(const void *a, const void *b);
int pointer_uint_pair_search_element_cmp(const void *a, const void *b);
int pointer_uint_pair_search_key_cmp(const void *a, const void *b);

struct uint_offset_size_pair
{
    uint32_t uint,size;
    long offset;
};
int uint_offset_size_pair_sort_element_cmp(const void *a, const void *b);
int uint_offset_size_pair_search_element_cmp(const void *a, const void *b);
int uint_offset_size_pair_search_key_cmp(const void *a, const void *b);

struct uint_pair
{
    uint32_t first,second;
};
int uint_pair_sort_by_first_element_cmp(const void *a, const void *b);
int uint_pair_search_by_first_element_cmp(const void *a, const void *b);
int uint_pair_search_by_first_key_cmp(const void *a, const void *b);


// FIFO Q
struct fifo_q_element
{
    void *val;
    struct fifo_q_element *next;
};

struct fifo_q
{
    struct fifo_q_element *head, *tail;
};

void fifo_q_push(struct fifo_q **q, void *val);
void *fifo_q_pop(struct fifo_q **q);
void *fifo_q_peek(struct fifo_q *q);

// BYTE ARRAY
struct byte_array
{
    char *data;
    uint32_t num, size;
};

struct byte_array *byte_array_init(uint32_t init_size);
void byte_array_destory(struct byte_array **ba);
void byte_array_append(struct byte_array *ba, void *data, uint32_t size);
void byte_array_append_zeros(struct byte_array *ba, uint32_t size);


// UINT32 ARRAY
struct uint32_t_array
{
    uint32_t num, size, *data;
};

struct uint32_t_array *uint32_t_array_init(uint32_t init_size);
void uint32_t_array_destroy(struct uint32_t_array **ua);
uint32_t uint32_t_array_add(struct uint32_t_array *ua, uint32_t val);
uint32_t uint32_t_array_set(struct uint32_t_array *ua,
                            uint32_t val,
                            uint32_t index);
uint32_t *uint32_t_array_get(struct uint32_t_array *ua, uint32_t index);

// UINT64 ARRAY
struct uint64_t_array
{
    uint64_t num, size, *data;
};

struct uint64_t_array *uint64_t_array_init(uint64_t init_size);
void uint64_t_array_destroy(struct uint64_t_array **ua);
uint64_t uint64_t_array_add(struct uint64_t_array *ua, uint64_t val);
uint64_t *uint64_t_array_get(struct uint64_t_array *ua, uint64_t index);

#endif
