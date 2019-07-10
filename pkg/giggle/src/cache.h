#ifndef __CACHE_H__
#define __CACHE_H__

#include <stdint.h>
#include <pthread.h>
#include <htslib/khash.h>
#include "disk_store.h"
#include "lists.h"


struct cache_handler
{
    uint64_t (*serialize)(void *deserialized, void **serialized);
    uint64_t (*deserialize)(void *serialized,
                            uint64_t serialized_size,
                            void **deserialized);
    void (*free_mem)(void **deserialized);
};

struct cache_handler uint32_t_cache_handler;
uint64_t uint32_t_serialize(void *deserialized, void **serialized);
uint64_t uint32_t_deserialize(void *serialized,
                             uint64_t serialized_size,
                             void **deserialized);
void uint32_t_free_mem(void **deserialized);

struct cache_def {
    void *(*init)(uint32_t size,
                  uint32_t num_domains,
                  char **file_names);
    uint32_t (*seen)(uint32_t domain);
    void (*add)(uint32_t domain,
                uint32_t key,
                void *data,
                uint64_t data_size,
                struct cache_handler *handler);
    void *(*get)(uint32_t domain,
                 uint32_t key,
                 struct cache_handler *handler);
    void (*store)(uint32_t domain,
                  uint32_t *disk_id_order);
    void (*remove)(uint32_t domain, uint32_t key);
    void (*destroy)();
};

struct cache_def cache;

//void *_cache;
void *_cache[10];
uint32_t CACHE_NAME_SPACE;

struct value_cache_handler_pair
{
    void *value;
    struct cache_handler *handler;
    struct lru_ll_element *lru_node;
};

struct cache_def simple_cache_def;

struct lru_ll_element
{
    uint32_t domain, key;
    uint64_t size;
    struct lru_ll_element *prev, *next;
};

struct simple_cache
{
    struct indexed_list **ils;
    struct hash_list **hashls;
    uint32_t *sizes, *nums, *seens, num_domains;
    char **index_file_names, **data_file_names;
    struct disk_store **dss;
    pthread_mutex_t mutex;

    struct lru_ll_element *lru_head, *lru_tail;
    uint64_t max_bytes, curr_bytes;
    uint32_t dirty;
};

void *simple_cache_init(uint32_t size,
                        uint32_t num_domains,
                        char **file_names);
uint32_t simple_cache_seen(uint32_t domain);
void *simple_cache_get(uint32_t domain,
                       uint32_t key,
                       struct cache_handler *handler);
void simple_cache_remove(uint32_t key);
void simple_cache_add(uint32_t domain,
                      uint32_t key,
                      void *value,
                      uint64_t value_size,
                      struct cache_handler *handler);
void simple_cache_store(uint32_t domain,
                        uint32_t *disk_id_order);
void simple_cache_destroy();

void free_wrapper(void **v);


#if 0
struct cc_hash
{
    uint32_t num, sizes, *keys[2];
    void **values[2];
    uint32_t (*hashes[2])(uint32_t x, uint32_t limit);
};

uint32_t hash_A(uint32_t x, uint32_t limit);
uint32_t hash_B(uint32_t x, uint32_t limit);

struct cc_hash *cc_hash_init(uint32_t size);
int cc_hash_add(struct cc_hash *hash, uint32_t key, void *value);
void *cc_hash_get(struct cc_hash *hash, uint32_t key);
void *cc_hash_remove(struct cc_hash *hash, uint32_t key);
void cc_hash_destroy(struct cc_hash **hash);

struct linked_list_node
{
    void *value; 
    uint32_t key;
    struct linked_list_node *prev, *next;
    void (*free_value)(void **data);
};

struct lru_cache
{
    struct cc_hash *hash_table;
    struct linked_list_node *head, *tail;
    uint32_t size, num, seen;
};

void *lru_cache_init(uint32_t init_size, FILE *fp);
uint32_t lru_cache_seen(void *lruc); 
void *lru_cache_get(void *cache, uint32_t key);
void lru_cache_remove(void *cache, uint32_t key);
void lru_cache_add(void *cache,
                   uint32_t key,
                   void *value,
                   void (*free_value)(void **data));
void lru_cache_destroy(void **lruc);

struct cache_def lru_cache_def;
#endif
#endif
