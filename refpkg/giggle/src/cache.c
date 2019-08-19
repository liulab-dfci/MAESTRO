#define _GNU_SOURCE

#include <stdlib.h>
#include <stdio.h>
#include <stdint.h>
#include <err.h>
#include <sysexits.h>
#include <string.h>
#include <unistd.h>
#include <htslib/khash.h>

#include "cache.h"
#include "util.h"
#include "lists.h"

uint32_t CACHE_NAME_SPACE = 0;
void *_cache[10] =
    {NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL};
struct cache_def cache;


struct cache_handler uint32_t_cache_handler = {uint32_t_serialize, 
                                               uint32_t_deserialize,
                                               uint32_t_free_mem};
//{{{ uint32_t_cache_handler
//{{{uint64_t uint32_t_serialize(void *deserialized, void **serialized)
uint64_t uint32_t_serialize(void *deserialized, void **serialized)
{
    uint32_t *de = (uint32_t *)deserialized;

    //uint8_t *data = (uint8_t *)calloc(1, sizeof(uint32_t));
    uint8_t *data = (uint8_t *)calloc(sizeof(uint32_t), sizeof(uint8_t));
    if (data == NULL)
        err(1, "calloc error in uint32_t_serialize().");

    memcpy(data, de, sizeof(uint32_t));

    *serialized = (void *)data;
    return sizeof(uint32_t);
}
//}}}

//{{{ uint64_t uint32_t_deserialize(void *serialized,
uint64_t uint32_t_deserialize(void *serialized,
                             uint64_t serialized_size,
                             void **deserialized)
{
    uint8_t *data = (uint8_t *)serialized;

    uint32_t *u = (uint32_t *)malloc(sizeof(uint32_t));
    if (u == NULL)
        err(1, "malloc error in uint32_t_deserialize().");

    memcpy(u, data, sizeof(uint32_t));

    *deserialized = (void *)u;

    return sizeof(uint32_t);
}
//}}}

//{{{void uint32_t_free_mem(void **deserialized)
void uint32_t_free_mem(void **deserialized)
{
    uint32_t **de = (uint32_t **)deserialized;
    free(*de);
    *de = NULL;
}
//}}}
//}}}

//{{{ simple_cache
struct cache_def simple_cache_def = {
    simple_cache_init,
    simple_cache_seen,
    simple_cache_add,
    simple_cache_get,
    simple_cache_store,
    NULL,
    simple_cache_destroy
};

//{{{ void *simple_cache_init(uint32_t size,
void *simple_cache_init(uint32_t size,
                        uint32_t num_domains,
                        char **file_names) 
{
    struct simple_cache *sc = (struct simple_cache *)
            malloc(sizeof(struct simple_cache));
    if (sc == NULL)
        err(1, "malloc error in simple_cache_init().");

    _cache[CACHE_NAME_SPACE] = sc;
    cache = simple_cache_def;

    (void) pthread_mutex_init (&(sc->mutex), NULL);

    sc->num_domains = num_domains;

    sc->dirty = 0;

    char *max_bytes_env = getenv("GIGGLE_MAX_BYTES");
    if (max_bytes_env != NULL) {
        char *pend;
        sc->max_bytes = strtoull(max_bytes_env, &pend, 10);
    } else {
        sc->max_bytes = 1000000000;
    }

    sc->curr_bytes = 0;
    sc->lru_head = NULL;
    sc->lru_tail = NULL;

    //sc->ils = (struct indexed_list **)calloc(num_domains,
                                             //sizeof(struct indexed_list *));
    sc->hashls = (struct hash_list **)calloc(num_domains,
                                           sizeof(struct hash_list *));
    //if (sc->ils == NULL)
        //err(1, "calloc error in simple_cache_init().");
    if (sc->hashls == NULL)
        err(1, "calloc error in simple_cache_init().");

    sc->nums = (uint32_t *)calloc(num_domains, sizeof(uint32_t));
    if (sc->nums == NULL)
        err(1, "calloc error in simple_cache_init().");

    sc->seens = (uint32_t *)calloc(num_domains, sizeof(uint32_t));
    if (sc->seens == NULL)
        err(1, "calloc error in simple_cache_init().");

    sc->dss = NULL;
    sc->data_file_names = NULL;
    sc->index_file_names = NULL;

    uint32_t i;
    // opend attached files
    if (file_names != NULL) {
        sc->dss = (struct disk_store **)calloc(num_domains, 
                                               sizeof(struct disk_store *));
        if (sc->dss == NULL)
            err(1, "calloc error in simple_cache_init().");

        sc->index_file_names = (char **)calloc(num_domains, 
                                               sizeof(char *));
        if (sc->index_file_names == NULL)
            err(1, "calloc error in simple_cache_init().");

        sc->data_file_names = (char **)calloc(num_domains, 
                                              sizeof(char *));
        if (sc->data_file_names == NULL)
            err(1, "calloc error in simple_cache_init().");

        uint32_t ret;
        char *index_file_name, *data_file_name;
        for ( i = 0; i < num_domains; ++i) {
            ret = asprintf(&index_file_name, "%s.idx", file_names[i]);
            ret = asprintf(&data_file_name, "%s.dat", file_names[i]);
            sc->index_file_names[i] = index_file_name;
            sc->data_file_names[i] = data_file_name;
            // test to see if these files are in place
            if( (access(index_file_name, F_OK) != -1 ) && 
                (access(data_file_name, F_OK) != -1) )
                sc->dss[i] = disk_store_load(NULL,
                                             index_file_name,
                                             NULL,
                                             data_file_name);
            else 
                sc->dss[i] = NULL;
                /*
                sc->dss[i] = disk_store_init(size,
                                             NULL,
                                             index_file_name,
                                             NULL,
                                             data_file_name);
                */
        }
    }

    // set size
    sc->sizes = (uint32_t *)calloc(num_domains, sizeof(uint32_t));
    if (sc->sizes == NULL)
        err(1, "calloc error in simple_cache_init().");

    for ( i = 0; i < num_domains; ++i) {
        sc->sizes[i] = size;
        sc->nums[i] = 0;
        sc->seens[i] = 0;

        if ((sc->dss != NULL) && (sc->dss[i] != NULL)) {
            sc->nums[i] = sc->dss[i]->num;
            sc->seens[i] = sc->dss[i]->num;

            while (sc->sizes[i] < sc->dss[i]->num)
                sc->sizes[i] = sc->sizes[i] * 2;
        }

        //sc->ils[i] = indexed_list_init(sc->sizes[i],
                                       //sizeof(struct value_cache_handler_pair));
        uint64_t max_bytes = 10000000000; // 10GB
        sc->hashls[i] = hash_list_init();
    }

    return sc;
}
//}}}

//{{{uint32_t simple_cache_seen(void *_sc)
uint32_t simple_cache_seen(uint32_t domain)
{
    if (_cache[CACHE_NAME_SPACE] == NULL)
        errx(1, "Cache has not been initialized.");

    struct simple_cache *sc = (struct simple_cache *)_cache[CACHE_NAME_SPACE];
    return sc->seens[domain];
}
//}}}

//{{{void simple_cache_add(void *_sc,
void simple_cache_add(uint32_t domain,
                      uint32_t key,
                      void *value,
                      uint64_t value_size,
                      struct cache_handler *handler)
{
    if (_cache[CACHE_NAME_SPACE] == NULL)
        errx(1, "Cache has not been initialized.");
    struct simple_cache *sc = (struct simple_cache *)_cache[CACHE_NAME_SPACE];

    while (sc->max_bytes < sc->curr_bytes + value_size) {
        struct lru_ll_element *lru_curr = sc->lru_head;

        if (lru_curr == NULL)
            errx(1, "Not enough memory to load current element into cache.");

        struct value_cache_handler_pair *vh_curr = 
                (struct value_cache_handler_pair *)
                hash_list_remove(sc->hashls[lru_curr->domain], lru_curr->key);

#ifdef DEBUG_CACHE
    fprintf(stderr,
            "simple_cache_add: remove element from cache domain:%u key:%u\n",
            lru_curr->domain,
            lru_curr->key);
#endif

        if (vh_curr != NULL) {
            if ((vh_curr->handler != NULL) && 
                (vh_curr->handler->free_mem != NULL))
                vh_curr->handler->free_mem(&(vh_curr->value));

            free(vh_curr);
        }

        sc->curr_bytes -= lru_curr->size;

        sc->lru_head = lru_curr->next;
        free(lru_curr);
    }

    struct value_cache_handler_pair vh;
    vh.value = value;
    vh.handler = handler;
    vh.lru_node = (struct lru_ll_element *)
            malloc(sizeof(struct lru_ll_element));
    vh.lru_node->domain = domain;
    vh.lru_node->key = key;
    vh.lru_node->size = value_size;
    vh.lru_node->prev = NULL;
    vh.lru_node->next = NULL;

    if (sc->lru_head == NULL) {
        sc->lru_head = vh.lru_node;
        sc->lru_tail = vh.lru_node;
    } else {
        sc->lru_tail->next = vh.lru_node;
        vh.lru_node->prev = sc->lru_tail;
        sc->lru_tail = vh.lru_node;
    }

    //indexed_list_add(sc->ils[domain], key, &vh);
    hash_list_add(sc->hashls[domain],
                  key,
                  &vh,
                  sizeof(struct value_cache_handler_pair));
    sc->nums[domain] += 1;
    sc->seens[domain] += 1;
}
//}}}

//{{{void *simple_cache_get(void *_sc, uint32_t key)
void *simple_cache_get(uint32_t domain,
                       uint32_t key,
                       struct cache_handler *handler)
{
    if (_cache[CACHE_NAME_SPACE] == NULL)
        errx(1, "Cache has not been initialized.");
    struct simple_cache *sc = (struct simple_cache *)_cache[CACHE_NAME_SPACE];
    //struct value_cache_handler_pair *vh = indexed_list_get(sc->ils[domain],
                                                           //key);
    struct value_cache_handler_pair *vh = hash_list_get(sc->hashls[domain],
                                                       key);

    if (vh == NULL) {
#ifdef DEBUG_CACHE
    fprintf(stderr,
            "simple_cache_get: is not in cache domain:%u key:%u\n",
            domain,
            key);
#endif

        (void) pthread_mutex_lock (&(sc->mutex));
        //vh = indexed_list_get(sc->ils[domain],
                              //key);
        /*
        vh = hash_list_get(sc->hashls[domain],
                          key);
        if (vh != NULL) {
            (void) pthread_mutex_unlock (&(sc->mutex));
            return vh->value;
        }
        */

        if ((sc->dss != NULL) && (sc->dss[domain] != NULL)) {
            uint64_t size;
            void *raw = disk_store_get(sc->dss[domain], key, &size);
            if (raw == NULL)
                return NULL;

            sc->curr_bytes += size;

            void *v;
            uint64_t deserialized_size = handler->deserialize(raw,
                                                              size,
                                                              &v);

            simple_cache_add(domain, key, v, size, handler);
            free(raw);
            (void) pthread_mutex_unlock (&(sc->mutex));
            return v;
        } else {
            (void) pthread_mutex_unlock (&(sc->mutex));
            return NULL;
        }
#ifdef DEBUG_CACHE
    fprintf(stderr,
            "simple_cache_get: added to cache\n");
#endif

    } else { 
#ifdef DEBUG_CACHE
    fprintf(stderr,
            "simple_cache_get: is in cache domain:%u key:%u\n",
            domain,
            key);
#endif
        // move this to the tail
        if (sc->lru_tail != vh->lru_node) {
            // take out of the list
            if (sc->lru_head == vh->lru_node) {
                sc->lru_head = vh->lru_node->next;
            } else {
                vh->lru_node->next->prev = vh->lru_node->prev;
                vh->lru_node->prev->next = vh->lru_node->next;
            }

            vh->lru_node->prev = sc->lru_tail;
            sc->lru_tail->next = vh->lru_node;
            sc->lru_tail = vh->lru_node;
            sc->lru_tail->next = NULL;
        }

#ifdef DEBUG_CACHE
    fprintf(stderr,
            "simple_cache_get: LRU updatd\n");
#endif
        return vh->value;
    }
}
//}}}

//{{{void simple_cache_destroy()
void simple_cache_destroy()
{
    if (_cache[CACHE_NAME_SPACE] == NULL)
        errx(1, "Cache has not been initialized.");
    struct simple_cache *sc = (struct simple_cache *)_cache[CACHE_NAME_SPACE];

    (void) pthread_mutex_destroy (&(sc->mutex));

    uint32_t i,j;

    for (i = 0; i < sc->num_domains; ++i) {
#if 0
        for (j = 0; j <= sc->seens[i]; ++j) {
            //struct value_cache_handler_pair *vh =
                    //indexed_list_get(sc->ils[i], j);
            struct value_cache_handler_pair *vh =
                    lru_list_get(sc->lruls[i],
                                 j);
            if (vh != NULL) {
                if ((vh->handler != NULL) && (vh->handler->free_mem != NULL))
                    vh->handler->free_mem(&(vh->value));
            }
        }
#endif
        hash_list_value_cache_handler_pair_destroy(&(sc->hashls[i]));

        if (sc->data_file_names != NULL)
            free(sc->data_file_names[i]);

        if (sc->index_file_names != NULL)
            free(sc->index_file_names[i]);

        //indexed_list_destroy(&(sc->ils[i]));
    }

    if (sc->data_file_names != NULL)
        free(sc->data_file_names);
    if (sc->index_file_names != NULL)
        free(sc->index_file_names);

    //free(sc->ils);
    free(sc->hashls);
    free(sc->nums);
    free(sc->seens);
    if ( sc->dss != NULL) {
        for (i = 0; i < sc->num_domains; ++i)
            if ( sc->dss[i] != NULL) 
                disk_store_destroy(&(sc->dss[i]));
        free(sc->dss);
    }
    free(sc->sizes);


    struct lru_ll_element *lru_tmp, *lru_curr = sc->lru_head;
    while (lru_curr != NULL) {
        lru_tmp = lru_curr->next;
        free(lru_curr);
        lru_curr = lru_tmp;
    }

    free(sc);
    _cache[CACHE_NAME_SPACE] = NULL;
}
//}}}

//{{{void simple_cache_store(uint32_t domain,
void simple_cache_store(uint32_t domain,
                        uint32_t *disk_id_order)
{
    if (_cache[CACHE_NAME_SPACE] == NULL)
        errx(1, "Cache has not been initialized.");
    struct simple_cache *sc = (struct simple_cache *)_cache[CACHE_NAME_SPACE];

    if (sc->dss[domain] != NULL)
        errx(1, "Modifying and existing bpt is not currently supported.");

    /*
    fprintf(stderr, "%s %s\n", 
                                      sc->index_file_names[domain],
                                      sc->data_file_names[domain]);
    */
    sc->dss[domain] = disk_store_init(sc->seens[domain],
                                      NULL,
                                      sc->index_file_names[domain],
                                      NULL,
                                      sc->data_file_names[domain]);

    struct value_cache_handler_pair *vh;
    uint32_t mem_i, disk_i, ds_id;
    uint64_t serialized_size;
    void *v;
    for (disk_i = 0 ; disk_i < sc->seens[domain]; ++disk_i) {

        if (disk_id_order != NULL)
            mem_i = disk_id_order[disk_i];
        else
            mem_i = disk_i;

        //FIXME
        //vh = indexed_list_get(sc->ils[domain], mem_i);
        vh = hash_list_get(sc->hashls[domain], mem_i);
        if (vh == NULL)
            errx(1, "Value missing from cache.");

        if ((vh->handler == NULL) || (vh->handler->serialize == NULL))
            errx(1, "Cannot serialize given data without a valid handler.");

        serialized_size = vh->handler->serialize(vh->value,
                                                          &v);
        ds_id = disk_store_append(sc->dss[domain], v, serialized_size);

        if (disk_i != ds_id)
            errx(1, "Cache and disk are out of sync");
        free(v);
    }
}
//}}}
//}}}

//{{{void free_wrapper(void **v)
void free_wrapper(void **v)
{
    free(*v);
    *v = NULL;
}
//}}}

#if 0
//{{{ cc_hash

uint32_t hash_A(uint32_t x, uint32_t limit)
{
    x = x ^ (x>>4);
    x = (x^0xdeadbeef) + (x<<5);
    x = x ^ (x>>11);
    return x % limit;
}

uint32_t hash_B( uint32_t x, uint32_t limit)
{
    x = (x+0x7ed55d16) + (x<<12);
    x = (x^0xc761c23c) ^ (x>>19);
    x = (x+0x165667b1) + (x<<5);
    x = (x+0xd3a2646c) ^ (x<<9);
    x = (x+0xfd7046c5) + (x<<3);
    x = (x^0xb55a4f09) ^ (x>>16);
    return x % limit;
}

struct cc_hash *cc_hash_init(uint32_t size)
{
    struct cc_hash *hash = (struct cc_hash *)malloc(sizeof(struct cc_hash));
    hash->num = 0;
    hash->sizes = size / 2;
    hash->keys[0] = (uint32_t *) calloc(hash->sizes, sizeof(uint32_t));
    hash->keys[1] = (uint32_t *) calloc(hash->sizes, sizeof(uint32_t));
    hash->values[0] = (void **) calloc(hash->sizes, sizeof(void *));
    hash->values[1] = (void **) calloc(hash->sizes, sizeof(void *));
    
    hash->hashes[0] = hash_A;
    hash->hashes[1] = hash_B;

    return hash;
}

int cc_hash_add(struct cc_hash *hash, uint32_t key, void *value)
{
    uint32_t pos_0 = hash->hashes[0](key, hash->sizes);
    uint32_t pos_1 = hash->hashes[1](key, hash->sizes);

    if ( (hash->keys[0][pos_0] == key) || (hash->keys[1][pos_1] == key))
        return 1;

    uint32_t h_i = 0;
    uint32_t i, pos_i;

    for (i = 0; i < hash->sizes; ++i) {
        pos_i = hash->hashes[h_i](key, hash->sizes);
        if ( hash->values[h_i][pos_i] == NULL ) {
            hash->keys[h_i][pos_i] = key;
            hash->values[h_i][pos_i] = value;
            return 0;
        } else {
            uint32_t t_key = hash->keys[h_i][pos_i];
            void *t_value = hash->values[h_i][pos_i];
    
            hash->keys[h_i][pos_i] = key;
            hash->values[h_i][pos_i] = value;
    
            key = t_key;
            value = t_value;
            h_i = (h_i + 1) % 2;
        }
    }

    errx(1, "Could not place item\n");
}

void *cc_hash_get(struct cc_hash *hash, uint32_t key)
{
    uint32_t i, pos_i;
    for (i = 0; i < 2; ++i) {
        pos_i = hash->hashes[i](key, hash->sizes);
        if ((hash->keys[i][pos_i] == key) && (hash->values[i][pos_i] != NULL))
            return hash->values[i][pos_i];
    }

    return NULL;
}

void *cc_hash_remove(struct cc_hash *hash, uint32_t key)
{
    uint32_t i, pos_i;
    for (i = 0; i < 2; ++i) {
        pos_i = hash->hashes[i](key, hash->sizes);
        if ((hash->keys[i][pos_i] == key) && 
            (hash->values[i][pos_i] != NULL)) {
            void *r = hash->values[i][pos_i];
            hash->values[i][pos_i] = NULL;
            return r;
        }
    }

    return NULL;
}

void cc_hash_destroy(struct cc_hash **hash)
{
    free((*hash)->keys[0]);
    free((*hash)->keys[1]);
    free((*hash)->values[0]);
    free((*hash)->values[1]);
    free(*hash);
    *hash = NULL;
}
//}}}

//{{{ lru_cache
#if 0
struct cache_def lru_cache_def = {
    NULL,
    lru_cache_init,
    lru_cache_seen,
    lru_cache_add,
    lru_cache_get,
    lru_cache_remove,
    lru_cache_destroy
};

//{{{struct lru_cache *lru_cache_init(uint32_t init_size)
void *lru_cache_init(uint32_t init_size, FILE *fp)
{
    struct lru_cache *lruc = (struct lru_cache *)
            malloc(sizeof(struct lru_cache));
    lruc->size = init_size;
    lruc->num = 0;
    lruc->seen = 0;
    lruc->hash_table = cc_hash_init(init_size * 2.5);
    lruc->head = NULL;
    lruc->tail = NULL;
    return lruc;
}
//}}}

//{{{uint32_t lru_cache_seen(struct lru_cache *lruc)
//uint32_t lru_cache_seen(struct lru_cache *lruc)
uint32_t lru_cache_seen(void *_lruc)
{
    struct lru_cache *lruc = (struct lru_cache *)_lruc;
    return lruc->seen;
}
//}}}

//{{{void lru_cache_add(struct lru_cache *lruc, uint32_t key, void *value)
void lru_cache_add(void *_lruc,
                   uint32_t key,
                   void *value,
                   void (*free_value)(void **data))
{
    struct lru_cache *lruc = (struct lru_cache *)_lruc;
    if (cc_hash_get(lruc->hash_table, key) != NULL)
        return;

    if (lruc->num == lruc->size) {
        // the head node is the lru
        struct linked_list_node *to_rem_l = lruc->head;
        lruc->head = to_rem_l->next;
        lruc->head->prev = NULL;

        struct linked_list_node *to_rem_h = 
                (struct linked_list_node *)
                cc_hash_remove(lruc->hash_table, to_rem_l->key);

        if (to_rem_h != to_rem_l)
            errx(1, "Inconsistency in LRU cache");

        if (to_rem_l->free_value != NULL)
            to_rem_l->free_value(&(to_rem_l->value));

        free(to_rem_l);
        lruc->num -= 1;
    }

    struct linked_list_node *ll = (struct linked_list_node *)
            malloc(sizeof(struct linked_list_node));
    ll->key = key;
    ll->free_value = free_value;
    ll->prev = NULL;
    ll->next = NULL;
    ll->value = value;

    if (lruc->head == NULL) {
        lruc->head = ll;
    } else {
        ll->prev = lruc->tail;
        lruc->tail->next = ll;
    }

    lruc->tail = ll;

    int r = cc_hash_add(lruc->hash_table, key, ll);

    lruc->num += 1;
    lruc->seen += 1;
}
//}}}

//{{{void *lru_cache_get(struct lru_cache *lruc, uint32_t key)
void *lru_cache_get(void *_lruc, uint32_t key)
{
    struct lru_cache *lruc = (struct lru_cache *)_lruc;
    struct linked_list_node *ll =
        (struct linked_list_node *) cc_hash_get(lruc->hash_table, key);

    if (ll == NULL)
        return NULL;

    // move this to the tail
    if (lruc->tail != ll) {

        // take ll out of the list
        if (lruc->head == ll) 
            lruc->head = ll->next;
        else {
            ll->next->prev = ll->prev;
            ll->prev->next = ll->next;
        }

        ll->prev = lruc->tail;
        lruc->tail->next = ll;
        lruc->tail = ll;
        lruc->tail->next = NULL;
    }
        
    return ll->value;
}
//}}}

//{{{void lru_cache_remove(struct lru_cache *lruc, uint32_t key)
void lru_cache_remove(void *_lruc, uint32_t key)
{
    struct lru_cache *lruc = (struct lru_cache *)_lruc;
    struct linked_list_node *to_rem = 
                (struct linked_list_node *)
                cc_hash_remove(lruc->hash_table, key);

    if (to_rem == NULL)
        return;

    // Take it out of the list
    if (to_rem == lruc->head) {
        lruc->head = NULL;
        lruc->tail = NULL;
    } else if (to_rem == lruc->tail) {
        lruc->tail->prev->next = NULL;
        lruc->tail = lruc->tail->prev;
    } else {
        to_rem->prev->next = to_rem->next;
    }


    if (to_rem->free_value != NULL)
        to_rem->free_value(&(to_rem->value));

    free(to_rem);
    lruc->num -= 1;
}
//}}}

//{{{void lru_cache_destroy(struct lru_cache **lruc)
void lru_cache_destroy(void **_lruc)
{
    struct lru_cache **lruc = (struct lru_cache **)_lruc;
    cc_hash_destroy(&((*lruc)->hash_table));

    struct linked_list_node *curr, *tmp;
    curr = (*lruc)->head;

    while (curr != NULL) {
        tmp = curr->next;;
        /*
        if ( (*lruc)->free_value != NULL)
            (*lruc)->free_value(&(curr->value));
        */
        if (curr->free_value != NULL)
            curr->free_value(&(curr->value));
        free(curr);
        curr = tmp;
    }
    free(*lruc);
    *lruc = NULL;
}
//}}}
#endif
//}}}
#endif
