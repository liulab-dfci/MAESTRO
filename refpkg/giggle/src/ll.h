#ifndef __LL_H__
#define __LL_H__

#include "giggle_index.h"

struct uint64_t_ll_bpt_non_leading_data
{
    struct uint64_t_ll *SA, *SE;
};

struct uint64_t_ll_bpt_leading_data
{
    struct uint64_t_ll *B;
};

struct uint64_t_ll_node
{
    uint64_t val;
    struct uint64_t_ll_node *next;
};

struct uint64_t_ll
{
    uint64_t len;
    struct uint64_t_ll_node *head, *tail;
};


void uint64_t_ll_append(struct uint64_t_ll **ll, uint64_t val);
void uint64_t_ll_uniq_append(struct uint64_t_ll **ll, uint64_t val);
void uint64_t_ll_remove(struct uint64_t_ll **ll, uint64_t val);
uint64_t uint64_t_ll_contains(struct uint64_t_ll *ll, uint64_t val);
void uint64_t_ll_free(void **ll);

struct long_ll_node
{
    long val;
    struct long_ll_node *next;
};

struct long_ll
{
    uint64_t len;
    struct long_ll_node *head, *tail;
};

void long_ll_append(struct long_ll **ll, long val);
void long_ll_uniq_append(struct long_ll **ll, long val);
void long_ll_remove(struct long_ll **ll, long val);
uint64_t long_ll_contains(struct long_ll *ll, long val);
void long_ll_free(void **ll);

struct long_uint_ll_node
{
    long long_val;
    uint64_t uint_val;
    struct long_uint_ll_node *next;
};

struct long_uint_ll
{
    uint64_t len;
    struct long_uint_ll_node *head, *tail;
};

void long_uint_ll_append(struct long_uint_ll **ll, long long_val, uint64_t uint_val);
void long_uint_ll_uniq_append(struct long_uint_ll **ll, long val);
void long_uint_ll_remove(struct long_uint_ll **ll, long val);
uint64_t long_uint_ll_contains(struct long_uint_ll *ll, long val);
void long_uint_ll_free(void **ll);



// bpt_node_repair :: uint64_t_ll_leading_repair
void uint64_t_ll_leading_repair(uint32_t domain,
                                struct bpt_node *a,
                                struct bpt_node *b);

// giggle_data_handler :: int64_t_ll_giggle_data_handler
void *uint64_t_ll_new_non_leading(uint32_t domain);
void *uint64_t_ll_new_leading(uint32_t domain);
void uint64_t_ll_non_leading_SA_add_scalar(uint32_t domain,
                                           void *_nld,
                                           void *_id);
void uint64_t_ll_non_leading_SE_add_scalar(uint32_t domain,
                                           void *_nld,
                                           void *_id);
void uint64_t_ll_leading_B_add_scalar(uint32_t domain,
                                      void *_ld,
                                      void *_id);
void uint64_t_ll_leading_union_with_B(uint32_t domain,
                                      void **R,
                                      void *leading);
void uint64_t_ll_non_leading_union_with_SA_subtract_SE(uint32_t domain,
                                                       void **R,
                                                       void *d);
void uint64_t_ll_non_leading_union_with_SA(uint32_t domain, void **R, void *d);

struct giggle_def uint64_t_ll_giggle_data_handler;

// cache_handler :: uint64_t_ll_non_leading_cache_handler
uint64_t uint64_t_ll_non_leading_serialize(void *deserialized,
                                           void **serialized);
uint64_t uint64_t_ll_non_leading_deserialize(void *serialized,
                                             uint64_t serialized_size,
                                             void **deserialized);
void uint64_t_ll_non_leading_free(void **deserialized);
struct cache_handler uint64_t_ll_non_leading_cache_handler;

// cache_handler :: uint64_t_ll_leading_cache_handler
uint64_t uint64_t_ll_leading_serialize(void *deserialized,
                                       void **serialized);
uint64_t uint64_t_ll_leading_deserialize(void *serialized,
                                         uint64_t serialized_size,
                                         void **deserialized);
void uint64_t_ll_leading_free(void **deserialized);
struct cache_handler uint64_t_ll_leading_cache_handler;

void uint64_t_ll_giggle_set_data_handler();

struct cache_handler uint64_t_ll_wah_leading_cache_handler;
void uint64_t_ll_wah_giggle_set_data_handler();
struct cache_handler uint64_t_ll_wah_non_leading_cache_handler;
uint64_t uint64_t_ll_leading_serialize_to_wah(void *deserialized,
                                              void **serialized);

uint64_t uint64_t_ll_non_leading_serialize_to_wah(void *deserialized,
                                                  void **serialized);

void uint64_t_ll_map_intersection_to_offset_list(struct giggle_index *gi,
                                                 struct giggle_query_result *gqr,
                                                 void *_R);

#endif
