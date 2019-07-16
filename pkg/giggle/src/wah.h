#ifndef __WAH_H__
#define __WAH_H__

#include <stdint.h>
#include "bpt.h"
#include "giggle_index.h"

#define BYTE 8
#define WAH_LEN(W) ( ((uint32_t *)W)[0] )
#define WAH_VAL(W,S) ( ((W >> (S-1))&1) == 1 ?  0 : W)
#define WAH_NUM_WORDS(W,S) ( ((W >> (S-1))&1) == 1 ?  W & ~(1<< (S-1)) : 1)

uint32_t WAH_SIZE;
uint32_t WAH_MAX_FILL_WORDS;

//uint8_t *wah_init(uint32_t word_size,
                  //uint32_t val);

uint8_t *wah_copy(uint8_t *w);

uint8_t *wah_init(uint32_t val);
uint32_t wah_or(uint8_t *X, uint8_t *Y, uint8_t **R, uint32_t *R_size);
uint32_t wah_nand(uint8_t *X, uint8_t *Y, uint8_t **R, uint32_t *R_size);

uint32_t wah_get_ints_count(uint8_t *X);
uint32_t wah_get_ints(uint8_t *X, uint32_t **R);

void wah_uniq_append(uint8_t **w, uint32_t id);


struct wah_bpt_non_leading_data
{
    uint8_t *SA, *SE;
};

struct wah_bpt_leading_data
{
    uint8_t *B;
};

// bpt_node_repair :: wah_leading_repair
void wah_leading_repair(uint32_t domain,
                          struct bpt_node *a,
                          struct bpt_node *b);

// giggle_data_handler :: wah_giggle_data_handler
void *wah_new_non_leading(uint32_t domain);
void *wah_new_leading(uint32_t domain);
void wah_non_leading_SA_add_scalar(uint32_t domain,
                                     void *_nld,
                                     void *_id);
void wah_non_leading_SE_add_scalar(uint32_t domain,
                                     void *_nld,
                                     void *_id);
void wah_leading_B_add_scalar(uint32_t domain,
                                void *_ld,
                                void *_id);
void wah_leading_union_with_B(uint32_t domain,
                                void **R,
                                void *leading);
void wah_non_leading_union_with_SA_subtract_SE(uint32_t domain,
                                                 void **R,
                                                 void *d);
void wah_non_leading_union_with_SA(uint32_t domain, void **R, void *d);

struct giggle_def wah_giggle_data_handler;

// cache_handler :: wah_non_leading_cache_handler
uint64_t wah_non_leading_serialize(void *deserialized,
                                     void **serialized);
uint64_t wah_non_leading_deserialize(void *serialized,
                                       uint64_t serialized_size,
                                       void **deserialized);
void wah_non_leading_free(void **deserialized);
struct cache_handler wah_non_leading_cache_handler;

// cache_handler :: wah_leading_cache_handler
uint64_t wah_leading_serialize(void *deserialized,
                                 void **serialized);
uint64_t wah_leading_deserialize(void *serialized,
                                   uint64_t serialized_size,
                                   void **deserialized);
void wah_leading_free(void **deserialized);
struct cache_handler wah_leading_cache_handler;

void wah_giggle_set_data_handler();


void set_wah_i(uint8_t *W, void *v, uint32_t word_size, uint32_t i);
void get_wah_i(uint8_t *W, void *v, uint32_t word_size, uint32_t i);

uint8_t *uints_to_wah(uint32_t *D, uint32_t D_num);
#endif
