#ifndef __BPT_H__
#define __BPT_H__

#include <stdint.h>
#include <stdio.h>
#include <stdbool.h>
#include "lists.h"
#include "cache.h"

uint32_t ORDER;

#define BPT_NODE_NUM_ELEMENTS (2*ORDER+10)
#define BPT_NODE_ELEMENT_SIZE sizeof(uint32_t)
#define BPT_NODE_SIZE (BPT_NODE_NUM_ELEMENTS*BPT_NODE_ELEMENT_SIZE)

#define BPT_ID(node)            ((node)->data[0])
#define BPT_PARENT(node)        ((node)->data[1])
#define BPT_IS_LEAF(node)       ((node)->data[2])
#define BPT_LEADING(node)       ((node)->data[3])
#define BPT_NEXT(node)          ((node)->data[4])
#define BPT_NUM_KEYS(node)      ((node)->data[5])
#define BPT_KEYS(node)          ((node)->data + 6)
#define BPT_POINTERS_BLOCK(node) (((node)->data + (6+ORDER+1))[0])
#define BPT_POINTERS(node)      ((node)->data + (6+ORDER+2))

struct ordered_set *id_to_offset_map;

struct bpt_node 
{
    // 0 parent
    // 1 is_leaf
    // 2 leading
    // 3 next
    // 4 num_keys
    // [5 ... 5+(ORDER+1)]
    // [(5+(ORDER+1)) ... ((5+ORDER+1)+1)+(ORDER+2))
    //
    // Total size is 5+(ORDER+1)+(ORDER+2) = 2*ORDER + 9
    
    uint32_t *data;
};

uint64_t bpt_node_serialize(void *deserialized, void **serialized);
uint64_t bpt_node_deserialize(void *serialized,
                              uint64_t serialized_size,
                              void **deserialized);
void bpt_node_free_mem(void **deserialized);

struct cache_handler bpt_node_cache_handler;

struct bpt_node *bpt_new_node(uint32_t domain);

uint32_t bpt_find_leaf(uint32_t domain, uint32_t curr, uint32_t key);

uint32_t bpt_place_new_key_value(uint32_t domain,
                                 uint32_t root_id,
                                 uint32_t *target_id,
                                 int *target_key_pos,
                                 uint32_t key,
                                 uint32_t value_id,
                                 struct cache_handler *handler);

uint32_t bpt_split_node(uint32_t domain,
                        uint32_t root_id,
                        uint32_t bpt_node_id,
                        uint32_t *lo_result_id,
                        uint32_t *hi_result_id,
                        int *lo_hi_split_point,
                        void (*repair)(uint32_t domain,
                                       struct bpt_node *,
                                       struct bpt_node *));

void (*bpt_node_repair)(uint32_t domain, struct bpt_node *, struct bpt_node *);

//uint64_t (*serialize_leading)(void *deserialized, uint8_t **serialized);
//uint64_t (*serialize_pointer)(void *deserialized, uint8_t **serialized);
//uint64_t serialize_uint32_t(void *deserialized, uint8_t **serialized);

void (*append)(uint32_t domain,
               uint32_t new_value_id,
               uint32_t existing_value_id,
               struct cache_handler *handler);

//struct bpt_node *bpt_to_node(void *n);

uint32_t bpt_insert(uint32_t domain,
                    uint32_t root_id,
                    uint32_t key,
                    uint32_t value_id,
                    struct cache_handler *handler,
                    uint32_t *leaf_id,
                    int *pos);

uint32_t bpt_insert_new_value(uint32_t domain,
                              uint32_t root_id,
                              uint32_t key,
                              void *value,
                              struct cache_handler *handler,
                              uint32_t *value_id,
                              uint32_t *leaf_id,
                              int *pos);

//void bpt_print_tree(struct bpt_node *curr, int level);
//void pbt_print_node(struct bpt_node *bpt_node);

int b_search(uint32_t key, const uint32_t *D, uint32_t D_size);

int bpt_find_insert_pos(struct bpt_node *leaf, uint32_t key);

/**
 * @brief Search for a posistion within a bplus tree
 *
 * @param domain cache domain, chrom id for GIGGLE
 * @param root_id id of the tree root
 * @param leaf_id set to the id of leaf that does (or should) contain the key
 * @param pos set to possition of the key if found
 * @param key search key
 *
 * @retval 0 if the key was not found, the id to the pointer if found
 *
 */
uint32_t bpt_find(uint32_t domain,
                  uint32_t root_id,
                  uint32_t *leaf_id,
                  int *pos,
                  uint32_t key);

void bpt_destroy_tree(struct bpt_node **root);

bool bpt_write_tree(uint32_t domain, uint32_t root_id);
void bpt_print_node(struct bpt_node *node);
#endif
