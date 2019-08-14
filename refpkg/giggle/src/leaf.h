#ifndef __LEAF_H__
#define __LEAF_H__

#include "bpt.h"

#define LEAF_DATA_LEADING_START(node) (0)
#define LEAF_DATA_LEADING_END(node) (node->num_leading)

uint64_t LEAF_POINTERS_SIZE;
uint64_t LEAF_NUMS_SIZE;
uint64_t LEAF_LEADING_STARTS_ENDS_SIZE;

struct leaf_data {
    uint64_t num_leading, num_starts, num_ends;
    uint32_t *starts_pointers, *ends_pointers;
    uint64_t *leading, *starts, *ends, *data;
};

struct leaf_data_result {
    uint64_t len;
    uint64_t *data;
    struct leaf_data_result *next;
};

void leaf_data_free_mem(void **deserialized);
uint64_t leaf_data_deserialize(void *serialized,
                               uint64_t serialized_size,
                               void **deserialized);
uint64_t leaf_data_serialize(void *deserialized, void **serialized);

uint32_t leaf_data_starts_start(struct leaf_data *ld,
                                struct bpt_node *ln,
                                int i);

uint32_t leaf_data_starts_end(struct leaf_data *ld,
                              struct bpt_node *ln,
                              int i);

uint32_t leaf_data_ends_start(struct leaf_data *ld,
                              struct bpt_node *ln,
                              int i);

uint32_t leaf_data_ends_end(struct leaf_data *ld,
                            struct bpt_node *ln,
                            int i);

void leaf_data_print(struct leaf_data *ld);

#endif
