#ifndef __OFFSET_H__
#define __OFFSET_H__

#include <stdint.h>
#include <htslib/kstring.h>

char *OFFSET_INDEX_FILE_NAME;

uint32_t offset_data_size;

void (*offset_data_append_data)(uint8_t *dest, kstring_t *line);

struct file_id_offset_pair
{
    uint32_t file_id;
    long offset;
};
void *file_id_offset_pair_load(FILE *f, char *file_name);
void file_id_offset_pair_store(void *v, FILE *f, char *file_name);

struct file_id_offset_pairs
{
    uint64_t num,size;
    struct file_id_offset_pair *vals;
};

struct offset_index
{
    struct file_id_offset_pairs *index; //<! file_index/offse pair list
    char *file_name; //<! offset_index file name
    uint32_t width;
    FILE *f;
    enum {in_memory, on_disk} type;
};

struct offset_index *offset_index_init(uint32_t init_size, char *file_name);
void offset_index_destroy(struct offset_index **oi);
void offset_index_store(struct offset_index *oi);
uint64_t offset_index_add(struct offset_index *oi,
                          long offset,
                          kstring_t *line,
                          uint32_t file_id);
struct offset_index *offset_index_load(char *file_name);
struct file_id_offset_pair offset_index_get(struct offset_index *oi,
                                            uint64_t id);


#define OFFSET_INDEX_PAIR(offset_index, i) \
    ( (struct file_id_offset_pair *) (((uint8_t *)(offset_index->index->vals)) \
        + ((offset_index->type == in_memory) \
            ? 0 : sizeof(uint64_t) + sizeof(uint32_t)) \
        + ((uint64_t)i * (uint64_t)offset_index->width)) )

#define OFFSET_INDEX_DATA(offset_index, i) \
    ( (void *) (((uint8_t *)(offset_index->index->vals)) \
        + ((offset_index->type == in_memory) \
            ? 0 : sizeof(uint64_t) + sizeof(uint32_t)) \
        + ((uint64_t)i * (uint64_t)offset_index->width + sizeof(struct file_id_offset_pair))) )

#endif
