#define _GNU_SOURCE

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <sysexits.h>
#include <inttypes.h>
#include <htslib/kstring.h>
#include <err.h>
#include <sys/mman.h>
#include <unistd.h>

#include "offset_index.h"
#include "util.h"

char *OFFSET_INDEX_FILE_NAME = "offset_index.dat";

uint32_t offset_data_size = 0;

void (*offset_data_append_data)(uint8_t *dest, kstring_t *line) = NULL;

//{{{struct offset_index *offset_index_init(uint32_t init_size, char
struct offset_index *offset_index_init(uint32_t init_size, char *file_name)
{
    struct offset_index *oi = 
            (struct offset_index *) malloc(sizeof(struct offset_index));
    if (oi == NULL)
        err(1, "malloc error in offset_index_init()");

    oi->width = sizeof(struct file_id_offset_pair) + offset_data_size;
    oi->index = (struct file_id_offset_pairs *)
            malloc(sizeof(struct file_id_offset_pairs));
    if (oi->index == NULL)
        err(1, "malloc error in offset_index_init()");

    oi->index->num = 0;
    oi->index->size = init_size;
    /*
    oi->index->vals = (struct file_id_offset_pair *)
            calloc(oi->index->size, oi->width);
    */

    oi->file_name = NULL;
    if (file_name != NULL) {
        oi->file_name = strdup(file_name);
    }

    //oi->f = NULL;

    //oi->type = in_memory;
    
    oi->f = fopen(oi->file_name, "w+");
    if (oi->f == NULL)
        err(1, "Could not open %s.", oi->file_name);

    int ret = ftruncate(fileno(oi->f),
                        sizeof(uint64_t) + sizeof(uint32_t) +
                        (oi->index->size * oi->width));
    if ( ret != 0 )
        err(1, "Could not extend %s.", oi->file_name);

    oi->index->vals = (struct file_id_offset_pair *)
            mmap(0,
                 (oi->index->size * oi->width),
                 PROT_WRITE | PROT_READ,
                 MAP_SHARED,
                 fileno(oi->f),
                 0);

    if (oi->index->vals == MAP_FAILED)
        err(1, "Error mmapping file.");
 
    oi->type = on_disk;

    return oi;
}
//}}}

//{{{void offset_index_destroy(struct offset_index **oi);
void offset_index_destroy(struct offset_index **oi)
{

    if ((*oi)->file_name != NULL) {
        free((*oi)->file_name);
        (*oi)->file_name = NULL;
    }

    if ((*oi)->f != NULL) {
        munmap((*oi)->index->vals, 
               sizeof(uint64_t) + sizeof(uint32_t) + 
               (*oi)->index->size * (*oi)->width);
        fclose((*oi)->f);
    } else {
        free((*oi)->index->vals);
    }

    free((*oi)->index);

    free(*oi);
    *oi = NULL;
}
//}}}

//{{{uint32_t offset_index_add(struct offset_index *oi)
uint64_t offset_index_add(struct offset_index *oi,
                          long offset,
                          kstring_t *line,
                          uint32_t file_id)
{
    uint64_t id = oi->index->num;
    oi->index->num = oi->index->num + 1;

    if (oi->index->num == oi->index->size) {
        oi->index->size = oi->index->size * 2;

        int ret = munmap(oi->index->vals, 
                         sizeof(uint64_t) + sizeof(uint32_t) + 
                         oi->index->num * oi->width);
        if ( ret != 0 )
            err(1, "offset_index_add: Could not munmap %s.", oi->file_name);

        ret = ftruncate(fileno(oi->f), 
                        sizeof(uint64_t) + sizeof(uint32_t) +
                        oi->index->size * oi->width);
        if ( ret != 0 )
            err(1, "offset_index_add:Could not extend %s.", oi->file_name);

        oi->index->vals = (struct file_id_offset_pair *)
                mmap(0,
                     sizeof(uint64_t) + sizeof(uint32_t) +
                        oi->index->size * oi->width,
                     PROT_WRITE | PROT_READ,
                     MAP_SHARED,
                     fileno(oi->f),
                     0);

        if (oi->index->vals == MAP_FAILED)
            err(1, "Error mmapping file.");

        /*
        oi->index->vals = (struct file_id_offset_pair *)
                realloc(oi->index->vals, oi->index->size * oi->width);
        if (oi->index->vals == NULL)
            err(1, "realloc error in offset_index_add().\n");
        memset((uint8_t *)oi->index->vals + (oi->index->num * oi->width),
               0,
               oi->index->num * oi->width);
        */
    } 
    OFFSET_INDEX_PAIR(oi, id)->offset = offset;
    OFFSET_INDEX_PAIR(oi, id)->file_id = file_id;

    if (offset_data_append_data != NULL) 
        offset_data_append_data((uint8_t *)OFFSET_INDEX_DATA(oi, id), line);

   return id;
}
//}}}

//{{{void offset_index_store(struct offset_index *oi)
void offset_index_store(struct offset_index *oi)
{
    
    ( (uint64_t *) oi->index->vals)[0] = oi->index->num;
    ( (uint32_t *) oi->index->vals)[2] = oi->width;
    int ret = ftruncate(fileno(oi->f), 
                        sizeof(uint64_t) + sizeof(uint32_t) + 
                        oi->index->num * oi->width);
    if ( ret != 0 )
        err(1, "offset_index_store: Could not truncate.");
 
/*
    if (oi->file_name == NULL)
        errx(1,"No output file given for offset_index.");

    FILE *f = fopen(oi->file_name, "wb");
    if (f == NULL)
        err(1, "Could not open %s.", oi->file_name);

    if (fwrite(&(oi->index->num),
               sizeof(uint64_t),1, f) != 1)
        err(EX_IOERR, "Error writing offset_index num to '%s'.",
            oi->file_name);

    if (fwrite(&(oi->width),
               sizeof(uint32_t),1, f) != 1)
        err(EX_IOERR, "Error writing offset_index width to '%s'.",
            oi->file_name);

    if (fwrite(oi->index->vals, 
               oi->width,
               oi->index->num, f) != oi->index->num)
        err(EX_IOERR, "Error writing file_id offset pairs to '%s'.",
            oi->file_name);
    fclose(f);
*/
}
//}}}

//{{{struct offset_index *offset_index_load(char *file_name)
struct offset_index *offset_index_load(char *file_name)
{
    struct offset_index *oi = (struct offset_index *)
            malloc(sizeof(struct offset_index));
    if (oi == NULL)
        err(1, "malloc error in offset_index_load().\n");

    oi->file_name = strdup(file_name);

    oi->f = fopen(file_name, "rb");
    if (oi->f == NULL)
        err(1, "Could not open %s.", oi->file_name);

    oi->index = (struct file_id_offset_pairs *)
            malloc(sizeof(struct file_id_offset_pairs));
    if (oi->index == NULL)
        err(1, "malloc error in offset_index_load().\n");

    size_t fr = fread(&(oi->index->num), sizeof(uint64_t), 1, oi->f);
    check_file_read(oi->file_name, oi->f, 1, fr);

    fr = fread(&(oi->width), sizeof(uint32_t), 1, oi->f);
    check_file_read(oi->file_name, oi->f, 1, fr);

    oi->index->size = oi->index->num;

    off_t filesize = lseek(fileno(oi->f), 0, SEEK_END);
    rewind(oi->f);

    oi->index->vals = (struct file_id_offset_pair *)
            mmap(0,
                 filesize,
                 PROT_READ,
                 MAP_SHARED,
                 fileno(oi->f),
                 0);

    if (oi->index->vals == MAP_FAILED)
        err(1, "Error mmapping file.");
    
    oi->type = on_disk;

    return oi;
}
//}}}

//{{{struct file_id_offset_pair offset_index_get(struct offset_index *oi,
struct file_id_offset_pair offset_index_get(struct offset_index *oi,
                                            uint64_t id)
{
    return *(OFFSET_INDEX_PAIR(oi, id));
}
//}}}
