#ifndef __DATA_REG__
#define __DATA_REG__

#include <stdint.h>

struct data_reg
{
    void **data;
    uint32_t num_data, data_size;
};

struct data_reg *data_reg_init(uint32_t init_size);

void data_reg_destroy(struct data_reg **dr);

uint32_t data_reg_add(struct data_reg *dr,
                      void *data);

void *data_reg_get(struct data_reg *dr, uint32_t i);

//data_reg_save

//data_reg_open

//data_reg_sort

#endif
