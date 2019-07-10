#include <stdlib.h>
#include <stdio.h>
#include <err.h>

#include "wah.h"
#include "timer.h"
#include "ll.h"

int main(int argc, char **argv)
{
    WAH_SIZE = 32;
    WAH_MAX_FILL_WORDS = (1<<(WAH_SIZE-1)) - 1;

    uint32_t size = atoi(argv[1]);

    uint32_t *W_1 = (uint32_t *)calloc(size, sizeof(uint32_t));

    uint32_t *W_2 = (uint32_t *)calloc(size, sizeof(uint32_t));

    uint32_t i;

    if (argv[2][0] == 'w') {
        uint8_t *w_2 = NULL;
        uint8_t *w_1 = NULL;

        start();
        for (i = 0; i < size; ++i) {
            W_1[i] = rand();
            wah_uniq_append(&w_1, W_1[i]);
        }
        stop();
        fprintf(stderr, "%lu\t", report());

        start();
        for (i = 0; i < size; ++i) {
            W_2[i] = rand();
            wah_uniq_append(&w_2, W_2[i]);
        }
        stop();
        fprintf(stderr, "%lu\t", report());


        uint8_t *r = NULL;
        uint32_t r_size = 0;

        start();
        uint32_t resize = wah_or(w_1, w_2, &r, &r_size);
        stop();
        fprintf(stderr, "%lu\n", report());
    } else if (argv[2][0] == 'l') {
     
        struct uint64_t_ll *l_1 = NULL;
        struct uint64_t_ll *l_2 = NULL;

        start();
        for (i = 0; i < size; ++i) {
            W_1[i] = rand();
            uint64_t_ll_append(&l_1, W_1[i]);
        }
        stop();
        fprintf(stderr, "%lu\t", report());

        start();
        for (i = 0; i < size; ++i) {
            W_2[i] = rand();
            uint64_t_ll_append(&l_2, W_2[i]);
        }
        stop();
        fprintf(stderr, "%lu\t", report());

        start();
        struct uint64_t_ll_node *curr = l_1->head;
        while (curr != NULL) {
            uint64_t_ll_uniq_append(&l_2, curr->val);
            curr = curr->next;
        }
        stop();
        fprintf(stderr, "%lu\n", report());
    }
}
