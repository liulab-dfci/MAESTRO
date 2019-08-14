#define _GNU_SOURCE

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <math.h>
#include <inttypes.h>
#include <string.h>
#include <err.h>


#include "wah.h"
#include "util.h"
#include "ll.h"

#define MIN(X, Y) (((X) < (Y)) ? (X) : (Y))
#define MAX(X, Y) (((X) > (Y)) ? (X) : (Y))

uint32_t WAH_SIZE = 32;
uint32_t WAH_MAX_FILL_WORDS = 2147483647; //(1<<(WAH_SIZE-1)) - 1

//{{{void set_wah_i(uint8_t *W, void *v, uint32_t word_size, uint32_t i)
void set_wah_i(uint8_t *W, void *v, uint32_t word_size, uint32_t i)
{
    memcpy(W + sizeof(uint32_t) + i*(word_size/BYTE), v, word_size/BYTE);
}
//}}}

//{{{void get_wah_i(uint8_t *W, void *v, uint32_t word_size, uint32_t i)
void get_wah_i(uint8_t *W, void *v, uint32_t word_size, uint32_t i)
{
    memcpy(v,
           W + sizeof(uint32_t) + i*(word_size/BYTE),
           (word_size/BYTE));
}
//}}}

//{{{uint8_t *wah_init(uint32_t val)
uint8_t *wah_init(uint32_t val)
{
    uint32_t bits_per_word = WAH_SIZE - 1;
    uint32_t num_words = (val + bits_per_word - 1) / bits_per_word;
    // the max number of words 8-bit fill word and represent is 
    // 2**7 - 1 = 127
    // LEN, and WAH_LEN is the number of words, it is independent of word size
    uint32_t len = 1 + (num_words > 1 ? 
            (num_words + WAH_MAX_FILL_WORDS - 1)/WAH_MAX_FILL_WORDS : 0);
    uint8_t *w = (uint8_t *)malloc(sizeof(uint32_t) + 
                    (len * (WAH_SIZE/BYTE)  * sizeof(uint8_t)));
    WAH_LEN(w) = len;

    uint32_t v, i = 0;
    uint32_t saved_words;
    while (val > bits_per_word) {
        saved_words = MIN(num_words - 1, WAH_MAX_FILL_WORDS);
        //WAH_I(w, WAH_SIZE, i) = (1 << (bits_per_word-10)) | (saved_words);
        v = (1 << (bits_per_word)) + (saved_words);
        //fprintf(stderr, "%u\n", v);
        set_wah_i(w, &v, WAH_SIZE, i);
        val -= saved_words * bits_per_word;
        num_words -= saved_words;
        i+=1;
    }

    if (val > 0) {
        //WAH_I(w, WAH_SIZE, i) =  1 << ( bits_per_word - val);
        v = 1 << ( bits_per_word - val);
        //fprintf(stderr, "%u\n", v);
        set_wah_i(w, &v, WAH_SIZE, i);
    } else {
        //WAH_I(w, WAH_SIZE,i) =  0;
        v = 0;
        //fprintf(stderr, "%u\n", v);
        set_wah_i(w, &v, WAH_SIZE, i);
    }

    return w;
}
//}}}

//{{{uint8_t *wah_copy(uint8_t *w)
uint8_t *wah_copy(uint8_t *w)
{
    if (w == NULL)
        return NULL;

    if (WAH_LEN(w) == 0)
        return NULL;

    uint32_t R_size = sizeof(uint32_t) + 
            (WAH_LEN(w) * (WAH_SIZE/BYTE) * sizeof(uint8_t));
    uint8_t *R = (uint8_t *)malloc(R_size);
    memcpy(R, w, R_size);

    return R;
}
//}}}

//{{{ uint32_t wah_or(uint8_t *X, uint8_t *Y, uint8_t **R, uint32_t *R_size)
uint32_t wah_or(uint8_t *X, uint8_t *Y, uint8_t **R, uint32_t *R_size)
{
    uint32_t R_i = 0, X_i = 0, Y_i = 0;
    uint32_t x, y;
    //uint8_t x, y;
    uint32_t x_size, y_size, r_size, y_done = 0, x_done = 0;
    uint32_t X_len = WAH_LEN(X), Y_len = WAH_LEN(Y);
    uint32_t R_len = X_len + Y_len;
    uint32_t reset_R = 0;

    if (*R == NULL) {
        //fprintf(stderr, "reset_R A\n");
        *R_size = sizeof(uint32_t) + (R_len*(WAH_SIZE/BYTE)*sizeof(uint8_t));
        *R = (uint8_t *)malloc(*R_size);
        memset(*R, 0, *R_size);
        reset_R = 1;
    } else if (*R_size < sizeof(uint32_t) + 
            (R_len*(WAH_SIZE/BYTE)*sizeof(uint8_t))) {
        /*
        fprintf(stderr, "reset_R B\tR_size:%u\t%lu\n",
                *R_size,
                sizeof(uint32_t) + (R_len*(WAH_SIZE/BYTE)*sizeof(uint8_t)));
        */
        free(*R);
        *R_size = sizeof(uint32_t) + (R_len*(WAH_SIZE/BYTE)*sizeof(uint8_t));
        *R = (uint8_t *)malloc(*R_size);
        memset(*R, 0, *R_size);
        reset_R = 1;
    }

    get_wah_i(X, &x, WAH_SIZE, X_i);
    get_wah_i(Y, &y, WAH_SIZE, Y_i);

    x_size = WAH_NUM_WORDS(x, WAH_SIZE);
    y_size = WAH_NUM_WORDS(y, WAH_SIZE);

    uint32_t v;
    while (1) {
        r_size = MIN(x_size, y_size);

        if (r_size > 1)  {
            v = ((1<< (WAH_SIZE - 1)) + r_size);
        } else {
            v = WAH_VAL(x, WAH_SIZE) | WAH_VAL(y, WAH_SIZE);
        }

        // Grow R if we need to
        if (sizeof(uint32_t) + R_i*(WAH_SIZE/BYTE)*sizeof(uint8_t) == *R_size) {
            uint32_t old_len = R_len;
            reset_R = 1;
            R_len = R_len * 2;
            *R_size = sizeof(uint32_t) + 
                    (R_len*(WAH_SIZE/BYTE)*sizeof(uint8_t));
            *R = (uint8_t *) realloc(*R, *R_size);
            memset(*R + sizeof(uint32_t) + 
                    (old_len*(WAH_SIZE/BYTE)*sizeof(uint8_t)),
                   0,
                   old_len*(WAH_SIZE/BYTE)*sizeof(uint8_t) );
        }

        //WAH_I(*R, WAH_SIZE, R_i) = (uint8_t) v;
        set_wah_i(*R, &v, WAH_SIZE, R_i);
        R_i += 1;

        x_size -= r_size;
        y_size -= r_size;

        if ((x_size == 0) && (x_done == 0)) {
            X_i += 1;
            if (X_i == X_len) {
                x_done = 1;
                x = 0;
            } else {
                //x = WAH_I(X, 8, X_i);
                x = 0;
                get_wah_i(X, &x, WAH_SIZE, X_i);
                x_size = WAH_NUM_WORDS(x, WAH_SIZE);
            }
        }

        if ((y_size == 0) && (y_done == 0)) {
            Y_i += 1;
            if (Y_i == Y_len) {
                y_done = 1;
                y = 0;
            } else {
                //y = WAH_I(Y, WAH_SIZE, Y_i);
                y = 0;
                get_wah_i(Y, &y, WAH_SIZE, Y_i);
                y_size = WAH_NUM_WORDS(y, WAH_SIZE);
            }
        }

        if ((x_done == 1) && (y_done == 1))
            break;
        else if (x_done == 1)
            x_size = y_size;
        else if (y_done == 1)
            y_size = x_size;
    }

    R_len = R_i;
    WAH_LEN(*R) = R_len;
    if (reset_R == 1) {
        *R_size = sizeof(uint32_t) + (R_len*(WAH_SIZE/BYTE)*sizeof(uint8_t));
        *R = (uint8_t *)realloc(*R, *R_size);
    }

    return reset_R;
}
//}}}

//{{{ uint32_t wah_nand(uint8_t *X, uint8_t *Y, uint8_t **R, uint32_t *R_size)
uint32_t wah_nand(uint8_t *X, uint8_t *Y, uint8_t **R, uint32_t *R_size)
{
    uint32_t R_i = 0, X_i = 0, Y_i = 0;
    uint32_t x, y;
    uint32_t x_size, y_size, r_size, y_done = 0, x_done = 0;
    uint32_t X_len = WAH_LEN(X), Y_len = WAH_LEN(Y);
    uint32_t R_len = X_len + Y_len;
    uint32_t reset_R = 0;

    if (*R == NULL) {
        *R_size = sizeof(uint32_t) + (R_len*(WAH_SIZE/BYTE)*sizeof(uint8_t));
        *R = (uint8_t *)malloc(*R_size);
        memset(*R, 0, *R_size);
        reset_R = 1;
    } else if (*R_size < sizeof(uint32_t) + 
            (R_len*(WAH_SIZE/BYTE)*sizeof(uint8_t))) {
        free(*R);
        *R_size = sizeof(uint32_t) + (R_len*(WAH_SIZE/BYTE)*sizeof(uint8_t));
        *R = (uint8_t *)malloc(*R_size);
        memset(*R, 0, *R_size);
        reset_R = 1;
    }

    //x = WAH_I(X, WAH_SIZE, X_i);
    get_wah_i(X, &x, WAH_SIZE, X_i);
    //y = WAH_I(Y, WAH_SIZE, Y_i);
    get_wah_i(Y, &y, WAH_SIZE, Y_i);

    x_size = WAH_NUM_WORDS(x, WAH_SIZE);
    y_size = WAH_NUM_WORDS(y, WAH_SIZE);

    uint32_t v;
    while (1) {
        r_size = MIN(x_size, y_size);

        if (r_size > 1)  {
            v = ((1<< (WAH_SIZE - 1)) + r_size);
            //v = (uint8_t) ((1<< (WAH_SIZE - 1)) + r_size);
        } else {
            v = (WAH_VAL(x, WAH_SIZE) & ~(WAH_VAL(y, WAH_SIZE)));
            //v = (uint8_t) (WAH_VAL(x, WAH_SIZE) | WAH_VAL(y, WAH_SIZE));
        }

        // Grow R if we need to
        if (sizeof(uint32_t) + R_i*(WAH_SIZE/BYTE)*sizeof(uint8_t) == *R_size) {
            uint32_t old_len = R_len;
            reset_R = 1;
            R_len = R_len * 2;
            *R_size = sizeof(uint32_t) + 
                    (R_len*(WAH_SIZE/BYTE)*sizeof(uint8_t));
            *R = (uint8_t *) realloc(*R, *R_size);
            memset(*R + sizeof(uint32_t) + 
                    (old_len*(WAH_SIZE/BYTE)*sizeof(uint8_t)),
                   0,
                   old_len*(WAH_SIZE/BYTE)*sizeof(uint8_t) );
        }

        //WAH_I(*R, WAH_SIZE, R_i) = (uint8_t) v;
        set_wah_i(*R, &v, WAH_SIZE, R_i);
        R_i += 1;

        x_size -= r_size;
        y_size -= r_size;

        if ((x_size == 0) && (x_done == 0)) {
            X_i += 1;
            if (X_i == X_len) {
                x_done = 1;
                x = 0;
            } else {
                //x = WAH_I(X, 8, X_i);
                x = 0;
                get_wah_i(X, &x, WAH_SIZE, X_i);
                x_size = WAH_NUM_WORDS(x, WAH_SIZE);
            }
        }

        if ((y_size == 0) && (y_done == 0)) {
            Y_i += 1;
            if (Y_i == Y_len) {
                y_done = 1;
                y = 0;
            } else {
                //y = WAH_I(Y, WAH_SIZE, Y_i);
                y = 0;
                get_wah_i(Y, &y, WAH_SIZE, Y_i);
                y_size = WAH_NUM_WORDS(y, WAH_SIZE);
            }
        }

        if ((x_done == 1) && (y_done == 1))
            break;
        else if (x_done == 1)
            x_size = y_size;
        else if (y_done == 1)
            y_size = x_size;
    }

    R_len = R_i;
    WAH_LEN(*R) = R_len;
    if (reset_R == 1) {
        *R_size = sizeof(uint32_t) + (R_len*(WAH_SIZE/BYTE)*sizeof(uint8_t));
        *R = (uint8_t *)realloc(*R, *R_size);
    }

    return reset_R;
}
//}}}

//{{{ uint32_t wah_get_ints(uint8_t *X, uint32_t **R)
uint32_t wah_get_ints(uint8_t *X, uint32_t **R)
{
    /*
    if (X == NULL)
        return 0;
    */
    //uint8_t x;
    uint32_t x;
    uint32_t x_i_size, x_size = 0;
    uint32_t X_len = WAH_LEN(X);
    uint32_t R_len = 0;

    uint32_t i;
    for (i = 0; i < X_len; ++i) {
        //x = WAH_I(X, WAH_SIZE, i);
        x = 0;
        get_wah_i(X, &x, WAH_SIZE, i);
        x_i_size = WAH_NUM_WORDS(x, WAH_SIZE);
        if (x_i_size == 1)
            R_len +=  __builtin_popcount(x);

        x_size += x_i_size * (WAH_SIZE - 1);
    }

    //__builtin_clz(x) takes in a unsigned int, so on smaller
    //types it will count extra zeros, diff counts how many extra there are 
    uint32_t diff = ((sizeof(unsigned int)*BYTE)/WAH_SIZE - 1)*WAH_SIZE;
    uint32_t offset = 0;

    *R = (uint32_t*)calloc(R_len, sizeof(uint32_t));
    uint32_t R_i = 0;
    x_size = 0;
    for (i = 0; i < X_len; ++i) {
        x = 0;
        get_wah_i(X, &x, WAH_SIZE, i);
        x_i_size = WAH_NUM_WORDS(x, WAH_SIZE);
        if ( x_i_size == 1 ) {
            while (x != 0) {
                offset = __builtin_clz(x) - diff;
                (*R)[R_i] = offset + x_size;
                R_i += 1;
                x &= ~(1 << (WAH_SIZE-1-offset));
            }
        }
        x_size += x_i_size * (WAH_SIZE - 1);
    }
    return R_len;
}
//}}}

//{{{ uint32_t wah_get_ints_count(uint8_t *X)
uint32_t wah_get_ints_count(uint8_t *X)
{
    uint8_t x;
    uint32_t x_i_size;
    uint32_t X_len = WAH_LEN(X);
    uint32_t R_len = 0;

    uint32_t i;
    for (i = 0; i < X_len; ++i) {

        //x = WAH_I(X, WAH_SIZE, i);
        x = 0;
        get_wah_i(X, &x, WAH_SIZE, i);

        x_i_size = WAH_NUM_WORDS(x, WAH_SIZE);

        if (x_i_size == 1)
            R_len +=  __builtin_popcount(x);
    }

    return R_len;
}
//}}}

//{{{void wah_uniq_append(uint8_t **w, uint32_t id)
void wah_uniq_append(uint8_t **w, uint32_t id)
{
    uint8_t *w_id = wah_init(id);

    if (*w == NULL) {
        *w = w_id;
    } else {
        uint8_t *r = NULL;
        uint32_t r_size = 0;
        uint32_t resize = wah_or(*w, w_id, &r, &r_size);
        free(*w);
        free(w_id);
        *w = r;
    }
}
//}}}

// NOTE:  this function could easily be rewritten using the general funtions
//{{{ void wah_leading_repair(uint32_t domain, 
void wah_leading_repair(uint32_t domain, 
                          struct bpt_node *a,
                          struct bpt_node *b)
{
#if DEBUG
    fprintf(stderr, "START leading_repair\n");
#endif

    if ( (BPT_IS_LEAF(a) == 1) && (BPT_IS_LEAF(b) == 1) ) {
        // make a new leading value for the new right node
        struct wah_bpt_leading_data *d = 
            (struct wah_bpt_leading_data *)
            malloc(sizeof(struct wah_bpt_leading_data));

        d->B = NULL;

        // if the left node had leading data, grab all of it
        if (BPT_LEADING(a) != 0) {
            // get it from cache
            struct wah_bpt_leading_data *l =  
                    cache.get(domain,
                              BPT_LEADING(a) - 1,
                              &wah_leading_cache_handler);
#if DEBUG
            uint32_t *R = NULL, R_size;
            uint32_t j;
            fprintf(stderr, "LEADING\t");
            if (l->B != NULL) {
                R_size = wah_get_ints(l->B, &R);
                fprintf(stderr, "R_size:%u\t", R_size);
                for (j = 0; j < R_size; ++j)
                    fprintf(stderr, "R[%u]:%u\t", j, R[j]);
                free(R);
                R = NULL;
            }
            fprintf(stderr, "\n");
#endif
            d->B = wah_copy(l->B);
        }

        uint32_t i;
        for (i = 0 ; i < BPT_NUM_KEYS(a); ++i) {
            struct wah_bpt_non_leading_data *nl = 
                    cache.get(domain,
                              BPT_POINTERS(a)[i] - 1,
                              &wah_non_leading_cache_handler);

#if DEBUG
            uint32_t *R = NULL, R_size;
            uint32_t j;

            fprintf(stderr, "%u\tSA\t",BPT_KEYS(a)[i]);
            if (nl->SA != NULL) {
                R_size = wah_get_ints(nl->SA, &R);
                fprintf(stderr, "R_size:%u\t", R_size);
                for (j = 0; j < R_size; ++j)
                    fprintf(stderr, "R[%u]:%u\t", j, R[j]);
                free(R);
                R = NULL;
            }
            fprintf(stderr, "\n");
            fprintf(stderr, "SE\t");
            if (nl->SE != NULL) {
                R_size = wah_get_ints(nl->SE, &R);
                fprintf(stderr, "R_size:%u\t", R_size);
                for (j = 0; j < R_size; ++j)
                    fprintf(stderr, "R[%u]:%u\t", j, R[j]);
                free(R);
                R = NULL;
            }
            fprintf(stderr, "\n");
#endif

            wah_non_leading_union_with_SA_subtract_SE(domain,
                                                        (void **)(&d->B),
                                                        nl);


        }

        if (d->B != NULL) {
            uint32_t v_id = cache.seen(domain) + 1;
            cache.add(domain,
                      v_id - 1,
                      d,
                      sizeof(struct wah_bpt_leading_data),
                      &wah_leading_cache_handler);
            BPT_LEADING(b) = v_id;
        } else {
            free(d);
        }
    }

#if DEBUG
    fprintf(stderr, "END leading_repair\n");
#endif
}
//}}}

//{{{ giggle_data_handler :: wah_giggle_data_handler
void wah_giggle_set_data_handler()
{
    bpt_node_repair = wah_leading_repair;

    wah_giggle_data_handler.non_leading_cache_handler =
        wah_non_leading_cache_handler;
    wah_giggle_data_handler.leading_cache_handler = 
        wah_leading_cache_handler;
    wah_giggle_data_handler.new_non_leading = 
        wah_new_non_leading;
    wah_giggle_data_handler.new_leading = 
        wah_new_leading;
    wah_giggle_data_handler.non_leading_SA_add_scalar = 
        wah_non_leading_SA_add_scalar;
    wah_giggle_data_handler.non_leading_SE_add_scalar = 
        wah_non_leading_SE_add_scalar;
    wah_giggle_data_handler.leading_B_add_scalar = 
        wah_leading_B_add_scalar;
    wah_giggle_data_handler.leading_union_with_B = 
        wah_leading_union_with_B;
    wah_giggle_data_handler.non_leading_union_with_SA = 
        wah_non_leading_union_with_SA;
    wah_giggle_data_handler.non_leading_union_with_SA_subtract_SE = 
        wah_non_leading_union_with_SA_subtract_SE;
    wah_giggle_data_handler.write_tree = NULL;

   wah_giggle_data_handler.giggle_collect_intersection =
        giggle_collect_intersection_data_in_pointers;

    wah_giggle_data_handler.map_intersection_to_offset_list =
        uint64_t_ll_map_intersection_to_offset_list;

    giggle_data_handler = wah_giggle_data_handler;
}

//{{{void *wah_new_non_leading()
void *wah_new_non_leading(uint32_t domain)
{
    struct wah_bpt_non_leading_data *d = 
            (struct wah_bpt_non_leading_data *)
            malloc(sizeof( struct wah_bpt_non_leading_data));
    d->SA = NULL;
    d->SE = NULL;

    return (void *)d;
}
//}}}

//{{{void *wah_new_leading()
void *wah_new_leading(uint32_t domain)
{
    struct wah_bpt_leading_data *d = 
            (struct wah_bpt_leading_data *)
            malloc(sizeof( struct wah_bpt_leading_data));
    d->B = NULL;

    return (void *)d;
}
//}}}

//{{{void wah_non_leading_SA_add_scalar(void *d, void *v)
void wah_non_leading_SA_add_scalar(uint32_t domain,
                                     void *_nld,
                                     void *_id)
{
#if DEBUG
    fprintf(stderr, "uint64_t_ll_non_leading_SA_add_scalar\n");
#endif
    struct wah_bpt_non_leading_data *nld =
            (struct wah_bpt_non_leading_data *)_nld;
    uint32_t *id = (uint32_t *)_id;

#if DEBUG
    fprintf(stderr, "id:%u\n", *id);
#endif

    // We need to add 1 here b/c we cannot store zero
    wah_uniq_append(&(nld->SA), *id + 1);
}
//}}}

//{{{void wah_non_leading_SE_add_scalar(void *d, void *v)
void wah_non_leading_SE_add_scalar(uint32_t domain,
                                     void *_nld,
                                     void *_id)
{
#if DEBUG
    fprintf(stderr, "uint64_t_ll_non_leading_SA_add_scalar\n");
#endif
    struct wah_bpt_non_leading_data *nld =
            (struct wah_bpt_non_leading_data *)_nld;
    uint32_t *id = (uint32_t *)_id;

#if DEBUG
    fprintf(stderr, "id:%u\n", *id);
#endif

    // We need to add 1 here b/c we cannot store zero
    wah_uniq_append(&(nld->SE), *id + 1);
}
//}}}

//{{{void wah_leading_B_add_scalar(void *d, void *v)
void wah_leading_B_add_scalar(uint32_t domain,
                                void *_ld,
                                void *_id)
{
#if DEBUG
    fprintf(stderr, "uint64_t_ll_non_leading_SA_add_scalar\n");
#endif
    struct wah_bpt_leading_data *ld =
            (struct wah_bpt_leading_data *)_ld;
    uint32_t *id = (uint32_t *)_id;

#if DEBUG
    fprintf(stderr, "id:%u\n", *id);
#endif

    // We need to add 1 here b/c we cannot store zero
    wah_uniq_append(&(ld->B), *id + 1);
}
//}}}

//{{{void wah_leading_union_with_B(void **R, void *leading)
void wah_leading_union_with_B(uint32_t domain,
                                void **R,
                                void *leading)
{
    struct wah_bpt_leading_data *ld = 
            (struct wah_bpt_leading_data *)leading;

    if ((ld != NULL) && (ld->B != NULL)) {
        uint8_t **w = (uint8_t **)R;

        if (*w == NULL)
            *w = wah_init(0);

        uint8_t *r = NULL;
        uint32_t r_size = 0;
        uint32_t resize = wah_or(*w, ld->B, &r, &r_size);

        free(*w);
        *w = r;
    }
}
//}}}

//{{{void wah_non_leading_union_with_SA(void **R, void *d)
void wah_non_leading_union_with_SA(uint32_t domain, void **R, void *d)
{
    struct wah_bpt_non_leading_data *nld = 
            (struct wah_bpt_non_leading_data *) d;
    if (nld != NULL) {
        if ((nld->SA != NULL)) {
            uint8_t **w = (uint8_t **)R;

            if (*w == NULL)
                *w = wah_init(0);

            uint8_t *r = NULL;
            uint32_t r_size = 0;
            uint32_t resize = wah_or(*w, nld->SA, &r, &r_size);

            free(*w);
            *w = r;
        }
    }
}
//}}}

//{{{void wah_non_leading_union_with_SA_subtract_SE(uint32_t domain,
void wah_non_leading_union_with_SA_subtract_SE(uint32_t domain,
                                                 void **R,
                                                 void *d)
{
    struct wah_bpt_non_leading_data *nld = 
            (struct wah_bpt_non_leading_data *) d;

    if (nld != NULL) {
        if ((nld->SA != NULL)) {
            uint8_t **w = (uint8_t **)R;

            if (*w == NULL)
                *w = wah_init(0);

            uint8_t *r = NULL;
            uint32_t r_size = 0;
            uint32_t resize = wah_or(*w, nld->SA, &r, &r_size);

            free(*w);
            *w = r;
        }
        if ((nld->SE != NULL)) {
            uint8_t **w = (uint8_t **)R;

            if (*w == NULL)
                *w = wah_init(0);

            uint8_t *r = NULL;
            uint32_t r_size = 0;
            uint32_t resize = wah_nand(*w, nld->SE, &r, &r_size);

            free(*w);
            *w = r;
        }
    }
}
//}}}
//}}}

//{{{ wah_non_leading_cache_handler
struct cache_handler wah_non_leading_cache_handler = {
        wah_non_leading_serialize,
        wah_non_leading_deserialize,
        wah_non_leading_free
};

//{{{uint64_t wah_non_leading_serialize(void *deserialized,
uint64_t wah_non_leading_serialize(void *deserialized,
                                     void **serialized)
{
    if (deserialized == NULL) {
        *serialized = NULL;
        return 0;
    }

    struct wah_bpt_non_leading_data *d =  
            (struct wah_bpt_non_leading_data *)deserialized;

    uint32_t SA_len = 0, SE_len = 0, serialized_len;

    if (d->SA != NULL)
        SA_len = sizeof(uint32_t) + 
            WAH_LEN(d->SA)*(WAH_SIZE/BYTE)*sizeof(uint8_t);

    if (d->SE != NULL)
        SE_len = sizeof(uint32_t) + 
            WAH_LEN(d->SE)*(WAH_SIZE/BYTE)*sizeof(uint8_t);

    serialized_len = 2*sizeof(uint32_t) + SA_len + SE_len;

    uint8_t *data = (uint8_t *)malloc(serialized_len);

    uint32_t *data_u = (uint32_t *)data;
    data_u[0] = SA_len;
    data_u[1] = SE_len;

    uint32_t data_i = 2*sizeof(uint32_t);


    if (d->SA != NULL) 
        memcpy(data + data_i, d->SA, SA_len);

    data_i += SA_len;


    if (d->SE != NULL) 
        memcpy(data + data_i, d->SE, SE_len);

    data_i += SE_len;

    if (data_i != serialized_len)
        errx(1,
             "Issue with wah_non_leading_serlize lengths. "
             "Expected:%u observed:%u.",
             serialized_len,
             data_i);

    *serialized = data;

    return serialized_len; 
}
//}}}

//{{{uint64_t wah_non_leading_deserialize(void *serialized,
uint64_t wah_non_leading_deserialize(void *serialized,
                                       uint64_t serialized_size,
                                       void **deserialized)
{
    if ((serialized_size == 0) || (serialized == NULL)) {
        *deserialized = NULL;
        return 0;
    }

    if (serialized_size < sizeof(uint32_t)*2)
        errx(1,
             "Malformed wah_non_leading serialized value. "
             "Too short");

    uint32_t *data_u32 = (uint32_t *)serialized;

    if ( 2*sizeof(uint32_t) + 
         (data_u32[0]+data_u32[1])*sizeof(uint8_t) != 
         serialized_size)
        errx(1,
             "Malformed wah_non_leading serialized value. "
             "Incorrect serialized_size.");

    struct wah_bpt_non_leading_data *d = 
            (struct wah_bpt_non_leading_data *)
                    calloc(1, sizeof(struct wah_bpt_non_leading_data));
    d->SA = NULL;
    d->SE = NULL;

    uint8_t *data_8 = (uint8_t *)(data_u32 + 2);

    if (data_u32[0] > 0)
        d->SA = wah_copy(data_8);

    data_8 = data_8 + data_u32[0];

    if (data_u32[1] > 0)
        d->SE = wah_copy(data_8);


    *deserialized = d;

    return sizeof(struct wah_bpt_non_leading_data);
}
//}}}

//{{{void wah_non_leading_free(void **deserialized)
void wah_non_leading_free(void **deserialized)
{
    struct wah_bpt_non_leading_data **d = 
            (struct wah_bpt_non_leading_data **)deserialized;
    if ((*d)->SA != NULL)
        free((*d)->SA);
    if ((*d)->SE != NULL)
        free((*d)->SE);
    free(*d);
    *d = NULL;
}
//}}}
//}}}

//{{{ wah_leading_cache_handler 
struct cache_handler wah_leading_cache_handler = {
        wah_leading_serialize,
        wah_leading_deserialize,
        wah_leading_free
};


//{{{uint64_t wah_leading_serialize(void *deserialized,
uint64_t wah_leading_serialize(void *deserialized,
                                 void **serialized)
{
    if (deserialized == NULL) {
        *serialized = NULL;
        return 0;
    }

    struct wah_bpt_leading_data *d =  
            (struct wah_bpt_leading_data *)deserialized;

    uint32_t B_len = 0, serialized_len;

    if (d->B != NULL)
        B_len = sizeof(uint32_t) + 
                WAH_LEN(d->B)*(WAH_SIZE/BYTE)*sizeof(uint8_t);

    serialized_len = sizeof(uint32_t) + B_len;

    uint8_t *data = (uint8_t *)malloc(serialized_len);

    uint32_t *data_u = (uint32_t *)data;
    data_u[0] = B_len;

    uint32_t data_i = sizeof(uint32_t);

    if (d->B != NULL) 
        memcpy(data + data_i, d->B, B_len);

    data_i += B_len;

    if (data_i != serialized_len)
        errx(1,
             "Issue with wah_leading_serlize lengths. "
             "Expected:%u observed:%u.",
             serialized_len,
             data_i);

    *serialized = data;

    return serialized_len; 
}
//}}}

//{{{uint64_t wah_leading_deserialize(void *serialized,
uint64_t wah_leading_deserialize(void *serialized,
                                   uint64_t serialized_size,
                                   void **deserialized)
{
    if ((serialized_size == 0) || (serialized == NULL)) {
        *deserialized = NULL;
        return 0;
    }

    if (serialized_size < sizeof(uint32_t))
        errx(1,
             "Malformed wah_leading serialized value. "
             "Too short");

    uint32_t *data_u32 = (uint32_t *)serialized;

    if ( sizeof(uint32_t) + 
         (data_u32[0])*sizeof(uint8_t) != serialized_size)
        errx(1,
             "Malformed wah_leading serialized value. "
             "Incorrect serialized_size.");

    struct wah_bpt_leading_data *d = 
            (struct wah_bpt_leading_data *)
                    calloc(1, sizeof(struct wah_bpt_non_leading_data));
    d->B = NULL;

    uint8_t *data_8 = (uint8_t *)(data_u32 + 1);

    if (data_u32[0] > 0)
        d->B = wah_copy(data_8);

    *deserialized = d;

    return sizeof(struct wah_bpt_leading_data);
}
//}}}

//{{{void wah_leading_free(void **deserialized)
void wah_leading_free(void **deserialized)
{
    struct wah_bpt_leading_data **d = 
            (struct wah_bpt_leading_data **)deserialized;
    if ( (*d)->B != NULL )
        free ((*d)->B);
    free(*d);
    *d = NULL;
}
//}}}
//}}}

//{{{uint8_t *uints_to_wah(uint32_t *D, uint32_t D_num)
uint8_t *uints_to_wah(uint32_t *D, uint32_t D_num)
{
    uint32_t bits_per_word = WAH_SIZE - 1;
    uint32_t curr_word = 0, // num of words previously considered
             curr_val = 0, // value at the current index
             word_i = 0, // index into the array of words
             dist, fill_size, first_val, last = 0;
    uint32_t val, i;


    uint32_t w_num = D_num*2;
    uint8_t *w = (uint8_t *)malloc(sizeof(uint32_t) + 
                                   (w_num * (WAH_SIZE/BYTE) * sizeof(uint8_t)));

    // loop over the sorted input
    for (i = 0 ; i < D_num; ++i) {
        // get the distance from the current value and the first value in the
        // current word
        val = D[i] - (curr_word * bits_per_word);

        // will the val fit in the current word?
        if (val <= bits_per_word) {
            curr_val |= 1 << ( bits_per_word - val);
        } else  {
            if (curr_val > 0) {
                set_wah_i(w, &curr_val, WAH_SIZE, word_i);
                //fprintf(stderr,"curr_val:%u\tword_i:%u\n", curr_val, word_i);

                curr_word += 1; // move to the next word
                word_i += 1;
                curr_val = 0;

                if (word_i > w_num) {
                    w_num *= 2;
                    w = (uint8_t *)realloc(w,
                                           sizeof(uint32_t) + 
                                           (w_num * 
                                           (WAH_SIZE/BYTE) * 
                                           sizeof(uint8_t)));
                }

                val = D[i] - (curr_word * bits_per_word);
            }

            uint32_t saved_words;
            while (val > bits_per_word) {
                fill_size = ((val + bits_per_word - 1) / bits_per_word) - 1;
                saved_words = MIN(fill_size, WAH_MAX_FILL_WORDS);
                curr_val = (1 << (bits_per_word)) + (saved_words);

                set_wah_i(w, &curr_val, WAH_SIZE, word_i);
                //fprintf(stderr,"curr_val:%u\tword_i:%u\n", curr_val, word_i);

                curr_word += saved_words; // move to the next word
                word_i += 1;
                curr_val = 0;
                if (word_i > w_num) {
                    w_num *= 2;
                    w = (uint8_t *)realloc(w,
                                           sizeof(uint32_t) + 
                                           (w_num * 
                                           (WAH_SIZE/BYTE) * 
                                           sizeof(uint8_t)));
                }
                
                val -= saved_words * bits_per_word;
            }

            if (val > 0) {
                curr_val = 1 << ( bits_per_word - val);
            } else {
                curr_val = 0;
            }
        }
    }

    if (curr_val > 0)  {
        set_wah_i(w, &curr_val, WAH_SIZE, word_i);
        //fprintf(stderr,"curr_val:%u\tword_i:%u\n", curr_val, word_i);
    }

    w = (uint8_t *)realloc(w,
                           sizeof(uint32_t) + 
                           ((word_i + 1) * 
                           (WAH_SIZE/BYTE) * 
                           sizeof(uint8_t)));

    WAH_LEN(w) = word_i + 1;
    //fprintf(stderr, "WAH_LEN:%u\n", WAH_LEN(w));
    return w;
}
//}}}
