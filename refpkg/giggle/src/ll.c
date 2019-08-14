#define _GNU_SOURCE

#include <stdlib.h>
#include <stdint.h>
#include <err.h>
#include "util.h"
#include "ll.h"
#include "giggle_index.h"

//{{{void uint64_t_ll_map_intersection_to_offset_list(struct giggle_index *gi,
void uint64_t_ll_map_intersection_to_offset_list(struct giggle_index *gi,
                                                 struct giggle_query_result *gqr,
                                                 void *_R)
{
#ifdef DEBUG
    fprintf(stderr,
            "uint64_t_map_intersection_to_offset_list\n");
#endif

    struct uint64_t_ll *R  = (struct uint64_t_ll *)_R;

    if (R != NULL) {
        struct uint64_t_ll_node *curr = R->head;

#ifdef DEBUG
        fprintf(stderr,
                "giggle_query R->len:%u\n",
                R->len);
#endif

        while (curr != NULL) {
            struct file_id_offset_pair fid_off = 
                    offset_index_get(gi->offset_idx, curr->val);
                    //gi->offset_index->vals[curr->val];
            //long_ll_append(&(gqr->offsets[fid_off.file_id]),fid_off.offset);
            long_uint_ll_append(&(gqr->offsets[fid_off.file_id]),
                                fid_off.offset,
                                curr->val);
            curr = curr->next;
        }

        uint64_t_ll_free((void **)&R);
        R=NULL;
    } 
}
//}}}

//{{{ uint64_t_ll_non_leading_cache_handler
struct cache_handler uint64_t_ll_non_leading_cache_handler = {
        uint64_t_ll_non_leading_serialize,
        uint64_t_ll_non_leading_deserialize,
        uint64_t_ll_non_leading_free
};

//{{{uint64_t uint64_t_ll_non_leading_serialize(void *deserialized,
uint64_t uint64_t_ll_non_leading_serialize(void *deserialized,
                                           void **serialized)
{
    if (deserialized == NULL) {
        *serialized = NULL;
        return 0;
    }

    struct uint64_t_ll_bpt_non_leading_data *d =  
            (struct uint64_t_ll_bpt_non_leading_data *)deserialized;

    uint64_t SA_len, SE_len, serialized_len;

    if (d->SA != NULL)
        SA_len = d->SA->len;
    else
        SA_len = 0;

    if (d->SE != NULL)
        SE_len = d->SE->len;
    else
        SE_len = 0;


    serialized_len = (2 + SA_len + SE_len) * sizeof(uint64_t);

    uint64_t *data = (uint64_t *)calloc(2 + SA_len + SE_len, sizeof(uint64_t));
    if (data == NULL)
        err(1, "calloc error in uint64_t_ll_non_leading_serialize.\n");


    uint64_t data_i = 0;
    data[data_i++] = SA_len;
    data[data_i++] = SE_len;

    uint64_t i;
    struct uint64_t_ll_node *curr;

    if (d->SA != NULL) {
        curr = d->SA->head;
        i = 0;
        while (curr != NULL) {
            data[data_i++] = curr->val;
            i += 1;
            curr = curr->next;
        }

        if (i != d->SA->len)
            errx(1, "Incorrect length in SA:non_leading_data.");
    }

    if (d->SE != NULL) {
        curr = d->SE->head;
        i = 0;
        while (curr != NULL) {
            data[data_i++] = curr->val;
            i += 1;
            curr = curr->next;
        }

        if (i != d->SE->len)
            errx(1, "Incorrect length in SE:non_leading_data.");
    }

    *serialized = data;

    return serialized_len; 
}
//}}}

//{{{uint64_t uint64_t_ll_non_leading_deserialize(void *serialized,
uint64_t uint64_t_ll_non_leading_deserialize(void *serialized,
                                             uint64_t serialized_size,
                                             void **deserialized)
{
    if ((serialized_size == 0) || (serialized == NULL)) {
        *deserialized = NULL;
        return 0;
    }

    if (serialized_size < 2)
        errx(1,
             "Malformed uint64_t_ll_non_leading serialized value. "
             "Too short");

    uint64_t *data = (uint64_t *)serialized;

    if ((2 + data[0] + data[1]) * sizeof(uint64_t) != serialized_size)
        errx(1,
             "Malformed uint64_t_ll_non_leading serialized value. "
             "Incorrect serialized_size.");

    struct uint64_t_ll_bpt_non_leading_data *d = 
            (struct uint64_t_ll_bpt_non_leading_data *)
                    calloc(1,
                           sizeof(struct uint64_t_ll_bpt_non_leading_data));
    if (d == NULL)
        err(1, "calloc error in uint64_t_ll_non_leading_deserialize");

    d->SA = NULL;
    d->SE = NULL;

    uint64_t i;
    for (i = 0; i < data[0]; ++i)
        uint64_t_ll_append(&(d->SA), data[2 + i]);

    for (i = 0; i < data[1]; ++i)
        uint64_t_ll_append(&(d->SE), data[2 + data[0]+ i]);

    *deserialized = d;

    return sizeof(struct uint64_t_ll_bpt_non_leading_data);
}
//}}}

//{{{void uint64_t_ll_non_leading_free(void **deserialized)
void uint64_t_ll_non_leading_free(void **deserialized)
{
    struct uint64_t_ll_bpt_non_leading_data **d = 
            (struct uint64_t_ll_bpt_non_leading_data **)deserialized;
    uint64_t_ll_free((void **)&((*d)->SA));
    uint64_t_ll_free((void **)&((*d)->SE));
    free(*d);
    *d = NULL;
}
//}}}
//}}}

//{{{ uint64_t_ll_leading_cache_handler 
struct cache_handler uint64_t_ll_leading_cache_handler = {
        uint64_t_ll_leading_serialize,
        uint64_t_ll_leading_deserialize,
        uint64_t_ll_leading_free
};

//{{{uint64_t uint64_t_ll_leading_serialize(void *deserialized,
uint64_t uint64_t_ll_leading_serialize(void *deserialized,
                                       void **serialized)
{
    if (deserialized == NULL) {
        *serialized = NULL;
        return 0;
    }

    struct uint64_t_ll_bpt_leading_data *d =  
            (struct uint64_t_ll_bpt_leading_data *)deserialized;

    uint64_t B_len, serialized_len;

    if (d->B != NULL)
        B_len = d->B->len;
    else
        B_len = 0;

    serialized_len = (1 + B_len) * sizeof(uint64_t);

    uint64_t *data = (uint64_t *)calloc(1 + B_len , sizeof(uint64_t));
    if (data == NULL)
        err(1, "calloc error in uint64_t_ll_leading_serialize().\n");

    uint64_t data_i = 0;
    data[data_i++] = B_len;

    uint64_t i;
    struct uint64_t_ll_node *curr;

    if (d->B != NULL) {
        curr = d->B->head;
        i = 0;
        while (curr != NULL) {
            data[data_i++] = curr->val;
            i += 1;
            curr = curr->next;
        }

        if (i != d->B->len)
            errx(1, "Incorrect length in B:leading_data.");
    }

    *serialized = data;

    return serialized_len; 
}
//}}}

//{{{uint64_t uint64_t_ll_leading_deserialize(void *serialized,
uint64_t uint64_t_ll_leading_deserialize(void *serialized,
                                         uint64_t serialized_size,
                                         void **deserialized)
{
    if ((serialized_size == 0) || (serialized == NULL)) {
        *deserialized = NULL;
        return 0;
    }

    uint64_t *data = (uint64_t *)serialized;

    if ((1 + data[0]) * sizeof(uint64_t) != serialized_size)
        errx(1,
             "Malformed uint64_t_ll_leading serialized value. "
             "Incorrect serialized_size.");

    struct uint64_t_ll_bpt_leading_data *d = 
            (struct uint64_t_ll_bpt_leading_data *)
                    calloc(1,
                           sizeof(struct uint64_t_ll_bpt_leading_data));
    if (d == NULL)
        err(1, "calloc error in uint64_t_ll_leading_deserialize.\n");

    d->B = NULL;

    uint64_t i;
    for (i = 0; i < data[0]; ++i)
        uint64_t_ll_append(&(d->B), data[1 + i]);

    *deserialized = d;
    return sizeof(struct uint64_t_ll_bpt_leading_data);
}
//}}}

//{{{void uint64_t_ll_leading_free(void **deserialized)
void uint64_t_ll_leading_free(void **deserialized)
{
    struct uint64_t_ll_bpt_leading_data **d = 
            (struct uint64_t_ll_bpt_leading_data **)deserialized;
    uint64_t_ll_free((void **)&((*d)->B));
    free(*d);
    *d = NULL;
}
//}}}
//}}}

//{{{ giggle_data_handler :: int64_t_ll_giggle_data_handler

void uint64_t_ll_giggle_set_data_handler()
{
    bpt_node_repair = uint64_t_ll_leading_repair;

    uint64_t_ll_giggle_data_handler.non_leading_cache_handler =
        uint64_t_ll_non_leading_cache_handler;
    uint64_t_ll_giggle_data_handler.leading_cache_handler = 
        uint64_t_ll_leading_cache_handler;
    uint64_t_ll_giggle_data_handler.new_non_leading = 
        uint64_t_ll_new_non_leading;
    uint64_t_ll_giggle_data_handler.new_leading = 
        uint64_t_ll_new_leading;
    uint64_t_ll_giggle_data_handler.non_leading_SA_add_scalar = 
        uint64_t_ll_non_leading_SA_add_scalar;
    uint64_t_ll_giggle_data_handler.non_leading_SE_add_scalar = 
        uint64_t_ll_non_leading_SE_add_scalar;
    uint64_t_ll_giggle_data_handler.leading_B_add_scalar = 
        uint64_t_ll_leading_B_add_scalar;
    uint64_t_ll_giggle_data_handler.leading_union_with_B = 
        uint64_t_ll_leading_union_with_B;
    uint64_t_ll_giggle_data_handler.non_leading_union_with_SA = 
        uint64_t_ll_non_leading_union_with_SA;
    uint64_t_ll_giggle_data_handler.non_leading_union_with_SA_subtract_SE = 
        uint64_t_ll_non_leading_union_with_SA_subtract_SE;

    uint64_t_ll_giggle_data_handler.write_tree = 
        giggle_write_tree_cache_dump;

    uint64_t_ll_giggle_data_handler.giggle_collect_intersection =
        giggle_collect_intersection_data_in_pointers;

    uint64_t_ll_giggle_data_handler.map_intersection_to_offset_list =
        uint64_t_ll_map_intersection_to_offset_list;

    giggle_data_handler = uint64_t_ll_giggle_data_handler;
}


//{{{void *uint64_t_ll_new_non_leading()
void *uint64_t_ll_new_non_leading(uint32_t domain)
{
    struct uint64_t_ll_bpt_non_leading_data *d = 
            (struct uint64_t_ll_bpt_non_leading_data *)
            malloc(sizeof( struct uint64_t_ll_bpt_non_leading_data));
    if (d == NULL)
        err(1, "malloc error in uint64_t_ll_bpt_non_leading_data.\n");

    d->SA = NULL;
    d->SE = NULL;

    return (void *)d;
}
//}}}

//{{{void *uint64_t_ll_new_leading()
void *uint64_t_ll_new_leading(uint32_t domain)
{
    struct uint64_t_ll_bpt_leading_data *d = 
            (struct uint64_t_ll_bpt_leading_data *)
            malloc(sizeof( struct uint64_t_ll_bpt_leading_data));
    if (d == NULL)
        err(1, "malloc error in uint64_t_ll_new_leading.\n");

    d->B = NULL;

    return (void *)d;
}
//}}}

//{{{void uint64_t_ll_non_leading_SA_add_scalar(void *d, void *v)
void uint64_t_ll_non_leading_SA_add_scalar(uint32_t domain,
                                           void *_nld,
                                           void *_id)
{
#if DEBUG
    fprintf(stderr, "uint64_t_ll_non_leading_SA_add_scalar\n");
#endif
    struct uint64_t_ll_bpt_non_leading_data *nld =
            (struct uint64_t_ll_bpt_non_leading_data *)_nld;
    uint64_t *id = (uint64_t *)_id;

#if DEBUG
    fprintf(stderr, "id:%u\n", *id);
#endif

    uint64_t_ll_uniq_append(&(nld->SA), *id);
}
//}}}

//{{{void uint64_t_ll_non_leading_SE_add_scalar(void *d, void *v)
void uint64_t_ll_non_leading_SE_add_scalar(uint32_t domain,
                                           void *_nld,
                                           void *_id)
{
#if DEBUG
    fprintf(stderr, "uint64_t_ll_non_leading_SE_add_scalar\n");
#endif
    struct uint64_t_ll_bpt_non_leading_data *nld =
            (struct uint64_t_ll_bpt_non_leading_data *)_nld;
    uint64_t *id = (uint64_t *)_id;

#if DEBUG
    fprintf(stderr, "id:%u\n", *id);
#endif

    uint64_t_ll_uniq_append(&(nld->SE), *id);
}
//}}}

//{{{void uint64_t_ll_leading_B_add_scalar(void *d, void *v)
void uint64_t_ll_leading_B_add_scalar(uint32_t domain,
                                      void *_ld,
                                      void *_id)
{
#if DEBUG
    fprintf(stderr, "uint64_t_ll_leading_B_add_scalar\n");
#endif
    struct uint64_t_ll_bpt_leading_data *ld =
            (struct uint64_t_ll_bpt_leading_data *)_ld;
    uint64_t *id = (uint64_t *)_id;

#if DEBUG
    fprintf(stderr, "id:%u\n", *id);
#endif

    uint64_t_ll_uniq_append(&(ld->B), *id);
}
//}}}

//{{{void uint64_t_ll_leading_union_with_B(void **R, void *leading)
void uint64_t_ll_leading_union_with_B(uint32_t domain,
                                      void **R,
                                      void *leading)
{
#if DEBUG
    fprintf(stderr, "uint64_t_ll_leading_union_with_B\n");
#endif

    struct uint64_t_ll_node *curr = 
        ((struct uint64_t_ll_bpt_leading_data *)(leading))->B->head;
    while (curr != NULL) {
#if DEBUG
        fprintf(stderr, "+%u\n", curr->val);
#endif
        uint64_t_ll_uniq_append((struct uint64_t_ll **)R, curr->val);
        curr = curr->next;
    }
}
//}}}

//{{{void uint64_t_ll_non_leading_union_with_SA(void **R, void *d)
void uint64_t_ll_non_leading_union_with_SA(uint32_t domain, void **R, void *d)
{
#if DEBUG
    fprintf(stderr, "uint64_t_ll_non_leading_union_with_SA\n");
#endif
    struct uint64_t_ll_bpt_non_leading_data *nld = 
            (struct uint64_t_ll_bpt_non_leading_data *) d;
    if (nld != NULL) {
        if ((nld->SA != NULL)) {
            struct uint64_t_ll_node *curr = nld->SA->head;
            while (curr != NULL) {
#if DEBUG
                fprintf(stderr, "+%u\n", curr->val);
#endif
                uint64_t_ll_uniq_append((struct uint64_t_ll **)R, curr->val);
                curr = curr->next;
            }
        }
    }
}
//}}}

//{{{void uint64_t_ll_non_leading_union_with_SA_subtract_SE(void **R, void *d)
void uint64_t_ll_non_leading_union_with_SA_subtract_SE(uint32_t domain,
                                                       void **R,
                                                       void *d)
{
#if DEBUG
    fprintf(stderr, "uint64_t_ll_non_leading_union_with_SA_subtract_SE\n");
#endif

    struct uint64_t_ll_bpt_non_leading_data *nld = 
            (struct uint64_t_ll_bpt_non_leading_data *) d;
    if (nld != NULL) {
        if ((nld->SA != NULL)) {
            struct uint64_t_ll_node *curr = nld->SA->head;
            while (curr != NULL) {
#if DEBUG
                fprintf(stderr, "+%u\n", curr->val);
#endif
                uint64_t_ll_uniq_append((struct uint64_t_ll **)R, curr->val);
                curr = curr->next;
            }
        }
        if ((nld->SE != NULL)) {
            struct uint64_t_ll_node *curr = nld->SE->head;
            while (curr != NULL) {
#if DEBUG
                fprintf(stderr, "-%u\n", curr->val);
#endif
                uint64_t_ll_remove((struct uint64_t_ll **)R, curr->val);
                curr = curr->next;
            }
        }
    }

}
//}}}
//}}}

//{{{void uint64_t_ll_append(struct uint64_t_ll **ll, uint64_t val)
void uint64_t_ll_append(struct uint64_t_ll **ll, uint64_t val)
{
    struct uint64_t_ll_node *n = (struct uint64_t_ll_node *)
        malloc(sizeof(struct uint64_t_ll_node));
    if (n == NULL)
        err(1, "malloc error in uint64_t_ll_append.\n");
    n->val = val;
    n->next = NULL;

    if (*ll == NULL) {
        *ll = (struct uint64_t_ll *)malloc(sizeof(struct uint64_t_ll));
        if (*ll == NULL)
            err(1, "malloc error in uint64_t_ll_append.\n");

        (*ll)->head = n;
        (*ll)->len = 1;
    } else {
        (*ll)->tail->next = n;
        (*ll)->len = (*ll)->len + 1;
    }

    (*ll)->tail = n;
}
//}}}

//{{{void uint64_t_ll_uniq_append(struct uint64_t_ll **ll, uint64_t val)
void uint64_t_ll_uniq_append(struct uint64_t_ll **ll, uint64_t val)
{
#if DEBUG
    fprintf(stderr, "uint64_t_ll_uniq_append\n");
    fprintf(stderr, "val:%u\n", val);
#endif
    struct uint64_t_ll_node *n = (struct uint64_t_ll_node *)
        malloc(sizeof(struct uint64_t_ll_node));
    if (n == NULL)
        err(1, "malloc error in uint64_t_ll_uniq_append.\n");

    n->val = val;
    n->next = NULL;

    if (*ll == NULL) {
        *ll = (struct uint64_t_ll *)malloc(sizeof(struct uint64_t_ll));
        if (*ll == NULL)
            err(1, "malloc error in uint64_t_ll_uniq_append.\n");

        (*ll)->head = n;
        (*ll)->len = 1;
    } else {

        struct uint64_t_ll_node *curr = (*ll)->head;
        while (curr != NULL) {
            if (curr->val == val) {
                free(n);
                return;
            }
            curr = curr->next;
        }

        (*ll)->tail->next = n;
        (*ll)->len = (*ll)->len + 1;
    }

    (*ll)->tail = n;
}
//}}}

//{{{void uint64_t_ll_remove(struct uint64_t_ll **ll, uint64_t val)
void uint64_t_ll_remove(struct uint64_t_ll **ll, uint64_t val)
{
    if (*ll != NULL) {
        struct uint64_t_ll_node *tmp, *last = NULL, *curr = (*ll)->head;
        while (curr != NULL) {
            if (curr->val == val) {
                if ((curr == (*ll)->head) && (curr == (*ll)->tail)) {
                    free(curr);
                    free(*ll);
                    *ll = NULL;
                    return;
                } else if (curr == (*ll)->head) {
                    tmp = curr->next;
                    free(curr);
                    (*ll)->head = tmp;
                    curr = tmp;
                    (*ll)->len = (*ll)->len - 1;
                } else if (curr == (*ll)->tail) {
                    free(curr);
                    curr = NULL;
                    (*ll)->tail = last;
                    last->next = NULL;
                    (*ll)->len = (*ll)->len - 1;
                } else {
                    tmp = curr;
                    last->next = tmp->next;
                    curr = tmp->next;
                    free(tmp);
                    (*ll)->len = (*ll)->len - 1;
                }
            } else {
                last = curr;
                curr = curr->next;
            }
        }
    }
}
//}}}

//{{{int uint64_t_ll_contains(struct uint64_t_ll *ll, uint64_t val)
uint64_t uint64_t_ll_contains(struct uint64_t_ll *ll, uint64_t val)
{
    if (ll == NULL)
       return 0; 

    uint64_t r = 0;
    
    struct uint64_t_ll_node *curr = ll->head;
    while (curr != NULL) {
        if (curr->val == val)
            r += 1;
        curr = curr->next;
    }

    return r;
}
//}}}

//{{{void uint64_t_ll_free(struct uint64_t_ll **ll)
//void uint64_t_ll_free(struct uint64_t_ll **ll)
void uint64_t_ll_free(void **_ll)
{
    struct uint64_t_ll **ll = (struct uint64_t_ll **)_ll;
    struct uint64_t_ll_node *curr, *tmp;

    if (*ll != NULL) {
        curr = (*ll)->head;
        while (curr != NULL) {
            tmp = curr->next;
            free(curr);
            curr = tmp;
        }

        free(*ll);
        *ll = NULL;
    }
}
//}}}

//{{{void leading_repair(struct bpt_node *a, struct bpt_node *b)
void uint64_t_ll_leading_repair(uint32_t domain, 
                                struct bpt_node *a,
                                struct bpt_node *b)
{
#if DEBUG
    fprintf(stderr, "leading_repair\n");
#endif

    if ( (BPT_IS_LEAF(a) == 1) && (BPT_IS_LEAF(b) == 1) ) {
        struct uint64_t_ll_bpt_leading_data *d = 
            (struct uint64_t_ll_bpt_leading_data *)
            malloc(sizeof(struct uint64_t_ll_bpt_leading_data));
        if (d == NULL)
            err(1, "malloc error in uint64_t_ll_leading_repair.\n");

        d->B = NULL;

        if (BPT_LEADING(a) != 0) {
            struct uint64_t_ll_bpt_leading_data *l =  
                    cache.get(domain,
                              BPT_LEADING(a) - 1,
                              &uint64_t_ll_leading_cache_handler);

            struct uint64_t_ll_node *curr = l->B->head;

            while (curr != NULL) {
#if DEBUG
                fprintf(stderr, "+ %u\n", curr->val);
#endif
                uint64_t_ll_uniq_append(&(d->B), curr->val);
                curr = curr->next;
            }
        }

        uint64_t i;
        for (i = 0 ; i < BPT_NUM_KEYS(a); ++i) {
            struct uint64_t_ll_bpt_non_leading_data *nl = 
                    cache.get(domain,
                              BPT_POINTERS(a)[i] - 1,
                              &uint64_t_ll_non_leading_cache_handler);

            if (nl->SA != NULL) {
                struct uint64_t_ll_node *curr = nl->SA->head;
                while (curr != NULL) {
                    uint64_t_ll_uniq_append(&(d->B), curr->val);
#if DEBUG
                    fprintf(stderr, "+ %u\n", curr->val);
#endif
                    curr = curr->next;
                }
            }

            if (nl->SE != NULL) {
                struct uint64_t_ll_node *curr = nl->SE->head;
                while (curr != NULL) {
                    uint64_t_ll_remove(&(d->B), curr->val);
#if DEBUG
                    fprintf(stderr, "- %u\n", curr->val);
#endif
                    curr = curr->next;
                }
            }
        }

        if (d->B != NULL) {
            uint64_t v_id = cache.seen(domain) + 1;
            cache.add(domain,
                      v_id - 1,
                      d,
                      sizeof(struct uint64_t_ll_bpt_leading_data),
                      &uint64_t_ll_leading_cache_handler);

            BPT_LEADING(b) = v_id;
        } else {
            free(d);
        }
    }
}
//}}}

//{{{void long_ll_append(struct long_ll **ll, long val)
void long_ll_append(struct long_ll **ll, long val)
{
    struct long_ll_node *n = (struct long_ll_node *)
        malloc(sizeof(struct long_ll_node));
    if (n == NULL)
        err(1, "malloc error in long_ll_append().\n");

    n->val = val;
    n->next = NULL;

    if (*ll == NULL) {
        *ll = (struct long_ll *)malloc(sizeof(struct long_ll));
        if (*ll == NULL)
            err(1, "malloc error in long_ll_append().\n");
        (*ll)->head = n;
        (*ll)->len = 1;
    } else {
        (*ll)->tail->next = n;
        (*ll)->len = (*ll)->len + 1;
    }

    (*ll)->tail = n;
}
//}}}

//{{{void long_ll_uniq_append(struct long_ll **ll, long val)
void long_ll_uniq_append(struct long_ll **ll, long val)
{
#if DEBUG
    fprintf(stderr, "long_ll_uniq_append\n");
    fprintf(stderr, "val:%lu\n", val);
#endif
    struct long_ll_node *n = (struct long_ll_node *)
        malloc(sizeof(struct long_ll_node));
    if (n == NULL)
        err(1, "malloc error in long_ll_uniq_append().\n");
    n->val = val;
    n->next = NULL;

    if (*ll == NULL) {
        *ll = (struct long_ll *)malloc(sizeof(struct long_ll));
        if (*ll == NULL)
            err(1, "malloc error in long_ll_uniq_append().\n");
        (*ll)->head = n;
        (*ll)->len = 1;
    } else {

        struct long_ll_node *curr = (*ll)->head;
        while (curr != NULL) {
            if (curr->val == val) {
                free(n);
                return;
            }
            curr = curr->next;
        }

        (*ll)->tail->next = n;
        (*ll)->len = (*ll)->len + 1;
    }

    (*ll)->tail = n;
}
//}}}

//{{{void long_ll_remove(struct long_ll **ll, long val)
void long_ll_remove(struct long_ll **ll, long val)
{
    if (*ll != NULL) {
        struct long_ll_node *tmp, *last = NULL, *curr = (*ll)->head;
        while (curr != NULL) {
            if (curr->val == val) {
                if ((curr == (*ll)->head) && (curr == (*ll)->tail)) {
                    free(curr);
                    free(*ll);
                    *ll = NULL;
                    return;
                } else if (curr == (*ll)->head) {
                    tmp = curr->next;
                    free(curr);
                    (*ll)->head = tmp;
                    curr = tmp;
                    (*ll)->len = (*ll)->len - 1;
                } else if (curr == (*ll)->tail) {
                    free(curr);
                    curr = NULL;
                    (*ll)->tail = last;
                    last->next = NULL;
                    (*ll)->len = (*ll)->len - 1;
                } else {
                    tmp = curr;
                    last->next = tmp->next;
                    curr = tmp->next;
                    free(tmp);
                    (*ll)->len = (*ll)->len - 1;
                }
            } else {
                last = curr;
                curr = curr->next;
            }
        }
    }
}
//}}}

//{{{int long_ll_contains(struct long_ll *ll, long val)
uint64_t long_ll_contains(struct long_ll *ll, long val)
{
    if (ll == NULL)
       return 0; 

    uint64_t r = 0;
    
    struct long_ll_node *curr = ll->head;
    while (curr != NULL) {
        if (curr->val == val)
            r += 1;
        curr = curr->next;
    }

    return r;
}
//}}}

//{{{void long_ll_free(struct long_ll **ll)
//void long_ll_free(struct long_ll **ll)
void long_ll_free(void **_ll)
{
    struct long_ll **ll = (struct long_ll **)_ll;
    struct long_ll_node *curr, *tmp;

    if (*ll != NULL) {
        curr = (*ll)->head;
        while (curr != NULL) {
            tmp = curr->next;
            free(curr);
            curr = tmp;
        }

        free(*ll);
        *ll = NULL;
    }
}
//}}}

//{{{void long_uint_ll_append(struct long_uint_ll **ll, long long_val, uint64_t
void long_uint_ll_append(struct long_uint_ll **ll, long long_val, uint64_t uint_val)
{
#if GIGGLE_QUERY_TRACE
    fprintf(stderr,
            "long_uint_ll_append\tlong_val:%lu\tuint_val:%llu\n",
            long_val,
            uint_val);
#endif
    struct long_uint_ll_node *n = (struct long_uint_ll_node *)
        malloc(sizeof(struct long_uint_ll_node));
#if GIGGLE_QUERY_TRACE
    fprintf(stderr,
            "long_uint_ll_append\tn:%p\n",
            n);
#endif
    if (n == NULL)
        err(1, "malloc error in long_uint_ll_append().\n");

    n->long_val = long_val;
    n->uint_val = uint_val;
    n->next = NULL;

    if (*ll == NULL) {
        *ll = (struct long_uint_ll *)malloc(sizeof(struct long_uint_ll));
        if (*ll == NULL)
            err(1, "malloc error in long_uint_ll_append().\n");
        (*ll)->head = n;
        (*ll)->len = 1;
    } else {
        (*ll)->tail->next = n;
        (*ll)->len = (*ll)->len + 1;
    }

    (*ll)->tail = n;
}
//}}}

//{{{void long_uint_ll_free(struct long_ll **ll)
void long_uint_ll_free(void **_ll)
{
    struct long_uint_ll **ll = (struct long_uint_ll **)_ll;
    struct long_uint_ll_node *curr, *tmp;

    if (*ll != NULL) {
        curr = (*ll)->head;
        while (curr != NULL) {
            tmp = curr->next;
            free(curr);
            curr = tmp;
        }

        free(*ll);
        *ll = NULL;
    }
}
//}}}
