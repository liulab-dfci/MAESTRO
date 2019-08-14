#define _GNU_SOURCE

#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <err.h>
#include <sysexits.h>
#include <stdbool.h>
#include "bpt.h"
#include "lists.h"
//#include "ll.h"
#include "util.h"
#include "disk_store.h"
#include "cache.h"

/*
 *  All bpt_nodes have an array of keys, and an array of values.
 *  In non-leaf bpt_nodes, values[i] points to the bpt_node that is either
 *  greater than or equal to key[i-1] or less than key[i].
 */

//uint32_t ORDER = 4;
uint32_t ORDER = 100;

void (*bpt_node_repair)(uint32_t domain,
                        struct bpt_node *,
                        struct bpt_node *) = NULL;

struct cache_handler bpt_node_cache_handler = {bpt_node_serialize, 
                                               bpt_node_deserialize,
                                               bpt_node_free_mem};

//{{{int b_search(int key, int *D, int D_size)
int b_search(uint32_t key, const uint32_t *D, uint32_t D_size)
{
    //fprintf(stderr, "D_size:%u\n", D_size);
    // This is a common case when incoming data is sorted.
    if(key > D[D_size-1]) {
        return D_size;
    }
    int lo = -1, hi = D_size, mid;
    while ( hi - lo > 1) {
        mid = (hi + lo) / 2;
        if ( D[mid] < key )
            lo = mid;
        else
            hi = mid;
    }

    return hi;
}
//}}}

//{{{ bpt_node_cache_handler
//{{{uint64_t bpt_node_serialize(void *deserialized, void **serialized)
uint64_t bpt_node_serialize(void *deserialized, void **serialized)
{
    struct bpt_node *de = (struct bpt_node *)deserialized;

    //uint8_t *data = (uint8_t *)calloc(2*ORDER+9, sizeof(uint32_t));
    uint8_t *data = (uint8_t *)calloc(BPT_NODE_NUM_ELEMENTS*sizeof(uint32_t),
                                      sizeof(uint8_t));
    if (data == NULL)
        err(1, "calloc error in bpt_node_serialize()");

    memcpy(data, de->data, BPT_NODE_NUM_ELEMENTS*sizeof(uint32_t));

    *serialized = (void *)data;
    return BPT_NODE_NUM_ELEMENTS*sizeof(uint32_t);
}
//}}}

//{{{ uint64_t bpt_node_deserialize(void *serialized,
uint64_t bpt_node_deserialize(void *serialized,
                             uint64_t serialized_size,
                             void **deserialized)
{
    uint8_t *data = (uint8_t *)serialized;

    struct bpt_node *n = (struct bpt_node *)malloc(sizeof(struct bpt_node));

    if (n == NULL)
        err(1, "malloc error in bpt_node_deserialize().");

    n->data = (uint32_t *)calloc(BPT_NODE_NUM_ELEMENTS,sizeof(uint32_t));
    if (n->data == NULL)
        err(1, "calloc error in bpt_node_deserialize().");

    memcpy(n->data, data, BPT_NODE_NUM_ELEMENTS*sizeof(uint32_t));

    *deserialized = (void *)n;

    return sizeof(struct bpt_node);
}
//}}}

//{{{void bpt_node_free_mem(void **deserialized)
void bpt_node_free_mem(void **deserialized)
{
    struct bpt_node **de = (struct bpt_node **)deserialized;
    free((*de)->data);
    free(*de);
    *de = NULL;
}
//}}}
//}}}

//{{{ struct bpt_node *bpt_new_node()
struct bpt_node *bpt_new_node(uint32_t domain)
{
    struct bpt_node *n = (struct bpt_node *)malloc(sizeof(struct bpt_node));
    if (n == NULL)
        err(1, "malloc error in bpt_new_node().");

    n->data = (uint32_t *)calloc(BPT_NODE_NUM_ELEMENTS, sizeof(uint32_t));
    if (n->data == NULL)
        err(1, "calloc error in bpt_new_node().");

    // zero is null for a bpt, so the ids start at one the cache number starts
    // at zero so add or subtract one to get the one/zero-based numbering 
    BPT_ID(n) = cache.seen(domain) + 1;
    cache.add(domain,
              BPT_ID(n) - 1,
              n,
              sizeof(struct bpt_node),
              &bpt_node_cache_handler);

    return n;
}
//}}}

//{{{ struct bpt_node *bpt_find_leaf(struct bpt_node *curr, int key)
uint32_t bpt_find_leaf(uint32_t domain, uint32_t curr_id, uint32_t key)
{
    //fprintf(stderr, "domain:%u\tcurr_id:%u\n", domain, curr_id);
    // cache is zero-based, while bpt is one-based
    struct bpt_node *curr = cache.get(domain,
                                      curr_id - 1,
                                      &bpt_node_cache_handler);
    if (curr == NULL)
        return 0;

    while (BPT_IS_LEAF(curr) != 1) {
        int i = bpt_find_insert_pos(curr, key);

        if ((i < BPT_NUM_KEYS(curr)) && (BPT_KEYS(curr)[i] == key))
            i+=1;

        // cache is zero-based, while bpt is one-based
        uint32_t next = BPT_POINTERS(curr)[i] - 1;
        curr = cache.get(domain,
                         next,
                         &bpt_node_cache_handler);
    }

    uint32_t id = BPT_ID(curr);
    return id;
}
//}}}

//{{{int bpt_find_insert_pos(struct bpt_node *leaf, int key)
int bpt_find_insert_pos(struct bpt_node *leaf, uint32_t key)
{
    return b_search(key, BPT_KEYS(leaf), BPT_NUM_KEYS(leaf));
}
//}}}

//{{{struct bpt_node *bpt_split_node(struct bpt_node *root, struct bpt_node
uint32_t bpt_split_node(uint32_t domain,
                        uint32_t root_id,
                        uint32_t bpt_node_id,
                        uint32_t *lo_result_id,
                        uint32_t *hi_result_id,
                        int *lo_hi_split_point,
                        void (*repair)(uint32_t domain,
                                       struct bpt_node *,
                                       struct bpt_node *))
{
#if DEBUG
        fprintf(stderr, "bpt_split_node\n");
#endif

#if DEBUG
    {
        fprintf(stderr, "bpt_split_node\n");
        fprintf(stderr, "root_id:%u\tbpt_node_id:%u\n", root_id, bpt_node_id);
    }
#endif

    // cache is zero-based, while bpt is one-based
    struct bpt_node *bpt_node = cache.get(domain,
                                          bpt_node_id - 1,
                                          &bpt_node_cache_handler);
    struct bpt_node *n = bpt_new_node(domain);

#if DEBUG
    {
        int i;
        fprintf(stderr, "keys\t");
        for (i = 0; i < BPT_NUM_KEYS(bpt_node); ++i)
            fprintf(stderr, "%u\t", BPT_KEYS(bpt_node)[i]);
        fprintf(stderr, "\n");
    }
#endif

    // set the split location
    int split_point = (ORDER + 1)/2;

    *lo_result_id = BPT_ID(bpt_node);
    *hi_result_id = BPT_ID(n);
    *lo_hi_split_point = split_point;

    // copy the 2nd 1/2 of the values over to the new bpt_node
    int bpt_node_i, new_bpt_node_i = 0;
    for (bpt_node_i = split_point;
         bpt_node_i < BPT_NUM_KEYS(bpt_node);
         ++bpt_node_i) {

        BPT_KEYS(n)[new_bpt_node_i] = BPT_KEYS(bpt_node)[bpt_node_i];
        BPT_POINTERS(n)[new_bpt_node_i] = BPT_POINTERS(bpt_node)[bpt_node_i];
        BPT_NUM_KEYS(n) += 1;
        new_bpt_node_i += 1;
    }

    // if the bpt_node is not a leaf, the far right pointer must be coppied too
    if (BPT_IS_LEAF(bpt_node) == 0) {
        BPT_POINTERS(n)[new_bpt_node_i] = BPT_POINTERS(bpt_node)[bpt_node_i];
        BPT_POINTERS(n)[0] = 0;
    }

    // set status of new bpt_node
    BPT_IS_LEAF(n) = BPT_IS_LEAF(bpt_node);
    BPT_PARENT(n) = BPT_PARENT(bpt_node);

    BPT_NUM_KEYS(bpt_node) = split_point;

    if (BPT_IS_LEAF(bpt_node) == 0) {
        // if the bpt_node is not a leaf, then update the parent pointer of the 
        // children
        for (bpt_node_i = 1; bpt_node_i <= BPT_NUM_KEYS(n); ++bpt_node_i) {
            // cache is zero-based, while bpt is one-based
            struct bpt_node *child = cache.get(domain,
                                               BPT_POINTERS(n)[bpt_node_i] - 1,
                                               &bpt_node_cache_handler);
            BPT_PARENT(child) = BPT_ID(n);
        }
    } else {
        // if the bpt_node is a leaf, then connect the two
        BPT_NEXT(n)= BPT_NEXT(bpt_node);
        BPT_NEXT(bpt_node) = BPT_ID(n);
    }

    if (bpt_node_repair != NULL) {
        bpt_node_repair(domain, bpt_node, n);
    }

    if (BPT_ID(bpt_node) == root_id) {
        // if we just split the root, create a new root witha single value
        struct bpt_node *new_root = bpt_new_node(domain);
        BPT_IS_LEAF(new_root) = 0;
        BPT_NUM_KEYS(new_root) += 1;
        BPT_KEYS(new_root)[0] = BPT_KEYS(n)[0];
        BPT_POINTERS(new_root)[0] = BPT_ID(bpt_node); 
        BPT_POINTERS(new_root)[1] = BPT_ID(n); 
        BPT_PARENT(bpt_node) = BPT_ID(new_root);
        BPT_PARENT(n) = BPT_ID(new_root);
        return BPT_ID(new_root);
    } else {
        // if we didnt split the root, place the new value in the parent
        // bpt_node
        int trash_pos;
        uint32_t parent_id =  BPT_PARENT(bpt_node);
        return bpt_place_new_key_value(domain,
                                       root_id,
                                       &parent_id,
                                       &trash_pos,
                                       BPT_KEYS(n)[0],
                                       BPT_ID(n),
                                       NULL);
    }
}
//}}}

//{{{ uint32_t bpt_place_new_key_value(uint32_t domain,
uint32_t bpt_place_new_key_value(uint32_t domain,
                                 uint32_t root_id,
                                 uint32_t *target_id,
                                 int *target_key_pos,
                                 uint32_t key,
                                 uint32_t value_id,
                                 struct cache_handler *handler)
{
#if DEBUG
    {
        fprintf(stderr, "bpt_place_new_key_value\n");
        fprintf(stderr,
                "root_id:%u\ttarget_id:%u\tkey:%u\n",
                root_id,
                *target_id,
                key);
    }
#endif

    // cache is zero-based, while bpt is one-based
    struct bpt_node *target_bpt_node = cache.get(domain,
                                                 *target_id - 1,
                                                 &bpt_node_cache_handler);

#if DEBUG
    {
        int i;
        fprintf(stderr, "keys\t");
        for (i = 0; i < BPT_NUM_KEYS(target_bpt_node); ++i)
            fprintf(stderr, "%u\t", BPT_KEYS(target_bpt_node)[i]);
        fprintf(stderr, "\n");
    }
#endif

    int bpt_insert_key_pos = bpt_find_insert_pos(target_bpt_node, key);

    int bpt_insert_value_pos = bpt_insert_key_pos;

    if (BPT_IS_LEAF(target_bpt_node) == 0)
        bpt_insert_value_pos += 1;

    if (BPT_IS_LEAF(target_bpt_node) == 1)
        *target_key_pos = bpt_insert_key_pos;


    if ((BPT_IS_LEAF(target_bpt_node) == 1) &&
         (*target_key_pos < (BPT_NUM_KEYS(target_bpt_node))) &&
        (BPT_KEYS(target_bpt_node)[*target_key_pos] == key )) {

        // If the append function is NULL assume overwrite
        if (append != NULL)
            append(domain,
                   value_id,
                   BPT_POINTERS(target_bpt_node)[*target_key_pos],
                   handler);
        else
            BPT_POINTERS(target_bpt_node)[*target_key_pos] = value_id;

        return root_id;
    }

    // move everything over
    int i;
    for (i = BPT_NUM_KEYS(target_bpt_node); i > bpt_insert_key_pos; --i) {
        BPT_KEYS(target_bpt_node)[i] = BPT_KEYS(target_bpt_node)[i-1];
    }

    if (BPT_IS_LEAF(target_bpt_node) == 1) {
        for (i = BPT_NUM_KEYS(target_bpt_node); i > bpt_insert_value_pos; --i) 
            BPT_POINTERS(target_bpt_node)[i] = 
                    BPT_POINTERS(target_bpt_node)[i-1];
    } else {
        for (i = BPT_NUM_KEYS(target_bpt_node)+1;
             i > bpt_insert_value_pos;
             --i) {
            BPT_POINTERS(target_bpt_node)[i] = 
                    BPT_POINTERS(target_bpt_node)[i-1];
        }
    }

#if DEBUG
    {
        fprintf(stderr,
                "bpt_insert_key_pos:%u\tbpt_insert_value_pos:%u\n",
                bpt_insert_key_pos,
                bpt_insert_value_pos);
    }
#endif

    BPT_KEYS(target_bpt_node)[bpt_insert_key_pos] = key;
    BPT_POINTERS(target_bpt_node)[bpt_insert_value_pos] = value_id;

    BPT_NUM_KEYS(target_bpt_node) += 1;

    // If there are now too many values in the bpt_node, split it
    if (BPT_NUM_KEYS(target_bpt_node) > ORDER) {
        uint32_t lo_result_id, hi_result_id;
        int lo_hi_split_point = 0;
        uint32_t new_root_id = bpt_split_node(domain,
                                              root_id,
                                              BPT_ID(target_bpt_node),
                                              &lo_result_id,
                                              &hi_result_id,
                                              &lo_hi_split_point,
                                              bpt_node_repair);

        // cache is zero-based, while bpt is one-based
        target_bpt_node = cache.get(domain,
                                    *target_id - 1,
                                    &bpt_node_cache_handler);

        if (BPT_IS_LEAF(target_bpt_node)) {
            if (bpt_insert_key_pos < lo_hi_split_point)
                *target_id = lo_result_id;
            else {
                *target_id = hi_result_id;
                *target_key_pos = bpt_insert_key_pos - lo_hi_split_point;
            }
        }

        return new_root_id;
    } else {
        return root_id;
    }
}
//}}}

//{{{ uint32_t bpt_insert(uint32_t domain,
uint32_t bpt_insert(uint32_t domain,
                    uint32_t root_id,
                    uint32_t key,
                    uint32_t value_id,
                    struct cache_handler *handler,
                    uint32_t *leaf_id,
                    int *pos)
{
#if DEBUG
    fprintf(stderr, "bpt_insert\n");
#endif

    if (root_id == 0) {
        struct bpt_node *root = bpt_new_node(domain);
        BPT_IS_LEAF(root) = 1;
        BPT_KEYS(root)[0] = key;
        BPT_POINTERS(root)[0] = value_id;
        BPT_NUM_KEYS(root) += 1;

        *leaf_id = BPT_ID(root);
        *pos = 0;

        return BPT_ID(root);
    } else {
        *leaf_id = bpt_find_leaf(domain, root_id, key);

#if DEBUG
        {
            fprintf(stderr, "root_id:%u\tleaf_id:%u\n", root_id, *leaf_id);
        }
#endif

        root_id = bpt_place_new_key_value(domain,
                                          root_id,
                                          leaf_id,
                                          pos,
                                          key,
                                          value_id,
                                          handler);
        return root_id;
    }
}
//}}}

//{{{uint32_t bpt_insert_new_value(uint32_t root_id,
uint32_t bpt_insert_new_value(uint32_t domain,
                              uint32_t root_id,
                              uint32_t key,
                              void *value,
                              struct cache_handler *handler,
                              uint32_t *value_id,
                              uint32_t *leaf_id,
                              int *pos)
{
#if DEBUG
    fprintf(stderr, "bpt_insert_new_value\n");
#endif

    *value_id = cache.seen(domain) + 1;
    cache.add(domain,
              *value_id - 1,
              value,
              sizeof(void *),
              handler);
    return bpt_insert(domain,
                      root_id,
                      key,
                      *value_id,
                      handler,
                      leaf_id,
                      pos);
}
//}}}

//{{{ uint32_t bpt_find(uint32_t root_id,
uint32_t bpt_find(uint32_t domain,
                  uint32_t root_id,
                  uint32_t *leaf_id,
                  int *pos,
                  uint32_t key) 
{
    if (root_id == 0)
        return 0;

    *leaf_id = bpt_find_leaf(domain, root_id, key);
    // cache is zero-based, while bpt is one-based
    struct bpt_node *leaf = cache.get(domain,
                                      *leaf_id - 1,
                                      &bpt_node_cache_handler);

    int bpt_insert_key_pos = bpt_find_insert_pos(leaf, key);


    *pos = bpt_insert_key_pos;
    if ((bpt_insert_key_pos + 1) > BPT_NUM_KEYS(leaf)) 
        return 0;
    else if (key != BPT_KEYS(leaf)[bpt_insert_key_pos])
        return 0;
    else {
        return BPT_POINTERS(leaf)[bpt_insert_key_pos];
    }
}
//}}}

//{{{void bpt_write_tree(uint32_t root_id, FILE *f, char *file_name)
bool bpt_write_tree(uint32_t domain, uint32_t root_id)
{
    if (root_id == 0)
        return false;

    cache.store(domain, NULL);

    return true;
    /*
     * For each domain, start by writing out all of the non-leaf nodes.
     * Next write the leaf-nodes, where each leaf is immediantly followed by a
     * block of the pointer values. In this scheme the leaf-node will have a
     * pointer (BPT_POINTERS_BLOCK) to the block of pointer values.  
     *
     * A pointer block contains all of the pointers for a leaf node.
     *
     * only the block will be tracked in the cache, and the full block must be
     * read to access any values in the block.
     *
     * Leaf node pointers are relative to the values in a block.
     *
     * The organization of these blocks can be applicaiton specific.
     *
     * For GIGGLE, each poiner holds a list of starts (SA) and ends (SE).  All
     * SAs can be grouped into a single list, and the pointers can be offsets
     * (end positions?) into that list.  Similarly for SEs.
     *
     * 
     */

#if 0
    if (root_id == 0)
        return;
    /*
     * The nodes may not exist in the cache in a way that makes sense to store,
     * so we will do a BFS traversal of the tree and write/renumber the IDs
     * accordingly and store the offsets in an indexed array.
     */

    // Make room in the file for the id -> file offset map so we can come back
    // after and lay it down.
    
    uint32_t num_seen =  cache.seen(cache.cache);

    struct disk_store *ds = disk_store_init(num_seen + 1, &f, file_name);

    struct fifo_q *node_q = NULL, *leaf_q = NULL;
    uint32_t *id;
    id = (uint32_t *)malloc(sizeof(uint32_t));
    *id = root_id;

    struct bpt_node *to_write_node = (struct bpt_node *)
            malloc(sizeof(struct bpt_node));
    to_write_node->data = (uint32_t *)calloc(BPT_NODE_NUM_ELEMENTS,
                                             sizeof(uint32_t));

    // We run through the code an renumber the nodes so that are laid our a
    // nicely

    /*
     * Use old_id_to_new_id_os to maintain the mapping between the IDs that
     * are in memory and those that will be written to disk.
     */
    struct ordered_set *old_id_to_new_id_os =
            ordered_set_init(cache.seen(cache.cache),
                             uint_pair_sort_by_first_element_cmp,
                             uint_pair_search_by_first_element_cmp,
                             uint_pair_search_by_first_key_cmp);

    struct uint_pair *p, *r;

    // put root into a map between the current id and the on-disk id
    p = (struct uint_pair *) malloc(sizeof(struct uint_pair));
    p->first = root_id;
    p->second = old_id_to_new_id_os->num + 1;
    r = ordered_set_add(old_id_to_new_id_os, p);

    fifo_q_push(&node_q, id);

    long node_start_offset = ftell(f);

    while (fifo_q_peek(node_q) != NULL) {
        // Zero out the node that we will write to disk
        memset(to_write_node->data, 0, BPT_NODE_SIZE);

        // Get the current node's id from the queue and data from the cache
        uint32_t *curr_idp = fifo_q_pop(&node_q);
        uint32_t curr_id = *curr_idp;
        free(curr_idp);
        // cache is zero-based, while bpt is one-based
        struct bpt_node *curr_node = cache.get(domain,
                                               curr_id - 1,
                                               &bpt_node_cache_handler);

        // Get the on-disk id
        uint32_t key = curr_id;
        r = ordered_set_get(old_id_to_new_id_os, &key);
        if (r == NULL)
            errx(1, "Node %u has not been seen yet.", curr_id);

        // Populate the node that we will write to disk
        BPT_ID(to_write_node) =  r->second;
        BPT_PARENT(to_write_node) = BPT_PARENT(curr_node);
        BPT_IS_LEAF(to_write_node) = BPT_IS_LEAF(curr_node);
        BPT_LEADING(to_write_node) = BPT_LEADING(curr_node);
        BPT_NEXT(to_write_node) = BPT_NEXT(curr_node);
        BPT_NUM_KEYS(to_write_node) = BPT_NUM_KEYS(curr_node);

        uint32_t i;
        for (i = 0; i <= BPT_NUM_KEYS(curr_node); ++i)
            BPT_KEYS(to_write_node)[i] = BPT_KEYS(curr_node)[i];


        // If the node is a leaf we need to deal with the leading values
        if (BPT_IS_LEAF(curr_node)) {
            if (BPT_LEADING(curr_node) != 0) {
                // put a map between the current id and the to disk id
                p = (struct uint_pair *) malloc(sizeof(struct uint_pair));
                p->first = BPT_LEADING(curr_node);
                p->second = old_id_to_new_id_os->num + 1;
                r = ordered_set_add(old_id_to_new_id_os, p);

                if (r->second != p->second)
                    errx(1, "%u has already been seen at %u\n",
                            p->first, r->first);

                BPT_LEADING(to_write_node) =  p->second;
            }
        }

        // loop over all the pointers and if the are nodes put them on the q
        uint32_t max_pointer;
        if (BPT_IS_LEAF(curr_node)) 
            max_pointer = BPT_NUM_KEYS(curr_node) - 1;
        else 
            max_pointer = BPT_NUM_KEYS(curr_node);

        for (i = 0; i <= max_pointer; ++i) {
            if (BPT_POINTERS(curr_node)[i] != 0) {
                // put a map between the current id and the to disk id
                p = (struct uint_pair *) malloc(sizeof(struct uint_pair));
                p->first = BPT_POINTERS(curr_node)[i];
                p->second = old_id_to_new_id_os->num + 1;
                r = ordered_set_add(old_id_to_new_id_os, p);

                if (r->second != p->second)
                    errx(1, "%u has already been seen at %u\n",
                            p->first, r->first);

                // if the node is a leaf, then its pointers are to data not
                // other nodes.  that data will be handled later
                if (! BPT_IS_LEAF(curr_node)) {
                    id = (uint32_t *)malloc(sizeof(uint32_t));
                    *id = BPT_POINTERS(curr_node)[i];
                    fifo_q_push(&node_q, id);
                }

                if (r->second != p->second)
                    errx(1, "%u has already been seen at %u\n",
                            p->first, r->first);

                BPT_POINTERS(to_write_node)[i] =  p->second;
            }
        }

        // we need to loop back over the leafs so we can write out the 
        // data in the pointers and leading
        if (BPT_IS_LEAF(curr_node)) {
            id = (uint32_t *)malloc(sizeof(uint32_t));
            *id = BPT_ID(curr_node);
            fifo_q_push(&leaf_q, id);
        }

        uint32_t ret = disk_store_append(ds,
                                         to_write_node->data,
                                         BPT_NODE_SIZE);

        if (ret + 1 != BPT_ID(to_write_node))
                errx(1,
                     "Disk write is out of sync.  Saw %u.  Expected %u.",
                     ret + 1, 
                     BPT_ID(to_write_node));
    }

    long node_end_offset = ftell(f);
    long data_start_offset = node_end_offset;

    // Write out the data to disk
    while (fifo_q_peek(leaf_q) != NULL) {
        // Get the current node's id from the queue and data from the cache
        uint32_t *curr_idp = fifo_q_pop(&leaf_q);
        uint32_t curr_id = *curr_idp;
        free(curr_idp);
        struct bpt_node *curr_node = cache.get(cache.cache, curr_id);

        //Write the leading value
        if (BPT_LEADING(curr_node) != 0) {
            // Get the on-disk id
            r = ordered_set_get(old_id_to_new_id_os, &(BPT_LEADING(curr_node)));
            if (r == NULL)
                errx(1,
                     "Node %u has not been seen yet.",
                     BPT_LEADING(curr_node));
            uint32_t on_disk_id = r->second;

            // Get the data
            void *curr_pointer = cache.get(cache.cache, BPT_LEADING(curr_node));

            uint8_t *serialized_data;
            uint64_t serialized_size = serialize_leading(curr_pointer,
                                                         &serialized_data);
            uint32_t ret = disk_store_append(ds,
                                             serialized_data,
                                             serialized_size);
            free(serialized_data);

            if (ret + 1 != on_disk_id)
                errx(1,
                     "Disk write is out of sync.  Saw %u.  Expected %u.",
                     ret + 1, 
                     on_disk_id);
        }

        //Write the pointer values
        uint32_t i;
        for (i = 0; i < BPT_NUM_KEYS(curr_node); ++i) {
            if (BPT_POINTERS(curr_node)[i] != 0) {
                // Get the on-disk id
                r = ordered_set_get(old_id_to_new_id_os,
                                    &(BPT_POINTERS(curr_node)[i]));
                if (r == NULL)
                    errx(1,
                         "Node %u has not been seen yet.",
                         BPT_POINTERS(curr_node)[i]);
                uint32_t on_disk_id = r->second;

                // get the pointer from cache
                void *curr_pointer = cache.get(cache.cache,
                                               BPT_POINTERS(curr_node)[i]);
                if (curr_pointer == NULL)
                    errx(1,
                         "Pointer data %u not in cache.", 
                         BPT_POINTERS(curr_node)[i]);

                uint8_t *serialized_data;
                uint64_t serialized_size = serialize_pointer(curr_pointer,
                                                             &serialized_data);
                uint32_t ret = disk_store_append(ds,
                                                 serialized_data,
                                                 serialized_size);
                free(serialized_data);

                if (ret + 1 != on_disk_id)
                    errx(1,
                         "Disk write is out of sync.  Saw %u.  Expected %u.",
                         ret + 1, 
                         on_disk_id);
            }
        }

    }

    disk_store_sync(ds);
    free(ds->file_name);
    free(ds->offsets);
    free(ds);
    ds = NULL;

    // Move back to the end of the file before we return
    fseek(f, 0, SEEK_END);

    ordered_set_destroy(&old_id_to_new_id_os, free_wrapper);
    free(to_write_node->data);
    free(to_write_node);
#endif
    return true;
}
//}}}

//{{{void bpt_print_node(struct bpt_node *node)
void bpt_print_node(struct bpt_node *node)
{
//BPT_NODE_NUM_ELEMENTS (2*ORDER+9)
//BPT_NODE_ELEMENT_SIZE sizeof(uint32_t)
//BPT_NODE_SIZE (BPT_NODE_NUM_ELEMENTS*BPT_NODE_ELEMENT_SIZE)

    printf("BPT_ID:%u\t"
           "BPT_PARENT:%u\t"
           "BPT_IS_LEAF:%u\t"
           "BPT_LEADING:%u\t"
           "BPT_NEXT:%u\t"
           "BPT_NUM_KEYS:%u\t"
           "\n",
            BPT_ID(node),
            BPT_PARENT(node),
            BPT_IS_LEAF(node),
            BPT_LEADING(node),
            BPT_NEXT(node),
            BPT_NUM_KEYS(node));

    uint32_t i;
    for (i = 0; i < BPT_NUM_KEYS(node); ++i) {
        if (i !=0 )
            printf("\t");
        printf("%u:%u", BPT_KEYS(node)[i], BPT_POINTERS(node)[i]);
    }
    printf("\n");

//BPT_KEYS(node)
//BPT_POINTERS(node)

}
//}}}
