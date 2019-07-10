#define _GNU_SOURCE

#include <err.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <dirent.h>
#include <glob.h>
#include <sysexits.h>
#include <inttypes.h>
#include <htslib/kstring.h>

#include "bpt.h"
#include "cache.h"
#include "giggle_index.h"
#include "offset_index.h"
#include "ll.h"
#include "lists.h"
#include "file_read.h"
#include "util.h"
#include "timer.h"
#include "fastlz.h"
#include "jsw_avltree.h"

char *CHRM_INDEX_FILE_NAME = "chrm_index.dat";
char *FILE_INDEX_FILE_NAME = "file_index.dat";
char *ROOT_IDS_FILE_NAME = "root_ids.dat";
char *CACHE_FILE_NAME_PREFIX = "cache.";

//{{{ void *file_id_offset_pair_load(FILE *f, char *file_name)
void *file_id_offset_pair_load(FILE *f, char *file_name)
{
    struct file_id_offset_pair *p = (struct file_id_offset_pair*)
            malloc(sizeof( struct file_id_offset_pair));
    if (p == NULL)
        err(1, "malloc error in file_id_offset_pair_load()");

    size_t fr = fread(&(p->file_id), sizeof(uint32_t), 1, f);
    check_file_read(file_name, f, 1, fr);

    fr = fread(&(p->offset), sizeof(long), 1, f);
    check_file_read(file_name, f, 1, fr);

    return p;
}
//}}}

//{{{void file_id_offset_pair_store(void *v, FILE *f, char *file_name)
void file_id_offset_pair_store(void *v, FILE *f, char *file_name)
{
    struct file_id_offset_pair *p = (struct file_id_offset_pair*)v;

    if (fwrite(&(p->file_id), sizeof(uint32_t), 1, f) != 1)
        err(EX_IOERR,
            "Error writing file_id_offset_pair file_id '%s'.", file_name);

    if (fwrite(&(p->offset), sizeof(long), 1, f) != 1)
        err(EX_IOERR,
            "Error writing file_id_offset_pair offset '%s'.", file_name);
}
//}}}

//{{{ void *c_str_load(FILE *f, char *file_name)
void *c_str_load(FILE *f, char *file_name)
{
    uint32_t size;
    size_t fr = fread(&size, sizeof(uint32_t), 1, f);
    check_file_read(file_name, f, 1, fr);

    char *c_str = (char *)calloc(size, sizeof(char));
    if (c_str == NULL)
        err(1, "calloc error in c_str_load()");

    fr = fread(c_str, sizeof(char), size, f);
    check_file_read(file_name, f, size, fr);

    return c_str;
}
//}}}

//{{{void c_str_store(void *v, FILE *f, char *file_name)
void c_str_store(void *v, FILE *f, char *file_name)
{
    char *c_str = (char *)v;
    uint32_t size = strlen(c_str);

    if (fwrite(&size, sizeof(uint32_t), 1, f) != 1)
        err(EX_IOERR,
            "Error writing c_str size '%s'.", file_name);

    if (fwrite(c_str, sizeof(char), size, f) != size)
        err(EX_IOERR,
            "Error writing file_id_offset_pair offset '%s'.", file_name);
}
//}}}

//{{{ uint32_t giggle_insert(struct bpt_node **root,
uint32_t giggle_insert(uint32_t domain,
                       uint32_t *root_id,
                       uint32_t start,
                       uint32_t end,
                       uint64_t id)
{
#if DEBUG
    fprintf(stderr, "giggle_insert\n");
#endif
    /*
     * BP: the set of all indexing points
     * t: is a time point
     * t-: the point in BP immediatly before t
     *   the point s.t. t- < t and there does not exist a point t_m in BP such
     *   that t- < t_m < t
     * t+: the point in BP immediatly after t
     *   the point s.t. t < t+ and there does not exist a point t_m in BP such
     *   that t < t_m < t+
     * B(t): is a bucket containing pointers to all object versions whose
     *   valid time contis the interval [t_i, t_i+ - 1] 
     * SA(t): is the set of object versions whose start time is t
     * SE(t): is the set of object versions whose end time is is t -1
     *
     * t_a <- e.start
     * t_b <- e.end + 1
     *
     * search for t_a
     *
     * if not found:
     *   add t_a
     *
     * if t_a is not a leading leaf node entry:
     *   add e into SA(t_a)
     *
     * search for t_b
     *
     * if not found:
     *   add t_b
     *
     * if t_b is not a leading leaf node entry:
     *   add e into SE(t_b)
     *
     * for each leading entry t_i = t_a...t_b of a leaf node
     *   add e to B(t_i)
     *
     *
     */

#if DEBUG
    fprintf(stderr, "%u %u\n", start, end);
#endif

    uint32_t end_leaf_id = 0;
    int end_pos;

    uint32_t end_id = bpt_find(domain,
                               *root_id,
                               &end_leaf_id,
                               &end_pos,
                               end + 1);

    if (end_id == 0) {
#if DEBUG
        fprintf(stderr ,"->%u\n", end);
#endif
        void *d = giggle_data_handler.new_non_leading(domain);
        giggle_data_handler.non_leading_SE_add_scalar(domain, d, &id);

        uint32_t v_id;
        int pos;
        *root_id = bpt_insert_new_value(
                domain,
                *root_id,
                end + 1,
                d,
                &giggle_data_handler.non_leading_cache_handler,
                &v_id,
                &end_leaf_id,
                &end_pos);
    } else {
#if DEBUG
        fprintf(stderr ,"||%u\n", end);
#endif
        void *end_v = cache.get(domain,
                                end_id - 1,
                                &giggle_data_handler.non_leading_cache_handler);
        giggle_data_handler.non_leading_SE_add_scalar(domain, end_v, &id);
    }

    uint32_t start_leaf_id = 0;
    int start_pos;

    uint32_t start_id = bpt_find(domain,
                                 *root_id,
                                 &start_leaf_id,
                                 &start_pos,
                                 start);

    if (start_id == 0) {
#if DEBUG
        fprintf(stderr ,"->%u\n", start);
#endif
        void *d = giggle_data_handler.new_non_leading(domain);
        giggle_data_handler.non_leading_SA_add_scalar(domain, d, &id);
        uint32_t v_id;
        int pos;
        *root_id = bpt_insert_new_value(
                domain,
                *root_id,
                start,
                d,
                &giggle_data_handler.non_leading_cache_handler,
                &v_id,
                &start_leaf_id,
                &start_pos);
    } else {
#if DEBUG
        fprintf(stderr ,"||%u\n", start);
#endif
        void *start_v = cache.get(
                domain,
                start_id - 1,
                &giggle_data_handler.non_leading_cache_handler);
        giggle_data_handler.non_leading_SA_add_scalar(domain, start_v, &id);
    }


#if DEBUG
    fprintf(stderr, "s_id:%u e_id:%u\t", start_leaf_id, end_leaf_id);
#endif


    // For now we need to search to see which leaf the values ended up in
    // because it is possible that the leaf split on the second insert but both
    // keys ended up on the same leaf.  If they are differnet we just double
    // check to see that this is not the case.
    
    //if (start_leaf_id != end_leaf_id) 
    start_leaf_id = bpt_find_leaf(domain, *root_id, start);
    end_leaf_id = bpt_find_leaf(domain, *root_id, end + 1);

#if DEBUG
    fprintf(stderr, ":::\ts_id:%u e_id:%u\n", start_leaf_id, end_leaf_id);
#endif

#if DEBUG
    struct bpt_node *test_leaf = cache.get(domain,
                                           start_leaf_id - 1,
                                           &bpt_node_cache_handler);
    uint32_t x;
    for (x = 0; x < BPT_NUM_KEYS(test_leaf); ++x) {
        fprintf(stderr, "%u ", BPT_KEYS(test_leaf)[x]);
    }
    fprintf(stderr, "\n");
#endif

    if (start_leaf_id != end_leaf_id) {
        struct bpt_node *curr_leaf = cache.get(domain,
                                               start_leaf_id - 1,
                                               &bpt_node_cache_handler);
        do {
            curr_leaf = cache.get(domain,
                                  BPT_NEXT(curr_leaf) - 1,
                                  &bpt_node_cache_handler);

            if (BPT_LEADING(curr_leaf) == 0) {
                void *d = giggle_data_handler.new_leading(domain);
                uint32_t v_id = cache.seen(domain) + 1;
                cache.add(domain,
                          v_id - 1,
                          d,
                          sizeof(void *),
                          &giggle_data_handler.leading_cache_handler);
                giggle_data_handler.leading_B_add_scalar(domain, d, &id);
                BPT_LEADING(curr_leaf) = v_id; 
            } else {
                void *d = cache.get(domain,
                                    BPT_LEADING(curr_leaf) - 1,
                                    &giggle_data_handler.leading_cache_handler);
                giggle_data_handler.leading_B_add_scalar(domain, d, &id);
            }
        } while (BPT_ID(curr_leaf) != end_leaf_id);
    }

    return 0;
}
//}}}

//{{{ void *giggle_search(uint32_t domain,
void *giggle_search(uint32_t domain,
                    uint32_t root_id,
                    uint32_t start,
                    uint32_t end)
{
#if GIGGLE_QUERY_TRACE
    fprintf(stderr,
	    "giggle_search\tdomain:%u\tstart:%u\tend:%u\n",
	    domain,
	    start,
	    end);
#endif

    if (root_id == 0)
        return 0;

    uint32_t leaf_start_id;
    int pos_start_id;

    uint32_t nld_start_id = bpt_find(domain,
                                     root_id,
                                     &leaf_start_id, 
                                     &pos_start_id,
                                     start);
    struct bpt_node *leaf_start = cache.get(domain,
                                            leaf_start_id - 1,
                                            &bpt_node_cache_handler);
    if ((pos_start_id == 0) && (BPT_KEYS(leaf_start)[0] != start))
        pos_start_id = -1;
    else if ( (pos_start_id >=0) && 
              (pos_start_id < BPT_NUM_KEYS(leaf_start)) &&
              (BPT_KEYS(leaf_start)[pos_start_id] > start))
        pos_start_id -= 1;


#if GIGGLE_QUERY_TRACE
    fprintf(stderr,
            "leaf_start_id:%u\tpos_start_id:%d\n",
            leaf_start_id,
            pos_start_id);
#endif

    uint32_t leaf_end_id;
    int pos_end_id;

    uint32_t nld_end_id = bpt_find(domain,
                                   root_id,
                                   &leaf_end_id, 
                                   &pos_end_id,
                                   end);

    struct bpt_node *leaf_end = cache.get(domain,
                                          leaf_end_id - 1,
                                          &bpt_node_cache_handler);

    if ((pos_end_id == 0) && (BPT_KEYS(leaf_end)[0] != end))
        pos_end_id = -1;
    else if ( (pos_end_id >=0) && 
              (pos_end_id < BPT_NUM_KEYS(leaf_end)) &&
              (BPT_KEYS(leaf_end)[pos_end_id] > end))
        pos_end_id -= 1;

#if GIGGLE_QUERY_TRACE
    fprintf(stderr,
            "leaf_end_id:%u\tpos_end_id:%u\t\n",
            leaf_end_id,
            pos_end_id);
#endif

#if GIGGLE_QUERY_TRACE
    fprintf(stderr, "pos_end_id:%d %u\n", pos_end_id,
            ( ((pos_end_id >=0)&&(pos_end_id<BPT_NUM_KEYS(leaf_end))) ?
              BPT_KEYS(leaf_end)[pos_end_id] : 0)
            );
#endif

    //if ((leaf_start_id == leaf_end_id) && (pos_start_id >= pos_end_id))
    if ((leaf_start_id == leaf_end_id) && (pos_start_id > pos_end_id))
        return NULL;

#if GIGGLE_QUERY_TRACE
    if (BPT_LEADING(leaf_start) == 0)
        fprintf(stderr, "BPT_LEADING(leaf_start) == 0\n");
#endif

    return giggle_data_handler.
            giggle_collect_intersection(leaf_start_id,
                                        pos_start_id,
                                        leaf_end_id,
                                        pos_end_id,
                                        domain,
                                        NULL); 
}
//}}}

//{{{void *giggle_collect_intersection_data_in_block(uint32_t leaf_start_id,
void *giggle_collect_intersection_data_in_block(uint32_t leaf_start_id,
                                                int pos_start_id,
                                                uint32_t leaf_end_id,
                                                int pos_end_id,
                                                uint32_t domain,
                                                void **r)
{
#if DEBUG
    fprintf(stderr, "giggle_collect_intersection_data_in_block\n");
#endif

    uint32_t I_size =
            giggle_leaf_data_get_intersection_size(leaf_start_id,
                                                   pos_start_id,
                                                   leaf_end_id,
                                                   pos_end_id,
                                                   domain);

#if GIGGLE_QUERY_TRACE
    fprintf(stderr,
	    "giggle_collect_intersection_data_in_block\t"
	    "I_size:%u\n",
	    I_size);
#endif

    uint64_t *I = (uint64_t *)calloc(I_size, sizeof(uint64_t));
    if (I == NULL)
        err(1, "calloc error in giggle_collect_intersection_data_in_block()");

    struct bpt_node *leaf_start = cache.get(domain,
                                            leaf_start_id - 1,
                                            &bpt_node_cache_handler);
    struct leaf_data *leaf_start_data = 
            cache.get(domain,
                      BPT_POINTERS_BLOCK(leaf_start) - 1,
                      &leaf_data_cache_handler);

    // get everything in the leading value

    // The first step is to take the leading and the starts up to 
    // and including pos_start_id and remove ends up to and including 
    // pos_start_id
    uint32_t buff_size = leaf_start_data->num_leading 
                         + leaf_data_starts_end(leaf_start_data,
                                                leaf_start,
                                                pos_start_id) 
                         + leaf_data_ends_end(leaf_start_data,
                                              leaf_start,
                                              pos_start_id);


#if GIGGLE_QUERY_TRACE
    fprintf(stderr,
	    "giggle_collect_intersection_data_in_block\t"
	    "leaf_start_data->num_leading:%u "
	    "leaf_start_data->num_starts:%u "
	    "leaf_start_data->num_ends:%u\n",
    	    leaf_start_data->num_leading,
    	    leaf_start_data->num_starts,
    	    leaf_start_data->num_ends);
#endif

#if GIGGLE_QUERY_TRACE
    fprintf(stderr,
	    "giggle_collect_intersection_data_in_block\t"
	    "leaf_start_data->num_leading:%u "
    	    "leaf_data_starts_end:%u "
    	    "leaf_data_ends_end:%u "
    	    "buff_size:%u\n",
    	    leaf_start_data->num_leading,
    	    leaf_data_starts_end(leaf_start_data, leaf_start, pos_start_id),
    	    leaf_data_ends_end(leaf_start_data, leaf_start,  pos_start_id),
            buff_size);
#endif

    uint64_t *buff = (uint64_t *)calloc(buff_size, sizeof(uint64_t));
    if (buff == NULL)
        err(1, "calloc error in giggle_collect_intersection_data_in_block()");

    memcpy(buff,
           leaf_start_data->leading,
           leaf_start_data->num_leading * sizeof(uint64_t));

    uint32_t j;
#if GIGGLE_QUERY_TRACE
    for (j = 0; j < leaf_start_data->num_leading; ++j)
        fprintf(stderr,
                "leading\t%llu\n",
                leaf_start_data->leading[j]);
#endif

    memcpy(buff + leaf_start_data->num_leading,
           leaf_start_data->starts,
           leaf_data_starts_end(leaf_start_data,
                                leaf_start,
                                pos_start_id)*sizeof(uint64_t));

#if GIGGLE_QUERY_TRACE
    for (j = 0;
         j < leaf_data_starts_end(leaf_start_data, leaf_start, pos_start_id);
         ++j)
        fprintf(stderr,
                "starts\t%llu\n",
                leaf_start_data->starts[j]);
#endif

    memcpy(buff + leaf_start_data->num_leading + 
                leaf_data_starts_end(leaf_start_data, leaf_start, pos_start_id),
           leaf_start_data->ends,
           leaf_data_ends_end(leaf_start_data,
                              leaf_start,
                              pos_start_id)*sizeof(uint64_t));

#if GIGGLE_QUERY_TRACE
    for (j = 0;
         j < leaf_data_ends_end(leaf_start_data,leaf_start, pos_start_id);
         ++j)
        fprintf(stderr,
                "ends\t%llu\n",
                leaf_start_data->ends[j]);
#endif

    qsort(buff, buff_size, sizeof(uint64_t), uint64_t_cmp);

    uint32_t i, I_i = 0;
    for (i = 0; i < buff_size; ++i) {
        if ( ((i + 1) == buff_size) || (buff[i] != buff[i+1])) {
            I[I_i++] =  buff[i];
            if (I_i > I_size) {
                bpt_print_node(leaf_start);
                leaf_data_print(leaf_start_data);
                errx(1, "Error in giggle_collect_intersection_data_in_block. Exceeded intersection size.");
            }
        } else {
            i+=1;
        }
    }
    free(buff);

    // now process everything in between the start and end
    struct bpt_node *leaf_curr = leaf_start;
    int pos_curr_id = pos_start_id + 1;
    struct leaf_data *ld;
    uint32_t curr_size;

    // any intermediate leaves
    while (BPT_ID(leaf_curr) != leaf_end_id) {
        // do from pos_curr to the last key
        ld = cache.get(domain,
                       BPT_POINTERS_BLOCK(leaf_curr) - 1,
                       &leaf_data_cache_handler);

        curr_size = leaf_data_starts_end(ld,
                                         leaf_curr,
                                         BPT_NUM_KEYS(leaf_curr) - 1) 
                                         - leaf_data_starts_start(ld,
                                                                  leaf_curr,
                                                                  pos_curr_id);
        memcpy(I + I_i,
               ld->starts + leaf_data_starts_start(ld, leaf_curr, pos_curr_id),
               curr_size * sizeof(uint64_t));

        I_i += curr_size;

        leaf_curr = cache.get(domain,
                              BPT_NEXT(leaf_curr) - 1,
                              &bpt_node_cache_handler);
        pos_curr_id = 0;
    }

    ld = cache.get(domain,
                   BPT_POINTERS_BLOCK(leaf_curr) - 1,
                   &leaf_data_cache_handler);


    curr_size = leaf_data_starts_end(ld, leaf_curr, pos_end_id) - 
                leaf_data_starts_start(ld, leaf_curr, pos_curr_id);

    memcpy(I + I_i,
           ld->starts + leaf_data_starts_start(ld, leaf_curr, pos_curr_id),
           curr_size * sizeof(uint64_t));

    struct leaf_data_result *ldr = (struct leaf_data_result *)
            calloc(1,sizeof(struct leaf_data_result));
    if (ldr == NULL)
        err(1, "calloc error in giggle_collect_intersection_data_in_block()");

    ldr->len = I_size;
    ldr->data = I;
    ldr->next = NULL;

#if GIGGLE_QUERY_TRACE
    for (i = 0; i < I_i; ++i)
        fprintf(stderr,
	        "giggle_collect_intersection_data_in_block\t"
                "i:%u I[i]:%llu\n",
                i, I[i]);
#endif

    if ((r == NULL) || (*r == NULL)) {
        return ldr;
    } else {
        ((struct leaf_data_result *)*r)->next = ldr;
        //*r = ldr;
    }

    return ldr;
}
//}}}

//{{{uint32_t giggle_leaf_data_get_intersection_size(uint32_t leaf_start_id,
uint32_t giggle_leaf_data_get_intersection_size(uint32_t leaf_start_id,
                                                int pos_start_id,
                                                uint32_t leaf_end_id,
                                                int pos_end_id,
                                                uint32_t domain)
{
#if DEBUG
    fprintf(stderr, "giggle_leaf_data_get_intersection_size\t"
            "leaf_start_id:%u\t"
            "pos_start_id:%d\t"
            "leaf_end_id:%u\t"
            "pos_end_id:%d\t"
            "domain:%u\n",
            leaf_start_id,
            pos_start_id,
            leaf_end_id,
            pos_end_id,
            domain);
#endif

    struct bpt_node *leaf_start = cache.get(domain,
                                            leaf_start_id - 1,
                                            &bpt_node_cache_handler);
    struct leaf_data *leaf_start_data = 
            cache.get(domain,
                      BPT_POINTERS_BLOCK(leaf_start) - 1,
                      &leaf_data_cache_handler);

    // Find sizes
    // get everything in the leading value
    uint32_t i_size = leaf_start_data->num_leading;
#if DEBUG
    fprintf(stderr,
            "leaf_start_data->num_leading:%llu\t"
            "->num_starts:%llu\t"
            "->num_ends:%llu\n",
            leaf_start_data->num_leading,
            leaf_start_data->num_starts,
            leaf_start_data->num_ends);
#endif

    uint32_t i;

    i_size += leaf_data_starts_end(leaf_start_data, leaf_start, pos_start_id);
    i_size -= leaf_data_ends_end(leaf_start_data, leaf_start,  pos_start_id);

#if DEBUG
    fprintf(stderr,
            "leaf_data_starts_end:%u\t"
            "leaf_data_ends_end:%u\n",
            leaf_data_starts_end(leaf_start_data, leaf_start, pos_start_id),
            leaf_data_ends_end(leaf_start_data, leaf_start,  pos_start_id));
            
    fprintf(stderr,
            "i_size:%u\n",
            i_size);
#endif

    // now process everything in between the start and end
    struct bpt_node *leaf_curr = leaf_start;
    struct leaf_data *leaf_curr_data = leaf_start_data;
    int pos_curr_id = pos_start_id + 1;

    // any intermediate leaves
    while (BPT_ID(leaf_curr) != leaf_end_id) {
#if DEBUG
        fprintf(stderr,
                "leaf_curr:%u\t"
                "BPT_NUM_KEYS(leaf_curr):%u\t"
                "pos_curr_id:%u\t"
                "leaf_data_starts_end:%u\t"
                "leaf_data_starts_start:%u\t%u\n",
                BPT_ID(leaf_curr),
                BPT_NUM_KEYS(leaf_curr),
                pos_curr_id,
                leaf_data_starts_end(leaf_curr_data,
                                     leaf_curr,
                                     BPT_NUM_KEYS(leaf_curr)-1),
                leaf_data_starts_start(leaf_curr_data, leaf_curr, pos_curr_id),
                leaf_data_starts_end(leaf_curr_data,
                                     leaf_curr,
                                     BPT_NUM_KEYS(leaf_curr)-1)
                - leaf_data_starts_start(leaf_curr_data,
                                         leaf_curr,
                                         pos_curr_id));
#endif

        // do from pos_curr to the last key
        i_size += leaf_data_starts_end(leaf_curr_data,
                                       leaf_curr,
                                       BPT_NUM_KEYS(leaf_curr) - 1)
                  - leaf_data_starts_start(leaf_curr_data,
                                           leaf_curr,
                                           pos_curr_id);

#if DEBUG
        fprintf(stderr,
                "i_size:%u\n",
                i_size);
#endif

        leaf_curr = cache.get(domain,
                              BPT_NEXT(leaf_curr) - 1,
                              &bpt_node_cache_handler);
        leaf_curr_data = 
            cache.get(domain,
                      BPT_POINTERS_BLOCK(leaf_curr) - 1,
                      &leaf_data_cache_handler);

        pos_curr_id = 0;
    }

    i_size += leaf_data_starts_end(leaf_curr_data, leaf_curr, pos_end_id)
              - leaf_data_starts_start(leaf_curr_data, leaf_curr, pos_curr_id);

    return i_size;
}
//}}}

//{{{void *giggle_collect_intersection_data_in_pointers(uint32_t
void *giggle_collect_intersection_data_in_pointers(uint32_t leaf_start_id,
                                                   int pos_start_id,
                                                   uint32_t leaf_end_id,
                                                   int pos_end_id,
                                                   uint32_t domain,
                                                   void **_r)
{
    void *r = NULL;
#if DEBUG
    fprintf(stderr, "giggle_collect_intersection_data_in_pointers\n");
#endif

    struct bpt_node *leaf_start = cache.get(domain,
                                            leaf_start_id - 1,
                                            &bpt_node_cache_handler);

    // get everything in the leading value
    if (BPT_LEADING(leaf_start) != 0) {
        void *ld = cache.get(domain,
                             BPT_LEADING(leaf_start) - 1,
                             &giggle_data_handler.leading_cache_handler);

        giggle_data_handler.leading_union_with_B(domain, &r, ld);
    }

    // add any SA and remove any that are an SE up to and including this point
    int i;
    for (i = 0; (i < BPT_NUM_KEYS(leaf_start)) && (i <= pos_start_id); ++i) {
#if DEBUG
        fprintf(stderr,
                "BPT_KEY(leaf_start)[%u] == %u\n",
                i,
                BPT_KEYS(leaf_start)[i]);
#endif
        void *nld = cache.get(domain,
                              BPT_POINTERS(leaf_start)[i] - 1,
                              &giggle_data_handler.non_leading_cache_handler);
        giggle_data_handler.
                non_leading_union_with_SA_subtract_SE(domain,&r, nld);
    }

    // now process everything in between the start and end
    struct bpt_node *leaf_curr = leaf_start;
    int pos_curr_id = pos_start_id + 1;

    // any intermediate leaves
    while (BPT_ID(leaf_curr) != leaf_end_id) {
        // do from pos_curr to the last key
        for (i = pos_curr_id; i < BPT_NUM_KEYS(leaf_curr); ++i) {
            void *nld = cache.get(
                    domain,
                    BPT_POINTERS(leaf_curr)[i] - 1,
                    &giggle_data_handler.non_leading_cache_handler);
            giggle_data_handler.non_leading_union_with_SA(domain, &r, nld);
        }

        leaf_curr = cache.get(domain,
                              BPT_NEXT(leaf_curr) - 1,
                              &bpt_node_cache_handler);
        pos_curr_id = 0;
    }

    if (BPT_ID(leaf_curr) == leaf_end_id) {
        // add all SA's from here to either the end point
        for ( i = pos_curr_id;
             (i < BPT_NUM_KEYS(leaf_curr)) && (i <= pos_end_id); 
              ++i) {
            void *nld = cache.get(
                    domain,
                    BPT_POINTERS(leaf_curr)[i] - 1,
                    &giggle_data_handler.non_leading_cache_handler);
            giggle_data_handler.non_leading_union_with_SA(domain, &r, nld);
        }
    }

    return r;
}
//}}}

//{{{struct giggle_index *giggle_init_index(uint32_t init_size);
struct giggle_index *giggle_init_index(uint32_t init_size,
                                       char *offset_file_name)
{
    struct giggle_index *gi = (struct giggle_index *)
            malloc(sizeof(struct giggle_index));
    if (gi == NULL)
        err(1, "malloc error in giggle_init_index()");

    gi->data_dir = NULL;
    gi->len = init_size;
    gi->num = 0;
    gi->root_ids = (uint32_t *)calloc(sizeof(uint32_t), gi->len);
    if (gi->root_ids == NULL)
        err(1, "calloc error in giggle_init_index()");

    gi->chrm_idx = chrm_index_init(init_size, NULL);
    gi->file_idx = file_index_init(3, NULL);

    gi->offset_idx = offset_index_init(1000, offset_file_name);
    gi->root_ids_file_name = NULL;

    return gi;
}
//}}}

//{{{int giggle_get_chrm_id(struct giggle_index *gi, char *chrm)
uint32_t giggle_get_chrm_id(struct giggle_index *gi, char *chrm)
{
    struct str_uint_pair *r = chrm_index_get(gi->chrm_idx, chrm);

    if (r == NULL) {
        uint32_t id = chrm_index_add(gi->chrm_idx, chrm);
        uint32_t size = gi->chrm_idx->index->num;

        gi->num += 1;

        if (gi->len < size) {
            gi->root_ids = realloc(gi->root_ids, (size)*sizeof(uint32_t));
            if (gi->root_ids == NULL)
                err(1, "realloc error in giggle_get_chrm_id()");
            gi->len = size;
            uint32_t i;
            for (i = gi->num; i < gi->len; ++i)
                gi->root_ids[i] = 0;
        }

        return id;
    }

    return r->uint;
}
//}}}

//{{{void giggle_index_destroy(struct giggle_index **gi)
void giggle_index_destroy(struct giggle_index **gi)
{
    if ((*gi)->root_ids_file_name != NULL)
        free((*gi)->root_ids_file_name);
    if ((*gi)->data_dir != NULL)
        free((*gi)->data_dir);

    free((*gi)->root_ids);

    file_index_destroy(&((*gi)->file_idx));

    offset_index_destroy(&((*gi)->offset_idx));

    chrm_index_destroy(&((*gi)->chrm_idx));

    free(*gi);
    *gi = NULL;
}
//}}}

//{{{uint32_t giggle_index_file(struct giggle_index *gi,
uint32_t giggle_index_file(struct giggle_index *gi,
                           char *file_name)
{
    //fprintf(stderr, "%s\n", file_name);
    struct input_file *i = input_file_init(file_name);
    if (i == NULL)
        errx(1, "Could not open %s.\n", file_name);

    int chrm_len = 10;
    char *chrm = (char *)malloc(chrm_len*sizeof(char));
    if (chrm == NULL)
        err(1, "realloc error in giggle_index_file()");

    uint32_t start, end;
    long offset;
    kstring_t line = {0, 0, NULL};

    uint32_t file_id = file_index_add(gi->file_idx, file_name);
    struct file_data *fd = file_index_get(gi->file_idx, file_id);

    uint32_t j = 0;

    struct file_id_offset_pair *p;
    uint64_t intrv_id;

    while (i->input_file_get_next_interval(i,
                                           &chrm,
                                           &chrm_len,
                                           &start,
                                           &end,
                                           &offset,
                                           &line) >= 0) {
        //fprintf(stderr, "%s %u %u\n", chrm, start, end); 
        intrv_id = offset_index_add(gi->offset_idx,
                                    offset,
                                    &line,
                                    file_id);

        uint32_t chrm_id = giggle_get_chrm_id(gi, chrm);
        uint32_t r = giggle_insert(chrm_id,
                                   &(gi->root_ids[chrm_id]),
                                   start,
                                   end,
                                   intrv_id);
        fd->mean_interval_size += end-start;
        fd->num_intervals += 1;
        j += 1;
    }
    fd->mean_interval_size = fd->mean_interval_size/fd->num_intervals;

    if (line.s != NULL)
        free(line.s);

    input_file_destroy(&i);
    free(chrm);
    return(j);
}
//}}}

//{{{void giggle_query_region(struct giggle_index *gi,
void *giggle_query_region(struct giggle_index *gi,
                          char *chrm,
                          uint32_t start,
                          uint32_t end)
{
    uint32_t off = 0;
    if (strncmp("chr", chrm, 3) == 0)
        off = 3;

    uint32_t chr_id = giggle_get_chrm_id(gi, chrm+off);
    return giggle_search(chr_id, gi->root_ids[chr_id], start, end);
}
//}}}

//{{{uint32_t giggle_index_directory(struct giggle_index *gi,
uint32_t giggle_index_directory(struct giggle_index *gi,
                                char *path_name,
                                int verbose)
{
    glob_t results;
    int ret = glob(path_name, 0, NULL, &results);
    if (ret != 0)
        fprintf(stderr,
                "Problem with %s (%s), stopping early\n",
                path_name,
                /* ugly: */ (ret == GLOB_ABORTED ? "filesystem problem" :
                ret == GLOB_NOMATCH ? "no match of pattern" :
                ret == GLOB_NOSPACE ? "no dynamic memory" :
                "unknown problem"));

    int i;
    uint32_t total = 0;
    for (i = 0; i < results.gl_pathc; i++) {
        if (verbose)
            fprintf(stderr, "%s\n", results.gl_pathv[i]);
        total += giggle_index_file(gi, results.gl_pathv[i]);
    }

    globfree(&results);

    return total;
}
//}}}

//{{{struct giggle_index *giggle_init(uint32_t num_chrms)
struct giggle_index *giggle_init(uint32_t num_chrms,
                                 char *data_dir,
                                 uint32_t force,
                                 void (*giggle_set_data_handler)(void))
{
    if (data_dir == NULL)
        errx(1,"giggle_init: data_dir cannot be NULL.");

    char **cache_names = NULL;

    struct stat st = {0};
    if (stat(data_dir, &st) == -1) {
        mkdir(data_dir, 0700);
    } else if (force == 1) {
        rmrf(data_dir);
        mkdir(data_dir, 0700);
    } else {
        fprintf(stderr,
                "The directory '%s' already exists. "
                "Use the force option to overwrite.\n",
                data_dir);
        return NULL;
    }

    cache_names = (char **)calloc(num_chrms, sizeof(char *));
    if (cache_names == NULL)
        err(1, "calloc error in giggle_init()");

    uint32_t i, ret;
    for (i = 0; i < num_chrms; ++i) {
        ret = asprintf(&(cache_names[i]),
                       "%s/%s%u",
                       data_dir,
                       CACHE_FILE_NAME_PREFIX,
                       i);
    }

    char *offset_idx_name = NULL;
    ret = asprintf(&offset_idx_name,
                   "%s/%s",
                   data_dir,
                   OFFSET_INDEX_FILE_NAME);

    struct giggle_index *gi = giggle_init_index(num_chrms, offset_idx_name);
    gi->data_dir = strdup(data_dir);

    struct simple_cache *sc = simple_cache_init(1000,
                                                num_chrms,
                                                cache_names);

    //uint64_t_ll_giggle_set_data_handler();
    giggle_set_data_handler();

    ret = asprintf(&(gi->chrm_idx->file_name),
                   "%s/%s",
                   data_dir,
                   CHRM_INDEX_FILE_NAME);

    ret = asprintf(&(gi->file_idx->file_name),
                   "%s/%s",
                   data_dir,
                   FILE_INDEX_FILE_NAME);

    ret = asprintf(&(gi->root_ids_file_name),
                   "%s/%s",
                   data_dir,
                   ROOT_IDS_FILE_NAME);

    for (i = 0; i < num_chrms; ++i)
        free(cache_names[i]);
    free(cache_names);

    free(offset_idx_name);

    return gi;
}
//}}}

//{{{ uint32_t giggle_store(struct giggle_index *gi)
uint32_t giggle_store(struct giggle_index *gi)
{
    if (gi->chrm_idx->file_name == NULL)
        return 1;

    uint32_t i;

    giggle_data_handler.write_tree(gi);

    FILE *f = fopen(gi->root_ids_file_name, "wb");

    if (fwrite(&(gi->len), sizeof(uint32_t), 1, f) != 1)
        err(EX_IOERR, "Error writing len for root_ids'%s'.",
            gi->root_ids_file_name);

    if (fwrite(&(gi->num), sizeof(uint32_t), 1, f) != 1)
        err(EX_IOERR, "Error writing num for root_ids'%s'.",
            gi->root_ids_file_name);

    if (fwrite(gi->root_ids, sizeof(uint32_t), gi->len, f) != gi->len)
        err(EX_IOERR, "Error writing root_ids '%s'.",
            gi->root_ids_file_name);
    fclose(f);

    chrm_index_store(gi->chrm_idx);

    file_index_store(gi->file_idx);

    offset_index_store(gi->offset_idx);
    return 0;
}
//}}}

//{{{struct giggle_index *giggle_load(char *data_dir,
struct giggle_index *giggle_load(char *data_dir,
                                 void (*giggle_set_data_handler)(void))
{
    //ORDER = 100;

    if (data_dir == NULL)
        return NULL;

    struct stat st = {0};
    if (stat(data_dir, &st) == -1)
        return NULL;

    struct giggle_index *gi = (struct giggle_index *)
            malloc(sizeof(struct giggle_index));
    if (gi == NULL)
        err(1, "calloc error in giggle_load()");

    gi->data_dir = strdup(data_dir);

    // root_ids
#ifdef TIME
    struct timeval read_root_ids = in();
#endif
    int ret = asprintf(&(gi->root_ids_file_name),
                       "%s/%s",
                       data_dir,
                       ROOT_IDS_FILE_NAME);

    FILE *f = fopen(gi->root_ids_file_name, "rb");
    if (f == NULL)
        err(1, "Could not open file '%s'.\n", gi->root_ids_file_name);

    size_t fr = fread(&(gi->len), sizeof(uint32_t), 1, f);
    check_file_read(gi->root_ids_file_name, f, 1, fr);

    fr = fread(&(gi->num), sizeof(uint32_t), 1, f);
    check_file_read(gi->root_ids_file_name, f, 1, fr);

    gi->root_ids = (uint32_t *)calloc(gi->len, sizeof(uint32_t));
    if (gi->root_ids == NULL)
        err(1, "calloc error in giggle_load()");

    fr = fread(gi->root_ids, sizeof(uint32_t), gi->len, f);
    check_file_read(gi->root_ids_file_name, f, gi->len, fr);

    fclose(f);
#ifdef TIME
    fprintf(stderr,
            "giggle_load\tread root_ids\t%lu\n",
            out(read_root_ids));
#endif

    // chrm_index
#ifdef TIME
    struct timeval read_chrm_index = in();
#endif
    char *chrm_index_file_name = NULL;

    ret = asprintf(&chrm_index_file_name,
                   "%s/%s",
                   data_dir,
                   CHRM_INDEX_FILE_NAME);

    gi->chrm_idx = chrm_index_load(chrm_index_file_name);
    free(chrm_index_file_name);
#ifdef TIME
    fprintf(stderr,
            "giggle_load\tread chrm_index\t%lu\n",
            out(read_chrm_index));
#endif

#ifdef TIME
    struct timeval read_file_index = in();
#endif

    char *file_index_file_name = NULL;
    ret = asprintf(&file_index_file_name,
                   "%s/%s",
                   data_dir,
                   FILE_INDEX_FILE_NAME);
    gi->file_idx = file_index_load(file_index_file_name);
    free(file_index_file_name);
#ifdef TIME
    fprintf(stderr,
            "giggle_load\tread file_index\t%lu\n",
            out(read_file_index));
#endif

#ifdef TIME
    struct timeval read_offset_index = in();
#endif
    char *offset_index_file_name = NULL;
    ret = asprintf(&offset_index_file_name,
                   "%s/%s",
                   data_dir,
                   OFFSET_INDEX_FILE_NAME);
    gi->offset_idx = offset_index_load(offset_index_file_name);
    free(offset_index_file_name);

#ifdef TIME
    fprintf(stderr,
            "giggle_load\tread offset_index\t%lu\n",
            out(read_offset_index));
#endif


    //start();
    char **cache_names = (char **)calloc(gi->len, sizeof(char *));
    if (cache_names == NULL)
        err(1, "calloc error in giggle_load()");

    uint32_t i;
    for (i = 0; i < gi->len; ++i) {
        ret = asprintf(&(cache_names[i]),
                          "%s/%s%u",
                          data_dir,
                          CACHE_FILE_NAME_PREFIX,
                          i);
    }

#ifdef TIME
    struct timeval load_simple_cache = in();
#endif
    struct simple_cache *sc = simple_cache_init(1000,
                                                gi->len,
                                                cache_names);
    for (i = 0; i < gi->len; ++i)
        free(cache_names[i]);
    free(cache_names);


#ifdef TIME
    fprintf(stderr,
            "giggle_load\tread load_simple_cache\t%lu\n",
            out(load_simple_cache));
#endif

    giggle_set_data_handler();

    return gi;
}
//}}}

//{{{struct giggle_query_result *giggle_query(struct giggle_index *gi,
struct giggle_query_result *giggle_query(struct giggle_index *gi,
                                        char *chrm,
                                        uint32_t start,
                                        uint32_t end,
                                        struct giggle_query_result *_gqr)
{
#if GIGGLE_QUERY_TRACE
    fprintf(stderr, "giggle_query\t%s\t%u\t%u\n", chrm, start, end);
#endif

    uint32_t off = 0;
    if (strncmp("chr", chrm, 3) == 0)
        off = 3;

    void *R = NULL;
    struct str_uint_pair *r = chrm_index_get(gi->chrm_idx, chrm + off);

#if GIGGLE_QUERY_TRACE
    fprintf(stderr,
	    "giggle_query\tstr_uint_pair->str:%s\tstr_uint_pair->uint:%u\n",
	    r->str,
	    r->uint);
#endif

    if (r != NULL) {
        uint32_t chr_id = giggle_get_chrm_id(gi, chrm + off);
#if GIGGLE_QUERY_TRACE
    fprintf(stderr,
	    "giggle_query\tchr_id:%u\n",
	    chr_id);
#endif
        // HERE R COULD BE A LIST
        R = giggle_search(chr_id,
                          gi->root_ids[chr_id],
                          start,
                          end);
    }


    uint32_t i,j;
    struct giggle_query_result *gqr;
    if (_gqr == NULL) {
        gqr = (struct giggle_query_result *) 
                malloc(sizeof(struct giggle_query_result));
        if (gqr == NULL)
            err(1, "malloc error in giggle_query()");

        gqr->gi = gi;
        gqr->num_files = gi->file_idx->index->num;
        gqr->num_hits = 0;
        gqr->offsets = (struct long_uint_ll **)
            calloc(gi->file_idx->index->num, sizeof(struct long_uint_ll *));
        if (gqr->offsets == NULL)
            err(1, "calloc error in giggle_query()");

        for (i = 0; i < gi->file_idx->index->num; ++i) {
            gqr->offsets[i] = NULL;
        }
    } else {
        gqr = _gqr;
    }

    giggle_data_handler.map_intersection_to_offset_list(gi, gqr, R);

    return gqr;
}
//}}}

//{{{void giggle_query_result_destroy(struct giggle_query_result **gqr)
void giggle_query_result_destroy(struct giggle_query_result **gqr)
{
    if (*gqr == NULL)
        return;
    uint32_t i;
    for (i = 0; i < (*gqr)->gi->file_idx->index->num; ++i) {
        long_uint_ll_free((void **)&((*gqr)->offsets[i]));
    }
    free((*gqr)->offsets);
    free(*gqr);
    *gqr = NULL;
}
//}}}

//{{{uint32_t giggle_get_query_len(struct giggle_query_result *gqr,
uint32_t giggle_get_query_len(struct giggle_query_result *gqr,
                              uint32_t file_id)
{
    if (gqr->offsets[file_id] == NULL)
        return 0;
    else 
        return gqr->offsets[file_id]->len;
}
//}}}
                      
//{{{struct giggle_query_iter *giggle_get_query_itr(struct giggle_query_result
struct giggle_query_iter *giggle_get_query_itr(struct giggle_query_result *gqr,
                                               uint32_t file_id)
{
#ifdef DEBUG
    fprintf(stderr ,"giggle_get_query_itr file_id:%u\n", file_id);
#endif

    struct giggle_query_iter *gqi = (struct giggle_query_iter *)
        malloc(sizeof(struct giggle_query_iter));
    if (gqi == NULL)
        err(1, "malloc error in giggle_get_query_itr()");

    gqi->gi = gqr->gi;
    gqi->curr = 0;
    gqi->num = 0;
    gqi->file_id = file_id;
    gqi->sorted_offsets = NULL;
    gqi->sorted_offset_id_pairs = NULL;
    gqi->ipf = NULL;

    if (gqr->offsets[file_id] == NULL)
        return gqi;

#ifdef DEBUG
    fprintf(stderr,
            "giggle_get_query_itr offsets->len:%u\n", 
            gqr->offsets[file_id]->len);
#endif

    gqi->num = gqr->offsets[file_id]->len;

    gqi->sorted_offsets = (long *)
            malloc(gqr->offsets[file_id]->len * sizeof(long));
    if (gqi->sorted_offsets == NULL)
        err(1, "malloc error in giggle_get_query_itr()");

    struct long_uint_ll_node *curr = gqr->offsets[file_id]->head;
    uint32_t i = 0;
    while (curr != NULL) {
        gqi->sorted_offsets[i++] = curr->long_val;
        curr = curr->next;
    }

    qsort(gqi->sorted_offsets,
          gqr->offsets[file_id]->len,
          sizeof(long),
          long_cmp);

    return gqi;
}
//}}}

//{{{struct giggle_query_iter *giggle_get_query_data_itr(struct
struct giggle_query_iter *giggle_get_query_data_itr(
        struct giggle_query_result *gqr,
        uint32_t file_id)
{
    struct giggle_query_iter *gqi = (struct giggle_query_iter *)
        malloc(sizeof(struct giggle_query_iter));
    if (gqi == NULL)
        err(1, "malloc error in giggle_get_query_data_itr()");

    gqi->gi = gqr->gi;
    gqi->curr = 0;
    gqi->num = 0;
    gqi->file_id = file_id;
    gqi->sorted_offsets = NULL;
    gqi->sorted_offset_id_pairs = NULL;
    gqi->ipf = NULL;

    if (gqr->offsets[file_id] == NULL)
        return gqi;

    gqi->num = gqr->offsets[file_id]->len;

    gqi->sorted_offset_id_pairs = (struct long_uint_pair *)
            malloc(gqr->offsets[file_id]->len * 
                   sizeof(struct long_uint_pair));
    if (gqi->sorted_offset_id_pairs == NULL)
        err(1, "malloc error in giggle_get_query_data_itr()");

    struct long_uint_ll_node *curr = gqr->offsets[file_id]->head;
    uint32_t i = 0;
    while (curr != NULL) {
        gqi->sorted_offset_id_pairs[i].long_val = curr->long_val;
        gqi->sorted_offset_id_pairs[i].uint_val = curr->uint_val;
        i+=1;
        curr = curr->next;
    }

    qsort(gqi->sorted_offset_id_pairs,
          gqr->offsets[file_id]->len,
          sizeof(struct long_uint_pair),
          long_uint_pair_cmp);

    return gqi;
}
//}}}

//{{{int giggle_query_next(struct giggle_query_iter *gqi,
int giggle_query_next(struct giggle_query_iter *gqi,
                      char **result)
{
#ifdef DEBUG
    fprintf(stderr,
            "giggle_query_next num:%u\tcurr:%u\n",
            gqi->num,
            gqi->curr);
#endif

    if ((gqi->num == 0) || (gqi->curr == gqi->num)) {
#ifdef DEBUG
        fprintf(stderr, "giggle_query_next -1\n");
#endif
        return -1; 
    }

    if (gqi->ipf == NULL) {
        struct file_data *fd = file_index_get(gqi->gi->file_idx,
                                              gqi->file_id); 
        gqi->ipf = input_file_init(fd->file_name);
        if (gqi->ipf == NULL)
            errx(1, "Could not open %s.\n", fd->file_name);
    }

    gqi->ipf->input_file_seek(gqi->ipf, gqi->sorted_offsets[gqi->curr]);
    gqi->ipf->input_file_get_next_line(gqi->ipf, result);

    gqi->curr += 1;

#ifdef DEBUG
    fprintf(stderr, "giggle_query_next 0\n");
#endif
    return 0;
}
//}}}

//{{{int giggle_query_next_data(struct giggle_query_iter *gqi,
int giggle_query_next_data(struct giggle_query_iter *gqi,
                           void **result)
{
    if ((gqi->num == 0) || (gqi->curr == gqi->num)) {
        return -1; 
    }

    *result = OFFSET_INDEX_DATA(gqi->gi->offset_idx, 
                               gqi->sorted_offset_id_pairs[gqi->curr].uint_val);

    gqi->curr += 1;

    return 0;
}
//}}}

//{{{void giggle_iter_destroy(struct giggle_query_iter **gqi)
void giggle_iter_destroy(struct giggle_query_iter **gqi)
{
    if ((*gqi)->ipf != NULL)
        input_file_destroy(&((*gqi)->ipf));
    if ((*gqi)->sorted_offsets != NULL)
        free((*gqi)->sorted_offsets);
    if ((*gqi)->sorted_offset_id_pairs != NULL)
        free((*gqi)->sorted_offset_id_pairs);
    free(*gqi);
    *gqi = NULL;
}
//}}}

//{{{ void giggle_write_tree_cache_dump(void *giggle_index)
void giggle_write_tree_cache_dump(void *giggle_index)
{
    struct giggle_index *gi = (struct giggle_index *)giggle_index;
    uint32_t domain;
    for (domain = 0; domain < gi->num; ++domain) {
        bpt_write_tree(domain, gi->root_ids[domain]);
    }
}
//}}}

//{{{ void giggle_write_tree_leaf_data(void *giggle_index)
void giggle_write_tree_leaf_data(void *giggle_index)
{
#if DEBUG
    fprintf(stderr, "giggle_write_tree_leaf_data\n");
#endif

    struct giggle_index *gi = (struct giggle_index *)giggle_index;
    struct simple_cache *sc = (struct simple_cache *)_cache[CACHE_NAME_SPACE];

    // we will use this node to fill in the new values for all the nodes that
    // are in cache
    struct bpt_node *to_write_node = (struct bpt_node *)
            malloc(sizeof(struct bpt_node));
    if (to_write_node == NULL)
        err(1, "malloc error in giggle_write_tree_leaf_data()");

    to_write_node->data = (uint32_t *)calloc(BPT_NODE_NUM_ELEMENTS,
                                             sizeof(uint32_t));
    if (to_write_node->data == NULL)
        err(1, "malloc error in giggle_write_tree_leaf_data()");

    uint32_t domain;
    for (domain = 0; domain < gi->num; ++domain) {
        if (sc->dss[domain] != NULL)
            errx(1, "Modifying and existing bpt is not currently supported.");

        // estimate the number of elements we want to write out, this is not
        // exactly right, but that is okay 
        uint32_t num_seen =  cache.seen(domain)/2;

        // Each new node or leaf data in this domain will be appended to this
        // disk store
        struct disk_store *ds = disk_store_init(sc->seens[domain],
                                                NULL,
                                                sc->index_file_names[domain],
                                                NULL,
                                                sc->data_file_names[domain]);


        // Use old_id_to_new_id_os to maintain the mapping between the IDs that
        // are in memory and those that will be written to disk.
        struct ordered_set *old_id_to_new_id_os =
            ordered_set_init(num_seen,
                             uint_pair_sort_by_first_element_cmp,
                             uint_pair_search_by_first_element_cmp,
                             uint_pair_search_by_first_key_cmp);

        // Start with the root node for this domain
        struct bpt_node *curr_node = cache.get(domain,
                                               gi->root_ids[domain] - 1,
                                               &bpt_node_cache_handler);

        struct uint_pair *p, *r;

        // put root into a map between the current id and the on-disk id
        // first will but current id 
        // second is the on-disk id
        p = (struct uint_pair *) malloc(sizeof(struct uint_pair));
        if (p == NULL)
            err(1, "malloc error in giggle_write_tree_leaf_data()");

        p->first = BPT_ID(curr_node);
        p->second = old_id_to_new_id_os->num + 1;
        r = ordered_set_add(old_id_to_new_id_os, p);

        uint32_t new_root_id = p->second;

        struct fifo_q *node_q = NULL, *leaf_q = NULL;
        uint32_t *id = (uint32_t *)malloc(sizeof(uint32_t));
        if (id == NULL)
            err(1, "malloc error in giggle_write_tree_leaf_data()");

        *id = BPT_ID(curr_node);
        fifo_q_push(&node_q, id);

        while (fifo_q_peek(node_q) != NULL) {
            // pop the next node off the queue and get it from cache
            uint32_t *curr_idp = fifo_q_pop(&node_q);
            uint32_t curr_id = *curr_idp;
            free(curr_idp);
            // cache is zero-based, while bpt is one-based
            curr_node = cache.get(domain,
                                  curr_id - 1,
                                  &bpt_node_cache_handler);
            // Basic steps:
            // - copy the values over to the temp node that we will write out
            // if the node is not a leaf:
            //   - map the pointer values using old_id_to_new_id_os
            //   - put the child nodes into the queue
            //   - set the pointer head to zero
            //   - put the node onto disk
            // if the node is a leaf:
            //   - get the leaf data
            //   - add the leaf data to the and the queue cache so we can write
            //     it out later
            //   - set the pointer head to the disk id of the leaf data
            //   -- each 32bit pointer will be split into 2 16bit offsets
            //      the first will be the start offset and the second the end

            // Zero out the node that we will write to disk
            memset(to_write_node->data, 0, BPT_NODE_SIZE);

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
            BPT_POINTERS_BLOCK(to_write_node) = 0;
            uint32_t i;
            for (i = 0; i <= BPT_NUM_KEYS(curr_node); ++i)
                BPT_KEYS(to_write_node)[i] = BPT_KEYS(curr_node)[i];


            if (BPT_IS_LEAF(curr_node) == 0) {

               for (i = 0; i <= BPT_NUM_KEYS(curr_node); ++i) {
                    if (BPT_POINTERS(curr_node)[i] != 0) {
                        // put a map between the current id and the to disk id
                        p = (struct uint_pair *)
                                malloc(sizeof(struct uint_pair));
                        if (p == NULL)
                            err(1,
                                "malloc error in "
                                "giggle_write_tree_leaf_data()");


                        p->first = BPT_POINTERS(curr_node)[i];
                        p->second = old_id_to_new_id_os->num + 1;
                        r = ordered_set_add(old_id_to_new_id_os, p);

                        if (r->second != p->second)
                            errx(1,
                                 "%u has already been seen at %u\n",
                                 p->first,
                                 r->first);

                        // update the node we are writing to disk with the new 
                        // id
                        BPT_POINTERS(to_write_node)[i] =  p->second;

                        // put the child on the queue
                        id = (uint32_t *)malloc(sizeof(uint32_t));
                        if (id == NULL)
                            err(1,
                                "malloc error in "
                                "giggle_write_tree_leaf_data()");

                        *id = BPT_POINTERS(curr_node)[i];
                        fifo_q_push(&node_q, id);
                    }
                }

                uint8_t *ds_node;
                uint64_t d_size = bpt_node_serialize(to_write_node,
                                                     (void **)&ds_node);

                // Write the mapped node to disk
                uint32_t ret = disk_store_append(ds,
                                                 ds_node,
                                                 d_size);

                // Make sure it gets the ID that we expect
                if (ret + 1 != BPT_ID(to_write_node))
                    errx(1,
                         "Disk write is out of sync.  Saw %u.  Expected %u.",
                         ret + 1, 
                         BPT_ID(to_write_node));

                free(ds_node);
            } else {

                // replace the next id with the on disk id
                if (BPT_NEXT(curr_node) != 0) {
                    key = BPT_NEXT(curr_node);
                    r = ordered_set_get(old_id_to_new_id_os, &key);
                    if (r == NULL)
                        errx(1, "Node %u has not been seen yet.", key);
                    BPT_NEXT(to_write_node) = r->second;
                }

                // get the leaf data, then add it to the cache so we can
                // grab it later
                struct leaf_data *lf = NULL;
                //uint16_t *starts_ends_offsets = NULL;
                uint32_t leaf_data_size = 
                        giggle_get_leaf_data(gi,
                                             domain,
                                             curr_id,
                                             &lf);
                                             //&starts_ends_offsets);
                if (leaf_data_size == 0)
                    errx(1, "Could not get leaf data.");

                //uint8_t *output = (uint8_t *)malloc(
                                   //2*leaf_data_size * sizeof(uint32_t));

                //int cs = fastlz_compress(lf->data,
                                    //leaf_data_size * sizeof(uint32_t),
                                    //output);

                uint32_t data_id = cache.seen(domain) + 1;
                cache.add(domain,
                          data_id - 1,
                          lf,
                          sizeof(struct leaf_data),
                          &leaf_data_cache_handler);

                p = (struct uint_pair *) malloc(sizeof(struct uint_pair));
                if (p == NULL)
                    err(1,
                        "malloc error in giggle_write_tree_leaf_data()");

                p->first = data_id;
                p->second = old_id_to_new_id_os->num + 1;
                r = ordered_set_add(old_id_to_new_id_os, p);

                id = (uint32_t *)malloc(sizeof(uint32_t));
                if (id == NULL)
                    err(1,
                        "malloc error in giggle_write_tree_leaf_data()");

                *id = data_id;
                fifo_q_push(&leaf_q, id);

                BPT_POINTERS_BLOCK(to_write_node) = p->second;

                uint8_t *ds_node;
                uint64_t d_size = bpt_node_serialize(to_write_node,
                                                     (void **)&ds_node);

                // Write the mapped node to disk
                uint32_t ret = disk_store_append(ds,
                                                 ds_node,
                                                 d_size);

                // Make sure it gets the ID that we expect
                if (ret + 1 != BPT_ID(to_write_node))
                    errx(1,
                         "Disk write is out of sync.  Saw %u.  Expected %u.",
                         ret + 1, 
                         BPT_ID(to_write_node));

                free(ds_node);
            }
        }
        while (fifo_q_peek(leaf_q) != NULL) {
            // pop the next node off the queue and get it from cache
            uint32_t *curr_idp = fifo_q_pop(&leaf_q);
            uint32_t curr_id = *curr_idp;
            free(curr_idp);
            // cache is zero-based, while bpt is one-based
            struct leaf_data *lf = cache.get(domain,
                                             curr_id - 1,
                                             &leaf_data_cache_handler);

            uint8_t *ds_data;
            uint64_t d_size = leaf_data_serialize(lf,
                                                 (void **)&ds_data);

            // Write the mapped node to disk
            uint32_t ret = disk_store_append(ds,
                                             ds_data,
                                             d_size);
            free(ds_data);
        }

        sc->dss[domain] = ds;
        ordered_set_destroy(&old_id_to_new_id_os, free_wrapper);
        gi->root_ids[domain] = new_root_id;
    }

    bpt_node_free_mem((void **)&to_write_node);
}
//}}}

//{{{ leaf_data_cache_handler

    /*
     * struct leaf_data {
     *   uint32_t num_leading, num_starts, num_ends;
     *   uint32_t *leading, *starts, *ends, *data;
     * };
     */
struct cache_handler leaf_data_cache_handler = {leaf_data_serialize, 
                                                leaf_data_deserialize,
                                                leaf_data_free_mem};


//}}}

//{{{ uint32_t giggle_get_leaf_data(struct giggle_index *gi,
uint32_t giggle_get_leaf_data(struct giggle_index *gi,
                              uint32_t domain,
                              uint32_t leaf_id,
                              struct leaf_data **lf)
{
    // cache is zero-based, while bpt is one-based
    struct bpt_node *curr_node = cache.get(domain,
                                           leaf_id - 1,
                                           &bpt_node_cache_handler);

    // If the node is a leaf we need to deal with the leading values
    if (BPT_IS_LEAF(curr_node)) {
        *lf = (struct leaf_data *) calloc(1, sizeof(struct leaf_data));
        if (*lf == NULL)
            err(1, "calloc error in giggle_get_leaf_data()");

        // Do one scan to find the sizes
        
        // Get the number of leading values
        if (BPT_LEADING(curr_node) != 0) {
            struct uint64_t_ll_bpt_leading_data *ld = 
                    cache.get(domain,
                              BPT_LEADING(curr_node) - 1,
                              &uint64_t_ll_leading_cache_handler);
            (*lf)->num_leading = ld->B->len;
        }

        uint32_t j, k;

        // Get the number of starts and ends
        for (j = 0; j <= BPT_NUM_KEYS(curr_node) - 1; ++j) {
            struct uint64_t_ll_bpt_non_leading_data *nld = 
                    cache.get(domain,
                              BPT_POINTERS(curr_node)[j] - 1,
                              &uint64_t_ll_non_leading_cache_handler);

            (*lf)->num_starts += (nld->SA == NULL ? 0 : nld->SA->len);
            (*lf)->num_ends += (nld->SE == NULL ? 0 : nld->SE->len);
        }

        // Allocate the memory and put in helper pointers
        (*lf)->data = (uint64_t *)calloc(
                (*lf)->num_leading + (*lf)->num_starts + (*lf)->num_ends,
                LEAF_LEADING_STARTS_ENDS_SIZE);

        if ((*lf)->data == NULL)
            err(1, "calloc error in giggle_get_leaf_data()");

        if ((*lf)->num_leading > 0)
            (*lf)->leading = (*lf)->data;
        else
            (*lf)->leading = NULL;

        if ((*lf)->num_starts > 0)
            (*lf)->starts = (*lf)->data + (*lf)->num_leading;
        else
            (*lf)->starts = NULL;

        if ((*lf)->num_ends > 0)
            (*lf)->ends = (*lf)->data + (*lf)->num_leading + (*lf)->num_starts;
        else
            (*lf)->ends = NULL;

        (*lf)->starts_pointers = (uint32_t *)calloc(ORDER, LEAF_POINTERS_SIZE);
        if ((*lf)->starts_pointers == NULL)
            err(1, "calloc error in giggle_get_leaf_data().\n");

        (*lf)->ends_pointers = (uint32_t *)calloc(ORDER, LEAF_POINTERS_SIZE);
        if ((*lf)->ends_pointers == NULL)
            err(1, "calloc error in giggle_get_leaf_data().\n");

        // Do a second scan to get the data into the array
        uint32_t leading_i = 0, starts_i = 0, ends_i = 0;

        // Copy data into leading buffer
        if (BPT_LEADING(curr_node) != 0) {
            struct uint64_t_ll_bpt_leading_data *ld = 
                    cache.get(domain,
                              BPT_LEADING(curr_node) - 1,
                              &uint64_t_ll_leading_cache_handler);
            k = leading_i;
            if (ld->B->len > 0) {
                struct uint64_t_ll_node *curr = ld->B->head;
                while (curr != NULL) {
                    (*lf)->leading[k] = curr->val;
                    k += 1;
                    curr = curr->next;
                }
                qsort((*lf)->leading + leading_i,
                      ld->B->len,
                      sizeof(uint64_t),
                      uint64_t_cmp);
            }
            leading_i = k;
        }

        // Copy data into starts and ends buffer
        for (j = 0; j <= BPT_NUM_KEYS(curr_node) - 1; ++j) {
            struct uint64_t_ll_bpt_non_leading_data *nld = 
                    cache.get(domain,
                              BPT_POINTERS(curr_node)[j] - 1,
                              &uint64_t_ll_non_leading_cache_handler);

            k = starts_i;
            if ((nld->SA != NULL) && (nld->SA->len > 0)) {
                struct uint64_t_ll_node *curr = nld->SA->head;
                while (curr != NULL) {
                    (*lf)->starts[k] = curr->val;
                    k += 1;
                    curr = curr->next;
                }
                qsort((*lf)->starts + starts_i,
                      nld->SA->len,
                      sizeof(uint64_t),
                      uint64_t_cmp);

                (*lf)->starts_pointers[j] = 
                    (j == 0) ? nld->SA->len 
                             : (*lf)->starts_pointers[j-1] + nld->SA->len;
                
            } else {
                (*lf)->starts_pointers[j] = 
                    (j == 0) ? 0 : (*lf)->starts_pointers[j-1];
            }

            starts_i = k;

            k = ends_i;
            if ((nld->SE != NULL) && (nld->SE->len > 0)) {
                struct uint64_t_ll_node *curr = nld->SE->head;
                while (curr != NULL) {
                    (*lf)->ends[k] = curr->val;
                    k += 1;
                    curr = curr->next;
                }
                qsort((*lf)->ends + ends_i,
                      nld->SE->len,
                      sizeof(uint64_t),
                      uint64_t_cmp);

                (*lf)->ends_pointers[j] = 
                    (j == 0) ? nld->SE->len 
                             : (*lf)->ends_pointers[j-1] + nld->SE->len;
            } else {
                (*lf)->ends_pointers[j] = 
                    (j == 0) ? 0 : (*lf)->ends_pointers[j-1];
            }

            ends_i = k;
        }
        return (*lf)->num_leading + (*lf)->num_starts + (*lf)->num_ends;
    } else {
        return 0;
    }
}
//}}}

//{{{void leaf_data_map_intersection_to_offset_list(struct giggle_index *gi,
void leaf_data_map_intersection_to_offset_list(struct giggle_index *gi,
                                            struct giggle_query_result *gqr,
                                            void *_R)
{
#ifdef DEBUG
    fprintf(stderr,
            "leaf_data_map_intersection_to_offset_list\n");
#endif
    struct leaf_data_result *R = (struct leaf_data_result *)_R;
    /*
    struct leaf_data_result {
        uint32_t len;
        uint32_t *data;
        struct leaf_data_result *next;
    };
    */

    if (R != NULL) {
        gqr->num_hits += R->len;

#ifdef GIGGLE_QUERY_TRACE
        fprintf(stderr,
                "leaf_data_map_intersection_to_offset_list\t"
                "R->len:%u\n",
                R->len);
#endif

        uint32_t i;
        for (i = 0; i < R->len; ++i) {
            struct file_id_offset_pair fid_off = 
                    offset_index_get(gi->offset_idx, R->data[i]);
            //long_ll_append(&(gqr->offsets[fid_off.file_id]),fid_off.offset);
            long_uint_ll_append(&(gqr->offsets[fid_off.file_id]),
                                fid_off.offset,
                                R->data[i]);
        }

        struct leaf_data_result *tmp_R = R->next;
        free(R->data);
        free(R);
        R = tmp_R;
    } 
}
//}}}

//{{{ chrm_index
//{{{struct chrm_index *chrm_index_init(uint32_t init_size,
struct chrm_index *chrm_index_init(uint32_t init_size,
                                   char *file_name)
{
    struct chrm_index *idx = 
        (struct chrm_index *)malloc(sizeof(struct chrm_index));
    if (idx == NULL)
        err(1, "malloc error in chrm_index_init()");

    idx->index = ordered_set_init(init_size,
                                  str_uint_pair_sort_element_cmp,
                                  str_uint_pair_search_element_cmp,
                                  str_uint_pair_search_key_cmp);
    if (file_name != NULL)
        idx->file_name = strdup(file_name);
    else
        idx->file_name = NULL;

    return idx;
}
//}}}

//{{{struct str_uint_pair *chrm_index_get(struct chrm_index *ci,
struct str_uint_pair *chrm_index_get(struct chrm_index *ci,
                                    char *chrm)
{
    return (struct str_uint_pair *) ordered_set_get(ci->index, chrm);
}
//}}}

//{{{uint32_t chrm_index_add(struct chrm_index *ci,
uint32_t chrm_index_add(struct chrm_index *ci,
                        char *chrm)
{
    struct str_uint_pair *p = (struct str_uint_pair *)
                malloc(sizeof(struct str_uint_pair));
    if (p == NULL)
        err(1, "malloc error in chrm_index_add()");

    p->uint = ci->index->num;
    p->str = strdup(chrm);

    struct str_uint_pair *r = (struct str_uint_pair *)
            ordered_set_add(ci->index, p);

    return r->uint;
}
//}}}

//{{{void chrm_index_destroy(struct chrm_index **ci)
void chrm_index_destroy(struct chrm_index **ci)
{
    if ((*ci)->file_name != NULL)
        free((*ci)->file_name);
    ordered_set_destroy(&((*ci)->index), str_uint_pair_free);
    free(*ci);
    *ci = NULL;
}
//}}}

//{{{struct chrm_index *chrm_index_load(char *file_name)
struct chrm_index *chrm_index_load(char *file_name)
{
    struct chrm_index *idx = 
        (struct chrm_index *)malloc(sizeof(struct chrm_index));
    if (idx == NULL)
        err(1, "malloc error in chrm_index_load()");

    idx->file_name = strdup(file_name);

    FILE *f = fopen(idx->file_name, "rb");
    idx->index = ordered_set_load(f,
                                  idx->file_name,
                                  str_uint_pair_load,
                                  str_uint_pair_sort_element_cmp,
                                  str_uint_pair_search_element_cmp,
                                  str_uint_pair_search_key_cmp);
    fclose(f);

    return idx;
}
//}}}

//{{{void chrm_index_store(struct chrm_index *ci)
void chrm_index_store(struct chrm_index *ci)
{
    if (ci->file_name == NULL)
        errx(1,"No output file given for chrm_index.");

    FILE *f = fopen(ci->file_name, "wb");
    ordered_set_store(ci->index,
                      f,
                      ci->file_name,
                      str_uint_pair_store);
    fclose(f);
}
//}}}
//}}}

//{{{ file_index
//{{{struct file_index *file_index_init(uint32_t init_size, char *file_name)
struct file_index *file_index_init(uint32_t init_size, char *file_name)
{
    struct file_index *fi = (struct file_index *)
            malloc(sizeof(struct file_index));
    if (fi == NULL)
        err(1, "malloc error in chrm_index_init()");

    fi->index = unordered_list_init(3);
    fi->file_name = NULL;
    if (file_name != NULL)
        fi->file_name = strdup(file_name);
    return fi;
}
//}}}

//{{{ void file_index_destroy(struct file_index **fi)
void file_index_destroy(struct file_index **fi)
{
    unordered_list_destroy(&((*fi)->index), file_data_free);
    if ((*fi)->file_name != NULL) {
        free((*fi)->file_name);
        (*fi)->file_name = NULL;
    }
    free(*fi);
    *fi = NULL;
}
//}}}

//{{{uint32_t file_index_add(struct file_index *fi, char *file_name)
uint32_t file_index_add(struct file_index *fi, char *file_name)
{
    struct file_data *fd = (struct file_data *)
            calloc(1, sizeof(struct file_data));
    if (fd == NULL)
        err(1, "calloc error in chrm_index_add()");

    fd->file_name = strdup(file_name);
    return unordered_list_add(fi->index, fd);
}
//}}}

//{{{void file_index_store(struct file_index *fi)
void file_index_store(struct file_index *fi)
{
    if (fi->file_name == NULL)
        errx(1,"No output file given for file_index.");

    FILE *f = fopen(fi->file_name, "wb");
    unordered_list_store(fi->index,
                         f,
                         fi->file_name,
                         file_data_store);
    fclose(f);
}
//}}}

//{{{struct file_index *file_index_load(char *file_name)
struct file_index *file_index_load(char *file_name)
{
    struct file_index *fi = (struct file_index *)
        malloc(sizeof(struct file_index));
    if (fi == NULL)
        err(1, "malloc error in chrm_index_load()");

    fi->file_name = strdup(file_name);
    FILE *f = fopen(file_name, "rb");
    fi->index = unordered_list_load(f,
                                    fi->file_name,
                                    file_data_load);
    fclose(f);

    return fi;
}
//}}}

//{{{struct file_data *file_index_get(struct file_index *fi, uint32_t id)
struct file_data *file_index_get(struct file_index *fi, uint32_t id)
{
    return (struct file_data *)unordered_list_get(fi->index, id);
}
//}}}
//}}}

//{{{ uint64_t giggle_bulk_insert(char *input_path_name,
uint64_t giggle_bulk_insert(char *input_path_name,
                            char *output_path_name,
                            uint32_t force)
{
    // Make the output directory
    struct stat st = {0};
    if (stat(output_path_name, &st) == -1) {
        mkdir(output_path_name, 0700);
    } else if (force == 1) {
        rmrf(output_path_name);
        mkdir(output_path_name, 0700);
    } else {
        fprintf(stderr,
                "The directory '%s' already exists. "
                "Use the force option to overwrite.\n",
                output_path_name);
        return 0;
    }

    struct giggle_index *gi = (struct giggle_index *)
            malloc(sizeof(struct giggle_index));
    if (gi == NULL)
        err(1, "malloc error in giggle_bulk_insert()");

    gi->data_dir = strdup(output_path_name);

    // open files
    gi->file_idx = NULL;
    struct input_file **i_files = NULL;
    uint32_t num_input_files = giggle_bulk_insert_open_files(input_path_name,
                                                             gi->data_dir, 
                                                             &i_files,
                                                             &(gi->file_idx));

    //init offset index
    char *offset_index_file_name = NULL;
    int ret = asprintf(&offset_index_file_name,
                       "%s/%s",
                       gi->data_dir,
                       OFFSET_INDEX_FILE_NAME);
    gi->offset_idx = offset_index_init(1000, offset_index_file_name);
    free(offset_index_file_name);
 
    //init chrm index
    char *chrm_index_file_name = NULL;
    ret = asprintf(&chrm_index_file_name,
                   "%s/%s",
                   gi->data_dir,
                   CHRM_INDEX_FILE_NAME);
    gi->chrm_idx = chrm_index_init(24, chrm_index_file_name);
    free(chrm_index_file_name);
 
    // prime pqs
    pri_queue pq_start = priq_new(num_input_files);
    // Since we know that each file will have at most one start in the priority
    // queue, we can reduce mallocs by reusing the array
    struct pq_data *pqd_starts = (struct pq_data *)
            malloc(num_input_files * sizeof(struct pq_data));
    if (pqd_starts == NULL)
        err(1, "malloc error in giggle_bulk_insert()");

    pri_queue pq_end = priq_new(num_input_files);

    giggle_bulk_insert_prime_pqs(gi,
                                 &pq_start,
                                 pqd_starts,
                                 &pq_end,
                                 i_files,
                                 num_input_files);
    
    // scan files
    giggle_bulk_insert_build_leaf_levels(gi,
                                         &pq_start,
                                         pqd_starts,
                                         &pq_end,
                                         i_files,
                                         num_input_files);

    // clean up
    priq_free(pq_end);
    priq_free(pq_start);
    free(pqd_starts);
    uint32_t i;
    for (i = 0; i < num_input_files; ++i)
        input_file_destroy(&(i_files[i]));
    free(i_files);

    // build trees
    giggle_bulk_insert_build_tree_on_leaves(gi);

    for (i = 0; i < gi->file_idx->index->num; ++i) {
        struct file_data *fd = file_index_get(gi->file_idx, i);
        fd->mean_interval_size = fd->mean_interval_size / fd->num_intervals;
    }

    // save 
    FILE *f = fopen(gi->root_ids_file_name, "wb");

    if (fwrite(&(gi->len), sizeof(uint32_t), 1, f) != 1)
        err(EX_IOERR, "Error writing len for root_ids'%s'.",
            gi->root_ids_file_name);

    if (fwrite(&(gi->num), sizeof(uint32_t), 1, f) != 1)
        err(EX_IOERR, "Error writing num for root_ids'%s'.",
            gi->root_ids_file_name);

    if (fwrite(gi->root_ids, sizeof(uint32_t), gi->len, f) != gi->len)
        err(EX_IOERR, "Error writing root_ids '%s'.",
            gi->root_ids_file_name);
    fclose(f);

    chrm_index_store(gi->chrm_idx);
    file_index_store(gi->file_idx);
    offset_index_store(gi->offset_idx);

    uint64_t num_intervals = gi->offset_idx->index->num;

    if (gi->root_ids_file_name != NULL)
        free(gi->root_ids_file_name);
    if (gi->data_dir != NULL)
        free(gi->data_dir);

    free(gi->root_ids);
    file_index_destroy(&(gi->file_idx));
    offset_index_destroy(&(gi->offset_idx));
    chrm_index_destroy(&(gi->chrm_idx));
    free(gi);
    gi = NULL;

    return num_intervals;
}
//}}}

//{{{void giggle_bulk_insert_build_tree_on_leaves(struct giggle_index *gi)
void giggle_bulk_insert_build_tree_on_leaves(struct giggle_index *gi)
{
    gi->len = gi->chrm_idx->index->num;
    gi->num = gi->len;
    gi->root_ids = (uint32_t *)malloc(gi->num * sizeof(uint32_t));
    if (gi->root_ids == NULL)
        err(1, "malloc error in giggle_bulk_insert_build_tree_on_leaves()");

    int ret = asprintf(&(gi->root_ids_file_name),
                       "%s/%s",
                       gi->data_dir,
                       ROOT_IDS_FILE_NAME);
 
    char *ds_curr_index_file_name = NULL, *ds_curr_data_file_name = NULL;
    struct disk_store *curr_ds;
    uint32_t curr_chrm_id;
    for (curr_chrm_id = 0;
         curr_chrm_id < gi->chrm_idx->index->num;
         ++curr_chrm_id) {

        ret = asprintf(&ds_curr_index_file_name,
                       "%s/%s%u.idx",
                       gi->data_dir,
                       CACHE_FILE_NAME_PREFIX,
                       curr_chrm_id);
        ret = asprintf(&ds_curr_data_file_name,
                       "%s/%s%u.dat",
                       gi->data_dir,
                       CACHE_FILE_NAME_PREFIX,
                       curr_chrm_id);
        curr_ds = disk_store_load(NULL,
                                  ds_curr_index_file_name,
                                  NULL,
                                  ds_curr_data_file_name);
        free(ds_curr_index_file_name);
        free(ds_curr_data_file_name);
        ds_curr_index_file_name = NULL;
        ds_curr_data_file_name = NULL;

        // Here we need to loop over each level of the tree until the current
        // level has just one element in which case that element is the root

        uint32_t num_leaf_node_leaf_data = curr_ds->num;
        uint32_t curr_level_num_nodes = num_leaf_node_leaf_data / 2;
        uint32_t curr_level_first_id = 1;
        uint32_t curr_level_is_leaf = 1;
        uint32_t new_level_first_id = 0;
        uint32_t new_level_len = 
                giggle_bulk_insert_add_tree_level(curr_ds,
                                                  curr_level_first_id,
                                                  curr_level_num_nodes,
                                                  curr_level_is_leaf,
                                                  &new_level_first_id);

        while(new_level_len > 1) {
            curr_level_num_nodes = new_level_len;
            curr_level_first_id = new_level_first_id;
            curr_level_is_leaf = 0;
            new_level_first_id = 0;
            new_level_len = 
                    giggle_bulk_insert_add_tree_level(curr_ds,
                                                      curr_level_first_id,
                                                      curr_level_num_nodes,
                                                      curr_level_is_leaf,
                                                      &new_level_first_id);
        }

        if (new_level_len == 1) {
            gi->root_ids[curr_chrm_id] = new_level_first_id;
        } else if (new_level_len == 0) {
            gi->root_ids[curr_chrm_id] = curr_level_first_id;
        } 

        disk_store_destroy(&curr_ds);
    }
}
//}}}

//{{{void giggle_bulk_insert_build_leaf_levels(struct giggle_index *gi,
void giggle_bulk_insert_build_leaf_levels(struct giggle_index *gi,
                                          pri_queue *pq_start,
                                          struct pq_data *pqd_starts,
                                          pri_queue *pq_end,
                                          struct input_file **i_files,
                                          uint32_t num_input_files)
{
    // Grab the top element on the start pq 
    priority pri_start;
    struct pq_data *pqd_start =
            (struct pq_data *)priq_top(*pq_start, &pri_start);
 
    // curr_pos and curr_chrm track the status of the indexing
    uint32_t curr_pos = pri_start.pos;
    char curr_chrm[50];
    strcpy(curr_chrm, pri_start.chrm);

    // register the chrom with chrom index
    uint32_t curr_chrm_id = chrm_index_add(gi->chrm_idx, curr_chrm);

    //init disk store, do this at the start of every chrom
    char *ds_curr_index_file_name = NULL, *ds_curr_data_file_name = NULL;
    uint32_t ret = asprintf(&ds_curr_index_file_name,
                            "%s/%s%u.idx",
                            gi->data_dir,
                            CACHE_FILE_NAME_PREFIX,
                            curr_chrm_id);
    ret = asprintf(&ds_curr_data_file_name,
                   "%s/%s%u.dat",
                   gi->data_dir,
                   CACHE_FILE_NAME_PREFIX,
                   curr_chrm_id);
    struct disk_store *curr_ds = disk_store_init(10,
                                                 NULL,
                                                 ds_curr_index_file_name,
                                                 NULL,
                                                 ds_curr_data_file_name);
    free(ds_curr_index_file_name);
    free(ds_curr_data_file_name);

    // Collect the values into this node, then write it and clear 
    struct bpt_node *bpn = (struct bpt_node *) malloc(sizeof(struct bpt_node));
    if (bpn == NULL)
        err(1, "malloc error in giggle_bulk_insert_build_leaf_levels()");

    bpn->data = (uint32_t *) malloc(BPT_NODE_NUM_ELEMENTS  * sizeof(uint32_t));
    if (bpn->data == NULL)
        err(1, "malloc error in giggle_bulk_insert_build_leaf_levels()");

    memset(bpn->data, 0, BPT_NODE_SIZE);

    BPT_ID(bpn) = curr_ds->num + 1;//1-based
    BPT_PARENT(bpn) = 0;
    BPT_IS_LEAF(bpn) = 1;
    BPT_LEADING(bpn) = 0;
    BPT_NEXT(bpn) = 0;
    BPT_NUM_KEYS(bpn) = 0;
    BPT_POINTERS_BLOCK(bpn) = 0;

    // These will be used to create the leaf data for each node
    uint32_t num_leading = 0, num_starts = 0, num_ends = 0;
    struct uint64_t_array *leading, *starts, *ends;
    leading = uint64_t_array_init(100);
    starts = uint64_t_array_init(100);
    ends = uint64_t_array_init(100);

    struct uint32_t_array *starts_pointers, *ends_pointers;
    starts_pointers = uint32_t_array_init(100);
    ends_pointers = uint32_t_array_init(100);

    // This tree will track intervals that have begun and not yet ended and
    // will be used to populate the leading value of nodes
    jsw_avltree_t *avl = jsw_avlnew(uint64_cmp_f, uint64_dup_f, uint64_rel_f);

    // add the current possition to the node
    ret = giggle_bulk_insert_append_bpt_key(bpn,
                                            curr_pos,
                                            curr_ds,
                                            avl,
                                            leading,
                                            starts,
                                            ends,
                                            starts_pointers,
                                            ends_pointers);

    // These will be used to read intervals from files
    int chrm_len = 50;
    char *chrm = (char *)malloc(chrm_len * sizeof(char));
    if (chrm == NULL)
        err(1, "malloc error in giggle_bulk_insert_build_leaf_levels()");

    uint32_t start, end;
    long offset;
    kstring_t line = {0, 0, NULL};

    priority pri_end;
    struct pq_data *pqd_end;
    // Loop over the start queue until it is empty
    while (priq_top(*pq_start, &pri_start) != NULL) {
        // Grab the top element
        pqd_start = (struct pq_data *)priq_pop(*pq_start, &pri_start);

        /* The posibilities for this start position are that:
         * 1) it has been seen before, in which case we will need to add the
         * interval id associated with that position to the starts leaf data
         * and leave the bp tree node alone
         * 2) it has not been seen before, so it will need to be eventually
         * added it to the tree we need to first let the ends catch up by
         * popping any end that is less than to the start that was just seen 3)
         * it is on a new chromosome and we need to do everything that is in 2)
         * as well as close out the disk store for the current chrom and start
         * a new one
         *
         * - Every start and end must be added to the starts and ends arrays.
         * - Any time a new node is created, we need to move the leading,
         *   starts, and ends arrays to a leaf node, and reset the arrays
         */

        if ((pri_start.pos == curr_pos) && 
            (strcmp(curr_chrm, pri_start.chrm) == 0)) {
            // The key didnt' change, so append the current
            // interval id to the end of the leaf data starts

            uint64_t idx = uint64_t_array_add(starts, pqd_start->interval_id);
            // bump starts 
            uint32_t_array_set(starts_pointers,
                               (uint32_t)(starts->num),
                               BPT_NUM_KEYS(bpn) - 1);
            uint32_t_array_set(ends_pointers,
                               (uint32_t)(ends->num),
                               BPT_NUM_KEYS(bpn) - 1);
            if (starts_pointers->num > ORDER)
                errx(1, "too many start pointers");
             
            // Add interval to tree to track intervals for leading value
#if DEBUG
            fprintf(stderr,
                    "-> %s %u %llu\n",
                    pri_start.chrm,
                    pri_start.pos,
                    pqd_start->interval_id);
#endif

            jsw_avlinsert(avl, &(pqd_start->interval_id));
        } else {
            //{{{ #2, we need to go through the ends to catch up
            pqd_end = (struct pq_data *)priq_top(*pq_end, &pri_end);

            // Since the key changed, flush out the ends to this or new keys
            // up to the value of the next start
            while ( (pqd_end != NULL) && //not empy
                    ((strcmp(pri_start.chrm, pri_end.chrm) != 0) || //same chr
                     (pri_end.pos < pri_start.pos)) ) { // < the start we saw

                pqd_end = (struct pq_data *)priq_pop(*pq_end, &pri_end);

                if (curr_pos == pri_end.pos)  {
                    // The key didnt' change, so append the current
                    // interval id to the end of the leaf data ends
                    uint64_t idx = 
                            uint64_t_array_add(ends,
                                               pqd_end->interval_id);
                    // bump ends
                    uint32_t_array_set(starts_pointers,
                                       (uint32_t)(starts->num),
                                       BPT_NUM_KEYS(bpn) - 1);
                    uint32_t_array_set(ends_pointers,
                                       (uint32_t)(ends->num),
                                       BPT_NUM_KEYS(bpn) - 1);


                    if (starts_pointers->num > ORDER)
                        errx(1, "too many start pointers");

                    // remove end from tree tracking leading values
#if DEBUG
                    fprintf(stderr,
                            "<- %s %u %u\n",
                            pri_end.chrm,
                            pri_end.pos,
                            pqd_end->interval_id);
#endif

                    ret = jsw_avlerase(avl, &(pqd_end->interval_id));
                    if (ret == 0)
                        errx(1, "Error removing element from tree.");
                } else {
                    ret = giggle_bulk_insert_append_bpt_key(bpn,
                                                            pri_end.pos,
                                                            curr_ds,
                                                            avl,
                                                            leading,
                                                            starts,
                                                            ends,
                                                            starts_pointers,
                                                            ends_pointers);


                    uint64_t idx =
                            uint64_t_array_add(ends,
                                               pqd_end->interval_id);
                    // bump ends
                    uint32_t_array_set(starts_pointers,
                                       (uint32_t)(starts->num),
                                       BPT_NUM_KEYS(bpn) - 1);
                    uint32_t_array_set(ends_pointers,
                                       (uint32_t)(ends->num),
                                       BPT_NUM_KEYS(bpn) - 1);
                    if (starts_pointers->num > ORDER)
                        errx(1, "too many start pointers");

                    // remove end from tree tracking leading values
#if DEBUG
                    fprintf(stderr,
                            "<- %s %u %u\n",
                            pri_end.chrm,
                            pri_end.pos,
                            pqd_end->interval_id);
#endif
                    ret = jsw_avlerase(avl, &(pqd_end->interval_id));
                    if (ret == 0)
                        errx(1, "Error removing element from tree.");

                    curr_pos = pri_end.pos;
                }

                free(pqd_end);
                pqd_end = (struct pq_data *)priq_top(*pq_end, &pri_end);
            }
            //}}}

            // If the chrom did change, we need to sync up the disk store and
            // open up a new one
            if (strcmp(curr_chrm, pri_start.chrm) != 0) {

                if (BPT_NUM_KEYS(bpn) > 0) {
                    BPT_POINTERS_BLOCK(bpn) = (curr_ds->num + 1) + 1;//1-based
                    BPT_NEXT(bpn) = 0;

                    giggle_bulk_insert_write_leaf_node(bpn,
                                                       curr_ds,
                                                       leading,
                                                       starts,
                                                       ends,
                                                       starts_pointers,
                                                       ends_pointers);
                    // Reset the bpt node
                    memset(bpn->data, 0, BPT_NODE_SIZE);
                    BPT_ID(bpn) =  1;
                    BPT_PARENT(bpn) = 0;
                    BPT_IS_LEAF(bpn) = 1;
                    BPT_LEADING(bpn) = 0;
                    BPT_NEXT(bpn) = 0;
                    BPT_NUM_KEYS(bpn) = 0;
                    BPT_POINTERS_BLOCK(bpn) = 0;
                }

                strcpy(curr_chrm, pri_start.chrm);
                //register the new chrom
                curr_chrm_id = chrm_index_add(gi->chrm_idx,
                                              curr_chrm);
                //{{{ fix up the disk store
                disk_store_sync(curr_ds);
                disk_store_destroy(&curr_ds);

                ret = asprintf(&ds_curr_index_file_name,
                               "%s/%s%u.idx",
                               gi->data_dir,
                               CACHE_FILE_NAME_PREFIX,
                               curr_chrm_id);
                ret = asprintf(&ds_curr_data_file_name,
                               "%s/%s%u.dat",
                               gi->data_dir,
                               CACHE_FILE_NAME_PREFIX,
                               curr_chrm_id);
                curr_ds = disk_store_init(10,
                                          NULL,
                                          ds_curr_index_file_name,
                                          NULL,
                                          ds_curr_data_file_name);
                free(ds_curr_index_file_name);
                free(ds_curr_data_file_name);
                //}}}
            }

            curr_pos = pri_start.pos;
            ret = giggle_bulk_insert_append_bpt_key(bpn,
                                                    curr_pos,
                                                    curr_ds,
                                                    avl,
                                                    leading,
                                                    starts,
                                                    ends,
                                                    starts_pointers,
                                                    ends_pointers);
            uint64_t idx = uint64_t_array_add(starts, pqd_start->interval_id);
            // bump starts
            uint32_t_array_set(starts_pointers,
                               (uint32_t)(starts->num),
                               BPT_NUM_KEYS(bpn) - 1);
            uint32_t_array_set(ends_pointers,
                               (uint32_t)(ends->num),
                               BPT_NUM_KEYS(bpn) - 1);

            if (starts_pointers->num > ORDER)
                errx(1, "too many start pointers");
            // add to tree tracking the leading values
#if DEBUG
            fprintf(stderr,
                    "-> %s %u %u\n",
                    pri_start.chrm,
                    pri_start.pos,
                    pqd_start->interval_id);
#endif

            jsw_avlinsert(avl, &(pqd_start->interval_id));
        }

        //{{{ put another interval from the file that just lost one 
        int ret = i_files[pqd_start->file_id]->
                    input_file_get_next_interval(i_files[pqd_start->file_id],
                                                 &chrm,
                                                 &chrm_len,
                                                 &start,
                                                 &end,
                                                 &offset,
                                                 &line);

        if (ret >= 0) {
            //FIXME
            uint64_t interval_id = offset_index_add(gi->offset_idx,
                                                    offset,
                                                    &line,
                                                    pqd_start->file_id);
            struct file_data *fd = file_index_get(gi->file_idx,
                                                  pqd_start->file_id);
            fd->mean_interval_size += end-start;
            fd->num_intervals += 1;

            pqd_starts[pqd_start->file_id].interval_id = interval_id;
            pri_start.pos = start;
            strcpy(pri_start.chrm, chrm);
            priq_push(*pq_start, &(pqd_starts[pqd_start->file_id]), pri_start);

            pqd_end = (struct pq_data *) malloc(sizeof(struct pq_data));
            if (pqd_end == NULL)
                err(1,
                    "malloc error in giggle_bulk_insert_build_leaf_levels()");

            pqd_end->file_id = pqd_start->file_id;
            pqd_end->interval_id = interval_id;
            pri_end.pos = end + 1;
            strcpy(pri_end.chrm, chrm);
            priq_push(*pq_end, pqd_end, pri_end);
        }
        //}}}
    }

    if (line.s != NULL)
        free(line.s);

    // Once the start queue is empty we need to drain the end queue
    while (priq_top(*pq_end, &pri_end) != NULL) {
        pqd_end = (struct pq_data *)priq_pop(*pq_end, &pri_end);

        if (curr_pos == pri_end.pos)  {
            uint64_t idx = uint64_t_array_add(ends, pqd_end->interval_id);
            // bump ends
            uint32_t_array_set(starts_pointers,
                               (uint32_t)(starts->num),
                               BPT_NUM_KEYS(bpn) - 1);
            uint32_t_array_set(ends_pointers,
                               (uint32_t)(ends->num),
                               BPT_NUM_KEYS(bpn) - 1);


            if (starts_pointers->num > ORDER)
                errx(1, "too many start pointers");

            // remove from tree tracking leading values
#if DEBUG
            fprintf(stderr,
                    "-> %s %u %u\n",
                    pri_end.chrm,
                    pri_end.pos,
                    pqd_end->interval_id);
#endif
            ret = jsw_avlerase(avl, &(pqd_end->interval_id));

            if (ret == 0)
                errx(1, "Error removing element from tree.");
        } else {
            curr_pos = pri_end.pos;
            ret = giggle_bulk_insert_append_bpt_key(bpn,
                                                    curr_pos,
                                                    curr_ds,
                                                    avl,
                                                    leading,
                                                    starts,
                                                    ends,
                                                    starts_pointers,
                                                    ends_pointers);
            uint64_t idx = uint64_t_array_add(ends, pqd_end->interval_id);
            // bump ends
            uint32_t_array_set(starts_pointers,
                               (uint32_t)(starts->num),
                               BPT_NUM_KEYS(bpn) - 1);
            uint32_t_array_set(ends_pointers,
                               (uint32_t)(ends->num),
                               BPT_NUM_KEYS(bpn) - 1);


            if (starts_pointers->num > ORDER)
                errx(1, "too many start pointers");
#if DEBUG
            fprintf(stderr,
                    "-> %s %u %u\n",
                    pri_end.chrm,
                    pri_end.pos,
                    pqd_end->interval_id);
#endif
            // remove from tree tracking leading values
            ret = jsw_avlerase(avl, &(pqd_end->interval_id));
            if (ret == 0)
                errx(1, "Error removing element from tree.");
        }

        free(pqd_end);
    }
    

    // Write out the data if there is a partially filled node left
    if (BPT_NUM_KEYS(bpn) > 0) {
        BPT_POINTERS_BLOCK(bpn) = (curr_ds->num + 1) + 1;//1-based
        BPT_NEXT(bpn) = 0;

        giggle_bulk_insert_write_leaf_node(bpn,
                                           curr_ds,
                                           leading,
                                           starts,
                                           ends,
                                           starts_pointers,
                                           ends_pointers);
    }

    disk_store_sync(curr_ds);
    disk_store_destroy(&curr_ds);

    jsw_avldelete(avl);

    free(bpn->data);
    free(bpn);
    uint64_t_array_destroy(&leading);
    uint64_t_array_destroy(&starts);
    uint64_t_array_destroy(&ends);
    uint32_t_array_destroy(&starts_pointers);
    uint32_t_array_destroy(&ends_pointers);
    free(chrm);
}
//}}}

//{{{void giggle_bulk_insert_prime_pqs(struct giggle_index *gi,
void giggle_bulk_insert_prime_pqs(struct giggle_index *gi,
                                  pri_queue *pq_start,
                                  struct pq_data *pqd_starts,
                                  pri_queue *pq_end,
                                  struct input_file **i_files,
                                  uint32_t num_input_files)
{
    // Priority queue of starts
    //pri_queue pq_start = priq_new(results.gl_pathc);
    priority pri_start;
    priority pri_end;

    // We cannot assume that there will be some set numberof ends per file
    // (contained intervals) so we must malloc on each insert
    struct pq_data *pqd_end;

    // Use these to read intervals from files
    int chrm_len = 10;
    char *chrm = (char *)malloc(chrm_len * sizeof(char));
    if (chrm == NULL)
        err(1, "malloc error in giggle_bulk_insert_prime_pqs()");

    uint32_t start, end;
    long offset;
    kstring_t line = {0, 0, NULL};

    uint64_t interval_id = 0;

    // add one interval from each file to the priority queue
    uint32_t i, ret;
    for (i = 0; i < num_input_files; i++) {
        ret = i_files[i]->input_file_get_next_interval(i_files[i],
                                                       &chrm,
                                                       &chrm_len,
                                                       &start,
                                                       &end,
                                                       &offset,
                                                       &line);

        struct file_data *fd = file_index_get(gi->file_idx, i);
        fd->mean_interval_size += end-start;
        fd->num_intervals += 1;

        // register the interval with the offset index
        interval_id = offset_index_add(gi->offset_idx, offset, &line, i);

        //fprintf(stderr, "%s %u %u %u\n", chrm, start, end, interval_id);

        //Update the pq data for the start, use the array to reduce mallocs
        pqd_starts[i].file_id = i;
        pqd_starts[i].interval_id = interval_id;
        pri_start.pos = start;
        strcpy(pri_start.chrm, chrm);
        priq_push(*pq_start, &(pqd_starts[i]), pri_start);

        //Update the pq data for the end
        pqd_end = (struct pq_data *) malloc(sizeof(struct pq_data));
        if (pqd_end == NULL)
            err(1, "malloc error in giggle_bulk_insert_prime_pqs()");
        
        pqd_end->file_id = i;
        pqd_end->interval_id = interval_id;
        pri_end.pos = end + 1; // use end + 1
        strcpy(pri_end.chrm, chrm);
        priq_push(*pq_end, pqd_end, pri_end);
    }

    if (line.s != NULL)
        free(line.s);

    free(chrm);
}
//}}}

//{{{uint32_t giggle_bulk_insert_open_files(char *input_path_name,
uint32_t giggle_bulk_insert_open_files(char *input_path_name,
                                       char *output_dir_name,
                                       struct input_file ***i_files,
                                       struct file_index **file_idx)
{
    glob_t results;
    int ret = glob(input_path_name, 0, NULL, &results);
    if (ret != 0) 
        errx(1,
             "Problem with %s (%s), stopping early\n",
             input_path_name,
             (ret == GLOB_ABORTED ? "filesystem problem" :
             ret == GLOB_NOMATCH ? "no match of pattern" :
             ret == GLOB_NOSPACE ? "no dynamic memory" :
             "unknown problem"));

    //Array of open pre-sorted input files
    *i_files = (struct input_file **)
            malloc(results.gl_pathc * sizeof(struct input_file *));
    if (i_files == NULL)
        err(1, "malloc error in giggle_bulk_insert_open_files()");

    char *file_index_file_name = NULL;
    ret = asprintf(&file_index_file_name,
                   "%s/%s",
                   output_dir_name,
                   FILE_INDEX_FILE_NAME);

    *file_idx = file_index_init(results.gl_pathc,
                                file_index_file_name);
    free(file_index_file_name);
 
    uint32_t i;
    for (i = 0; i < results.gl_pathc; i++) {
        (*i_files)[i] = input_file_init(results.gl_pathv[i]);
        if ((*i_files)[i] == NULL)
            errx(1, "Could not open %s.\n", results.gl_pathv[i]);
        // register the file with the file index
        uint32_t file_id = file_index_add(*file_idx, results.gl_pathv[i]);
        if (i != file_id)
            errx(1,
                 "Error with file_index synchronization. Saw %u, expected %u.",
                 file_id,
                 i);
    }
    uint32_t num_files = results.gl_pathc;
    globfree(&results);
 
    return num_files;
}
//}}}

//{{{int giggle_bulk_insert_append_bpt_key(struct bpt_node *bpn,
int giggle_bulk_insert_append_bpt_key(struct bpt_node *bpn,
                                      uint32_t key_val,
                                      struct disk_store *ds,
                                      jsw_avltree_t *avl,
                                      struct uint64_t_array *leading,
                                      struct uint64_t_array *starts,
                                      struct uint64_t_array *ends,
                                      struct uint32_t_array *starts_pointers,
                                      struct uint32_t_array *ends_pointers)
{
#if DEBUG
   fprintf(stderr,
           "giggle_bulk_insert_append_bpt_key: append %u %u\n",
           BPT_NUM_KEYS(bpn),
           key_val);
#endif

    int ret = 0;
    if (BPT_NUM_KEYS(bpn) == ORDER) {
        BPT_POINTERS_BLOCK(bpn) = (ds->num + 1) + 1;//1-based
        BPT_NEXT(bpn) = (ds->num + 2) + 1;//1-based

        giggle_bulk_insert_write_leaf_node(bpn,
                                           ds,
                                           leading,
                                           starts,
                                           ends,
                                           starts_pointers,
                                           ends_pointers);

        // Reset the bpt node
        memset(bpn->data, 0, BPT_NODE_SIZE);
        BPT_ID(bpn) =  ds->num + 1;//1-based
        BPT_PARENT(bpn) = 0;
        BPT_IS_LEAF(bpn) = 1;
        BPT_LEADING(bpn) = 0;
        BPT_NEXT(bpn) = 0;
        BPT_NUM_KEYS(bpn) = 0;
        BPT_POINTERS_BLOCK(bpn) = 0;

        // populate the leading values for the next leaf node
        jsw_avltrav_t *avl_t = jsw_avltnew();
        uint64_t *id = (uint64_t *)jsw_avltfirst( avl_t, avl);

        while (id != NULL) {
            uint64_t idx = uint64_t_array_add(leading, *id);
            id = (uint64_t *) jsw_avltnext(avl_t);
        }

        if (leading->num > 0)
            BPT_LEADING(bpn) = 1;

        jsw_avltdelete(avl_t);

        ret = 1;
    }

    BPT_KEYS(bpn)[BPT_NUM_KEYS(bpn)] = key_val;
    BPT_NUM_KEYS(bpn) = BPT_NUM_KEYS(bpn) + 1;

    return ret;
}
//}}}

//{{{void giggle_bulk_insert_write_leaf_node(struct bpt_node *bpn,
void giggle_bulk_insert_write_leaf_node(struct bpt_node *bpn,
                                        struct disk_store *ds,
                                        struct uint64_t_array *leading,
                                        struct uint64_t_array *starts,
                                        struct uint64_t_array *ends,
                                        struct uint32_t_array *starts_pointers,
                                        struct uint32_t_array *ends_pointers)
{
#ifdef LEAF_CHECK
//{{{
    uint64_t *_leading_starts =
            (uint64_t *)
            malloc((starts->num + ends->num)*sizeof(uint64_t));
    memcpy(_leading_starts,
           leading->data,
           leading->num * sizeof(uint64_t));
    memcpy(_leading_starts + leading->num,
           starts->data,
           starts->num * sizeof(uint64_t));

    qsort(_leading_starts,
          starts->num + ends->num,
          sizeof(uint64_t),
          uint64_t_cmp);

    uint32_t j, found = 0;
    for (j = 0; j < ends->num; ++j) {
        uint64_t *ret = bsearch(&(ends->data[j]),
                                _leading_starts,
                                starts->num + ends->num,
                                sizeof(uint64_t),
                                uint64_t_cmp);
        if (ret == NULL)
            errx(1, "Error: end with no matching start/n");
        else
            found += 1;
    }
    free(_leading_starts);
//}}}
#endif

#if DEBUG_GIGGLE_BULK_INSERT_WRITE_LEAF_NODE
    fprintf(stderr,
            "giggle_bulk_insert_write_leaf_node\t"
            "ORDER:%u "
            "leading->num:%llu "
            "starts->num:%llu "
            "ends->num:%llu "
            "starts_pointers->num:%u "
            "ends_pointers->num:%u\n",
            ORDER,
            leading->num,
            starts->num,
            ends->num,
            starts_pointers->num,
            ends_pointers->num);
#endif

    // Write the node out to disk
    void *v;
    uint64_t serialized_size = bpt_node_serialize((void *)bpn, &v);
    uint32_t ds_id = disk_store_append(ds,
                                       v,
                                       serialized_size);         
    free(v);

    // Write out the leaf data
    struct leaf_data *ld = (struct leaf_data *)
            malloc(sizeof(struct leaf_data));
    if (ld == NULL)
        err(1, "malloc error in giggle_bulk_insert_write_leaf_node()");

    ld->num_leading = leading->num;
    ld->num_starts = starts->num;
    ld->num_ends = ends->num;

    ld->starts_pointers = (uint32_t *)calloc(ORDER, LEAF_POINTERS_SIZE);
    if (ld->starts_pointers == NULL)
        err(1, "calloc error in giggle_bulk_insert_write_leaf_node.\n");
    memcpy(ld->starts_pointers,
           starts_pointers->data,
           starts_pointers->num*sizeof(uint32_t));

    ld->ends_pointers = (uint32_t *)calloc(ORDER, sizeof(uint32_t));
    memcpy(ld->ends_pointers,
           ends_pointers->data,
           ends_pointers->num*sizeof(uint32_t));

    ld->data = (uint64_t *) 
            malloc((ld->num_leading + ld->num_starts + ld->num_ends) * 
                    sizeof(uint64_t));
    if (ld->data == NULL)
        err(1, "malloc error in giggle_bulk_insert_write_leaf_node()");


    ld->leading = ld->data;
    ld->starts = ld->data + ld->num_leading;
    ld->ends = ld->data + ld->num_leading + ld->num_starts;
    memcpy(ld->leading, leading->data, ld->num_leading * sizeof(uint64_t));

    qsort(ld->leading, ld->num_leading, sizeof(uint64_t), uint64_t_cmp);

    memcpy(ld->starts, starts->data, ld->num_starts * sizeof(uint64_t));
    memcpy(ld->ends, ends->data, ld->num_ends * sizeof(uint64_t));

    serialized_size = leaf_data_serialize((void *)ld, &v);
    ds_id = disk_store_append(ds,
                              v,
                              serialized_size);         
    free(v);

    free(ld->starts_pointers);
    free(ld->ends_pointers);
    free(ld->data);
    free(ld);

    // Reset the leaf data
    leading->num = 0;
    starts->num = 0;
    ends->num = 0;

    memset(starts_pointers->data, 0, starts_pointers->num * sizeof(uint32_t));
    starts_pointers->num = 0;

    memset(ends_pointers->data, 0, ends_pointers->num * sizeof(uint32_t));
    ends_pointers->num = 0;
}
//}}}

//{{{uint32_t giggle_bulk_insert_add_tree_level(struct disk_store *curr_ds,
uint32_t giggle_bulk_insert_add_tree_level(struct disk_store *curr_ds,
                                           uint32_t curr_level_first_id,
                                           uint32_t curr_level_num_nodes,
                                           uint32_t curr_level_is_leaf,
                                           uint32_t *new_level_first_id)

{
    // If the level only has 1 node, then it will become the root
    if (curr_level_num_nodes == 1) {
        *new_level_first_id = curr_level_first_id;
        return 0;
    }

    // We will use this node to hold the data that will be written to disk
    struct bpt_node *new_bpn =
            (struct bpt_node *) malloc(sizeof(struct bpt_node));
    if (new_bpn == NULL)
        err(1, "malloc error in giggle_bulk_insert_add_tree_level()");

    new_bpn->data = 
            (uint32_t *) malloc(BPT_NODE_NUM_ELEMENTS  * sizeof(uint32_t));
    if (new_bpn->data == NULL)
        err(1, "malloc error in giggle_bulk_insert_add_tree_level()");

    memset(new_bpn->data, 0, BPT_NODE_SIZE);
    *new_level_first_id = curr_ds->num + 1;//1-based
    BPT_ID(new_bpn) =  curr_ds->num + 1;//1-based

    uint32_t curr_row_len = 0;

    // put the left most node input the pointers
    BPT_POINTERS(new_bpn)[0] = curr_level_first_id;
    //uint32_t num_leaf_nodes_leaf_data = curr_ds->num;

    uint32_t j, key_i = 0;
    uint64_t size;
    void *v;
    struct bpt_node *bpn_in;
    uint64_t deserialized_size;
    uint64_t curr_node_id = 0;

    for (j = 1; j < curr_level_num_nodes; j+=1) {
        // Read the current node from disk
        if (curr_level_is_leaf == 1)
            curr_node_id = curr_level_first_id + j*2;
        else
            curr_node_id = curr_level_first_id + j;
        v = disk_store_get(curr_ds, curr_node_id - 1, &size);
        deserialized_size = bpt_node_deserialize(v,
                                                 size,
                                                 (void **)&bpn_in); 
#if DEBUG
        fprintf(stderr,
                "%u(%u) ",
                BPT_KEYS(bpn_in)[0],
                BPT_ID(bpn_in));
#endif

        BPT_KEYS(new_bpn)[key_i] = BPT_KEYS(bpn_in)[0];
        BPT_NUM_KEYS(new_bpn) = BPT_NUM_KEYS(new_bpn) + 1;
        BPT_POINTERS(new_bpn)[key_i + 1] = BPT_ID(bpn_in);

        key_i += 1;

        if (key_i == ORDER) {
            // The node is full, write it to disk and reset
            void *v;
            uint64_t serialized_size = bpt_node_serialize((void *)new_bpn, &v);
            uint32_t ds_id = disk_store_append(curr_ds,
                                               v,
                                               serialized_size);         
            free(v);

            memset(new_bpn->data, 0, BPT_NODE_SIZE);
            BPT_ID(new_bpn) =  curr_ds->num + 1;//1-based
            key_i = 0;
            curr_row_len += 1;
        }
    
        free(bpn_in->data);
        free(bpn_in);
        free(v);
        v = NULL;
    }

    if (key_i > 0) {
        void *v;
        uint64_t serialized_size = bpt_node_serialize((void *)new_bpn, &v);
        uint32_t ds_id = disk_store_append(curr_ds,
                                           v,
                                           serialized_size);         
        free(v);
        curr_row_len += 1;
    }

    free(new_bpn->data);
    free(new_bpn);
    return curr_row_len;
}
//}}}

//{{{uint32_t giggle_get_indexed_files(char *index_dir_name,
uint32_t giggle_get_indexed_files(char *index_dir_name,
                                  char ***names,
                                  uint32_t **num_intervals,
                                  double **mean_interval_sizes)
{
    char *file_index_file_name = NULL;
    int ret = asprintf(&file_index_file_name,
                       "%s/%s",
                       index_dir_name,
                       FILE_INDEX_FILE_NAME);
    struct file_index *file_idx = file_index_load(file_index_file_name);
    free(file_index_file_name);

    uint32_t num = file_idx->index->num;
    *names = (char **)malloc(num * sizeof(char *));
    *num_intervals = (uint32_t *)malloc(num * sizeof(uint32_t));
    *mean_interval_sizes = (double *)malloc(num * sizeof(double));

    uint32_t i;
    for (i = 0; i < file_idx->index->num; ++i) {
        struct file_data *fd = file_index_get(file_idx, i);
        (*names)[i] = strdup(fd->file_name);
        (*num_intervals)[i] = fd->num_intervals;
        (*mean_interval_sizes)[i] = fd->mean_interval_size;
    }
    file_index_destroy(&file_idx);

    return num;
}
//}}}

//{{{void block_store_giggle_set_data_handler();
void block_store_giggle_set_data_handler()
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

    giggle_data_handler = uint64_t_ll_giggle_data_handler;

    giggle_data_handler.giggle_collect_intersection =
            giggle_collect_intersection_data_in_block;

    giggle_data_handler.map_intersection_to_offset_list =
            leaf_data_map_intersection_to_offset_list;

}
//

#if 0
//{{{ uint32_t giggle_merge_chrom(char *chrm_string,
uint32_t giggle_merge_chrom(char *chrm_string,
                            struct giggle_index *gi_0,
                            struct indexed_list *file_index_id_map_0,
                            uint32_t gi_0_cache_name_space,
                            struct giggle_index *gi_1,
                            struct indexed_list *file_index_id_map_1,
                            uint32_t gi_1_cache_name_space,
                            struct disk_store *ds,
                            struct file_id_offset_pairs **merged_offset_index)
{
    // Initialize values for tree 0
    CACHE_NAME_SPACE = gi_0_cache_name_space;
    uint32_t chr_id_0 = giggle_get_chrm_id(gi_0, chrm_string);
    uint32_t curr_leaf_id_0;
    int pos_start_id_0;

    // find the left-most leaf node for tree 0
    uint32_t nld_start_id_0 = bpt_find(chr_id_0,
                                       gi_0->root_ids[chr_id_0],
                                       &curr_leaf_id_0, 
                                       &pos_start_id_0,
                                       0);
    struct bpt_node *curr_leaf_0 = cache.get(chr_id_0,
                                             curr_leaf_id_0 - 1,
                                             &bpt_node_cache_handler);
    struct leaf_data *curr_leaf_data_0 = 
            cache.get(chr_id_0,
                      BPT_POINTERS_BLOCK(curr_leaf_0) - 1,
                      &leaf_data_cache_handler);

    // Initialize values for tree 1
    CACHE_NAME_SPACE = gi_1_cache_name_space;
    uint32_t chr_id_1 = giggle_get_chrm_id(gi_1, chrm_string);
    uint32_t curr_leaf_id_1;
    int pos_start_id_1;

    // find the left-most leaf node for tree 1
    uint32_t nld_start_id_1 = bpt_find(chr_id_1,
                                       gi_1->root_ids[chr_id_1],
                                       &curr_leaf_id_1, 
                                       &pos_start_id_1,
                                       0);
    struct bpt_node *curr_leaf_1 = cache.get(chr_id_1,
                                              curr_leaf_id_1 - 1,
                                              &bpt_node_cache_handler);
    struct leaf_data *curr_leaf_data_1 = 
            cache.get(chr_id_1,
                      BPT_POINTERS_BLOCK(curr_leaf_1) - 1,
                      &leaf_data_cache_handler);

    uint32_t i_0 = 0, i_1 = 0;

    // These trees will track the intervals that are currently in context
    jsw_avltree_t *context_tree_0 = jsw_avlnew(int_cmp_f, int_dup_f, int_rel_f);
    jsw_avltree_t *context_tree_1 = jsw_avlnew(int_cmp_f, int_dup_f, int_rel_f);

    /*
     * As the keys/pointers are scanned in each leaf node a new leaf node is
     * built with the merged data. At the same time a new leaf data is built
     * and any spanning nodes must be tracked for the leading data
     *
     * Must keep a list of in context interval ids
     *
     * Interval ids from each tree must be mapped to new ids for the new tree
     *
     * At each start we need to add all of those nodes to context
     *
     * At each end we must remove those from context
     *
     * Anything that is in context when a leaf node is full must be placed into
     * the leading node of the next leaf node
     *
     *
     *
     */

    // offset ids must be mapped from the values in the distisct trees
    // to merged values, offset_id_map_0 tracks the values from tree 0 and 
    // offset_id_map_1 from tree 1
    // the key will be the original id and value will be the merged id
    struct indexed_list *offset_id_map_0 = indexed_list_init(1000,
                                                             sizeof(uint64_t));
    struct indexed_list *offset_id_map_1 = indexed_list_init(1000,
                                                             sizeof(uint64_t));

    // These lists will become the leaf data for the merged node
    uint32_t merged_starts_size = 1000, merged_starts_num = 0;
    uint32_t merged_ends_size = 1000, merged_ends_num = 0;

    uint32_t *merged_starts =
            (uint32_t *)malloc(merged_starts_size * sizeof(uint32_t));
    if (merged_starts  == NULL)
        err(1, "calloc error in giggle_merge_chrom()");

    uint32_t *merged_ends = 
            (uint32_t *)malloc(merged_ends_size * sizeof(uint32_t));
    if (merged_ends  == NULL)
        err(1, "calloc error in giggle_merge_chrom()");

    // Collect the values into this node, then write it and clear 
    struct bpt_node *to_write_node = (struct bpt_node *)
            malloc(sizeof(struct bpt_node));
    if (to_write_node  == NULL)
        err(1, "calloc error in giggle_merge_chrom()");

    to_write_node->data = (uint32_t *)
            malloc(BPT_NODE_NUM_ELEMENTS  * sizeof(uint32_t));
    if (to_write_node->data  == NULL)
        err(1, "calloc error in giggle_merge_chrom()");

    memset(to_write_node->data, 0, BPT_NODE_SIZE);

    BPT_ID(to_write_node) =  ds->num;
    BPT_PARENT(to_write_node) = 0;
    BPT_IS_LEAF(to_write_node) = 1;
    BPT_LEADING(to_write_node) = 0;
    BPT_NEXT(to_write_node) = 0;
    BPT_NUM_KEYS(to_write_node) = 0;
    BPT_POINTERS_BLOCK(to_write_node) = 0;


    uint32_t node_key_i = 0;

    // loop over all the elments in the chrom tree and merge into a single tree
    while (true) {

        if ( (curr_leaf_id_1 == 0 ) && (curr_leaf_id_0 == 0) ) 
            break;

        uint32_t bpt_key_value = 0, bpt_pointer_value = 0;

        // Choose the next lowest value to merge in hhhhhhhh 
        if ((curr_leaf_id_1 == 0 ) ||
            (BPT_KEYS(curr_leaf_0)[i_0]) < (BPT_KEYS(curr_leaf_1)[i_1])) {
            //{{{pick the value from 0 if it is the only one left, or it is
            //less

            bpt_key_value = BPT_KEYS(curr_leaf_0)[i_0];
            
            uint32_t orig_merged_starts_num = merged_starts_num;
            uint32_t orig_merged_ends_num = merged_ends_num;

            //fprintf(stderr,"0: ");
            giggle_merge_leaf_key(curr_leaf_0,
                                  curr_leaf_data_0,
                                  i_0,
                                  context_tree_0,
                                  offset_id_map_0,
                                  file_index_id_map_0,
                                  gi_0->offset_idx->index,
                                  merged_offset_index,
                                  &merged_starts, 
                                  &merged_starts_size, 
                                  &merged_starts_num, 
                                  &merged_ends,
                                  &merged_ends_size,
                                  &merged_ends_num);

            bpt_pointer_value = (merged_starts_num << 16) + merged_ends_num;

            /*
            fprintf(stderr,
                    "merged s:%u-%u,%u e:%u-%u,%u\n",
                    orig_merged_starts_num,
                    merged_starts_num,
                    (orig_merged_starts_num << 16) + merged_starts_num,
                    orig_merged_ends_num,
                    merged_ends_num,
                    (orig_merged_ends_num << 16) + merged_ends_num);
            */
            i_0+=1;
            //}}}
        } else if ((curr_leaf_id_0 == 0 ) ||
                   (BPT_KEYS(curr_leaf_0)[i_0] > BPT_KEYS(curr_leaf_1)[i_1])) {
            //{{{ pick the value from 1 if it is the only one left, or it is
            //less
            bpt_key_value = BPT_KEYS(curr_leaf_1)[i_1];

            uint32_t orig_merged_starts_num = merged_starts_num;
            uint32_t orig_merged_ends_num = merged_ends_num;

            //fprintf(stderr,"1: ");
            giggle_merge_leaf_key(curr_leaf_1,
                                  curr_leaf_data_1,
                                  i_1,
                                  context_tree_1,
                                  offset_id_map_1,
                                  file_index_id_map_1,
                                  gi_1->offset_idx->index,
                                  merged_offset_index,
                                  &merged_starts, 
                                  &merged_starts_size, 
                                  &merged_starts_num, 
                                  &merged_ends,
                                  &merged_ends_size,
                                  &merged_ends_num);

            bpt_pointer_value = (merged_starts_num << 16) + merged_ends_num;

            /*
            fprintf(stderr,
                    "merged s:%u-%u,%u e:%u-%u,%u\n",
                    orig_merged_starts_num,
                    merged_starts_num,
                    (orig_merged_starts_num << 16) + merged_starts_num,
                    orig_merged_ends_num,
                    merged_ends_num,
                    (orig_merged_ends_num << 16) + merged_ends_num);
            */

            i_1+=1;
            //}}}
        } else if ((BPT_KEYS(curr_leaf_0)[i_0]) == 
                   (BPT_KEYS(curr_leaf_1)[i_1])) {
            // {{{ they are equal

            uint32_t orig_merged_starts_num = merged_starts_num;
            uint32_t orig_merged_ends_num = merged_ends_num;

            bpt_key_value = BPT_KEYS(curr_leaf_0)[i_0];

            //fprintf(stderr,"0 1: ");
            giggle_merge_leaf_key(curr_leaf_0,
                                  curr_leaf_data_0,
                                  i_0,
                                  context_tree_0,
                                  offset_id_map_0,
                                  file_index_id_map_0,
                                  gi_0->offset_idx->index,
                                  merged_offset_index,
                                  &merged_starts, 
                                  &merged_starts_size, 
                                  &merged_starts_num, 
                                  &merged_ends,
                                  &merged_ends_size,
                                  &merged_ends_num);


            giggle_merge_leaf_key(curr_leaf_1,
                                  curr_leaf_data_1,
                                  i_1,
                                  context_tree_1,
                                  offset_id_map_1,
                                  file_index_id_map_1,
                                  gi_1->offset_idx->index,
                                  merged_offset_index,
                                  &merged_starts, 
                                  &merged_starts_size, 
                                  &merged_starts_num, 
                                  &merged_ends,
                                  &merged_ends_size,
                                  &merged_ends_num);

            bpt_pointer_value = (merged_starts_num << 16) + merged_ends_num;

            /*
            fprintf(stderr,
                    "merged s:%u-%u,%u e:%u-%u,%u\n",
                    orig_merged_starts_num,
                    merged_starts_num,
                    (orig_merged_starts_num << 16) + merged_starts_num,
                    orig_merged_ends_num,
                    merged_ends_num,
                    (orig_merged_ends_num << 16) + merged_ends_num);
            */

            i_0+=1;
            i_1+=1;
            //}}}
        } else {
            errx(1, "Not > < or ==");
        }

        //{{{ see if we are moving to the next leaf node on tree 0
        if ( (curr_leaf_id_0 > 0) &&
             (i_0 == BPT_NUM_KEYS(curr_leaf_0)) ) {

            curr_leaf_id_0 = BPT_NEXT(curr_leaf_0);

            if (curr_leaf_id_0 != 0) {
                CACHE_NAME_SPACE = gi_0_cache_name_space;
                curr_leaf_0 = cache.get(chr_id_0,
                                        curr_leaf_id_0 - 1,
                                        &bpt_node_cache_handler);
                i_0 = 0;

                curr_leaf_data_0 = 
                    cache.get(chr_id_0,
                              BPT_POINTERS_BLOCK(curr_leaf_0) - 1,
                              &leaf_data_cache_handler);
                /*
                fprintf(stderr,
                        "node_0: #:%u\n",
                        BPT_NUM_KEYS(curr_leaf_0));
                fprintf(stderr,
                        "leaf_0: l:%u s:%u e:%u\n",
                        curr_leaf_data_0->num_leading,
                        curr_leaf_data_0->num_starts,
                        curr_leaf_data_0->num_ends);
                */

            }
        }
        //}}}

        //{{{ see if we are moving to the next leaf node on tree 1
        if ( (curr_leaf_id_1 > 0) &&
             (i_1 == BPT_NUM_KEYS(curr_leaf_1)) ) {

            curr_leaf_id_1 = BPT_NEXT(curr_leaf_1);

            if (curr_leaf_id_1 != 0) {
                CACHE_NAME_SPACE = gi_1_cache_name_space;
                curr_leaf_1 = cache.get(chr_id_1,
                                        curr_leaf_id_1 - 1,
                                        &bpt_node_cache_handler);
                i_1 = 0;

                curr_leaf_data_1 = 
                    cache.get(chr_id_1,
                              BPT_POINTERS_BLOCK(curr_leaf_1) - 1,
                              &leaf_data_cache_handler);
            }
        }
        //}}}

        /*
        fprintf(stderr,
                "i:%u key:%u\tpointer:%u\t%u\t%u\n",
                node_key_i,
                bpt_key_value,
                bpt_pointer_value,
                bpt_pointer_value >> 16,
                bpt_pointer_value & 65535);
        */

        //BPT_ID(to_write_node) =  ds->num;
        //BPT_PARENT(to_write_node) = 0;
        //BPT_IS_LEAF(to_write_node) = 1;
        //BPT_LEADING(to_write_node) = 0;
        //BPT_NEXT(to_write_node) = 0;
        //BPT_NUM_KEYS(to_write_node) = 0;
        //BPT_POINTERS_BLOCK(to_write_node) = ds->num + 1;
        
        BPT_KEYS(to_write_node)[node_key_i] = bpt_pointer_value;
        BPT_POINTERS(to_write_node)[node_key_i] = bpt_pointer_value;

        node_key_i += 1;

        // The current node is full
        if (node_key_i == ORDER - 1) {

            uint32_t j;
            for (j = 0; j < node_key_i; ++j) {
                fprintf(stderr,
                        "%u\t%u %u %u %u\n",
                        j,
                        BPT_ID(to_write_node),
                        BPT_KEYS(to_write_node)[j],
                        (BPT_POINTERS(to_write_node)[j]) >> 16,
                        (BPT_POINTERS(to_write_node)[j])  & 65535);
            }

            // Add a leading node
            if ( ( jsw_avlsize(context_tree_0) > 0 ) || 
                 ( jsw_avlsize(context_tree_0) > 0 ) ) {
                BPT_LEADING(to_write_node) = 1;
                BPT_POINTERS_BLOCK(to_write_node) = ds->num + 1;
            }

            // Everything is stored through the disk_store struct ds
            // write the current node and the leaf node 
            // set up the next node

            fprintf(stderr, "tree_0 size:%zu\n", jsw_avlsize(context_tree_0));
            fprintf(stderr, "tree_1 size:%zu\n", jsw_avlsize(context_tree_1));

            //jsw_avltrav_t *trav = jsw_avltnew ( void );
            //void *jsw_avltfirst ( jsw_avltrav_t *trav, jsw_avltree_t *tree );
            //void *jsw_avltnext ( jsw_avltrav_t *trav )
            //void jsw_avltdelete ( jsw_avltrav_t *trav );

            // if void           jsw_avltdelete ( jsw_avltrav_t *trav );
        //if ( (curr_leaf_id_1 == 0 ) && (curr_leaf_id_0 == 0) ) 
        //Not going to be a next node

            node_key_i = 0;
        }
    }

    if (node_key_i > 0) {
        uint32_t j;
        for (j = 0; j < node_key_i; ++j) {
            fprintf(stderr,
                    "%u\t%u %u %u\n",
                    j,
                    BPT_KEYS(to_write_node)[j],
                    (BPT_POINTERS(to_write_node)[j]) >> 16,
                    (BPT_POINTERS(to_write_node)[j])  & 65535);
        }
    }



    jsw_avldelete(context_tree_0);
    jsw_avldelete(context_tree_1);

    indexed_list_destroy(&offset_id_map_0);
    indexed_list_destroy(&offset_id_map_1);
    return 0;
}
//}}}
//{{{ uint32_t giggle_merge_add_file_index(struct giggle_index *gi,
uint32_t giggle_merge_add_file_index(struct giggle_index *gi,
                                     struct indexed_list *file_index_id_map,
                                     struct unordered_list *merged_file_index)
{
    uint32_t i;
    for (i = 0 ; i < gi->file_idx->index->num; ++i) {
        struct file_data *fd = file_index_get(gi->file_idx, i);

        struct file_data *merged_fd = (struct file_data *)
                malloc(sizeof(struct file_data));
        if (merged_fd  == NULL)
            err(1, "calloc error in giggle_merge_add_file_index()");

        merged_fd->num_intervals = fd->num_intervals;
        merged_fd->mean_interval_size = fd->mean_interval_size;
        merged_fd->file_name = strdup(fd->file_name);

        uint32_t merged_id = merged_file_index->num;

        uint32_t r = unordered_list_add(merged_file_index, merged_fd);

        r = indexed_list_add(file_index_id_map, i, &merged_id);
    }

    return gi->file_idx->index->num;
}
//}}}
//{{{ uint32_t giggle_merge_get_chrm_index(struct giggle_index *gi_0,
uint32_t giggle_merge_chrm_union(struct giggle_index *gi_0,
                                 struct giggle_index *gi_1,
                                 char ***merged_chrm_set)
{
    uint32_t gi_0_num = gi_0->chrm_idx->index->num;
    uint32_t gi_1_num = gi_1->chrm_idx->index->num;
    // Find the union of the two chrom sets
    char **full_chrm_set = 
            (char **)malloc( (gi_0_num + gi_1_num) * sizeof (char *));
    if (*full_chrm_set  == NULL)
        err(1, "calloc error in giggle_merge_chrm_union()");

    uint32_t i;
    for (i = 0; i < gi_0_num; ++i) {
        struct str_uint_pair *p = 
                (struct str_uint_pair *)gi_0->chrm_idx->index->data[i];
        full_chrm_set[i] = strdup(p->str);
    }
    for (i = 0; i < gi_1_num; ++i) {
        struct str_uint_pair *p = 
                (struct str_uint_pair *)gi_1->chrm_idx->index->data[i];
        full_chrm_set[i + gi_0_num] = strdup(p->str);
    }

    qsort(full_chrm_set, gi_0_num + gi_1_num, sizeof(char *), char_p_cmp);

    uint32_t j, num_uniq = 0;

    for (i = 0; i < gi_0_num + gi_1_num; ) {
        num_uniq += 1;
        j = i + 1;
        while ((j < gi_0_num + gi_1_num) &&
               (strcmp(full_chrm_set[i], full_chrm_set[j]) == 0)) {
            j += 1;
        }
        i = j;
    }

    *merged_chrm_set = (char **)malloc(num_uniq * sizeof (char *));
    if (merged_chrm_set  == NULL)
        err(1, "malloc error in giggle_merge_chrm_union()");

    uint32_t merged_chrm_set_i = 0;
    for (i = 0; i < gi_0_num + gi_1_num; ) {
        (*merged_chrm_set)[merged_chrm_set_i] = strdup(full_chrm_set[i]);
        merged_chrm_set_i += 1;
        j = i + 1;
        while ((j < gi_0_num + gi_1_num) &&
               (strcmp(full_chrm_set[i], full_chrm_set[j]) == 0)) {
            j += 1;
        }
        i = j;
    }

    for (i = 0; i < gi_0_num + gi_1_num; ++i)
        free(full_chrm_set[i]);
    free(full_chrm_set);

    return num_uniq;
}
//}}}
//{{{void giggle_merge_leaf_key(struct bpt_node *node,
void giggle_merge_leaf_key(struct bpt_node *node,
                           struct leaf_data *data,
                           uint32_t key_i,
                           jsw_avltree_t *context_tree,
                           struct indexed_list *offset_id_map,
                           struct indexed_list *file_index_id_map,
                           struct file_id_offset_pairs *offset_index,
                           struct file_id_offset_pairs **merged_offset_index,
                           uint32_t **merged_starts, 
                           uint32_t *merged_starts_size, 
                           uint32_t *merged_starts_num, 
                           uint32_t **merged_ends,
                           uint32_t *merged_ends_size,
                           uint32_t *merged_ends_num)
{
    // put all of the offset ids in starts into the tree
    // get merged offset ids for each start
    // add the merged ids to the start for this position
    // take all of the offset ids out of the tree
    fprintf(stderr,
            "%u\t%u\t%u\ts:%u-%u\te:%u-%u"
            "\t",
            key_i,
            BPT_POINTERS(node)[key_i],
            BPT_KEYS(node)[key_i],
            LEAF_DATA_STARTS_START(node,key_i),
            LEAF_DATA_STARTS_END(node,key_i),
            LEAF_DATA_ENDS_START(node,key_i),
            LEAF_DATA_ENDS_END(node,key_i)
            );

    uint32_t *starts = NULL, *ends = NULL;
    uint32_t starts_size = 0, ends_size = 0;
    // get the list of starts and ends at this key
    uint32_t buff_size = leaf_data_get_starts_ends(node,
                                                   data,
                                                   key_i,
                                                   key_i,
                                                   &starts,
                                                   &starts_size,
                                                   &ends,
                                                   &ends_size);

    //grow the merged offset if we need to
    while ((*merged_offset_index)->size < 
           (*merged_offset_index)->num + starts_size) {

        (*merged_offset_index)->size = (*merged_offset_index)->size * 2;

        fprintf(stderr,
                "merged_offset_index size: %" PRIu64 "\n",
                (*merged_offset_index)->size);


        (*merged_offset_index)->vals = (struct file_id_offset_pair *)
            realloc((*merged_offset_index)->vals,
                    (*merged_offset_index)->size * 
                    sizeof(struct file_id_offset_pair));
        if ((*merged_offset_index)->vals  == NULL)
            err(1, "realloc error in giggle_merge_leaf_key()");

        memset((*merged_offset_index)->vals + (*merged_offset_index)->num,
               0,
               ((*merged_offset_index)->size - (*merged_offset_index)->num) *
               sizeof(struct file_id_offset_pair));
    }

    //grow the merged starts if we need to
    while (*merged_starts_size < *merged_starts_num + starts_size) {
        *merged_starts_size =  *merged_starts_size * 2;

        fprintf(stderr, "merged_starts: %u\n", *merged_starts_size);

        *merged_starts = (uint32_t *)
            realloc(*merged_starts, *merged_starts_size * sizeof(uint32_t)); 
        if (merged_starts == NULL)
            err(1, "realloc error in giggle_merge_leaf_key()");
    }
  
    //grow the merged ends if we need to
    while (*merged_ends_size < *merged_ends_num + ends_size) {
        *merged_ends_size =  *merged_ends_size * 2;

        fprintf(stderr, "merged_ends: %u\n", *merged_ends_size);

        *merged_ends = (uint32_t *)
            realloc(*merged_ends, *merged_ends_size * sizeof(uint32_t)); 
        if (merged_ends == NULL)
            err(1, "realloc error in giggle_merge_leaf_key()");
    }

    // loop over the starts
    // get the merged id
    // add to context tree
    // - add to new offset_index
    // - add to new starts list at this position
    uint32_t j;
    fprintf(stderr, "starts:\t");
    for (j = 0; j < starts_size; ++j) {
        int r = indexed_list_add(offset_id_map,
                                 starts[j],
                                 &((*merged_offset_index)->num));

        fprintf(stderr,
                "oid:%u,mid:%" PRIu64 " ",
                starts[j],
                (*merged_offset_index)->num);

        r = jsw_avlinsert(context_tree, starts + j);


        uint32_t *merged_file_id = 
                indexed_list_get(file_index_id_map, 
                                 offset_index->vals[starts[j]].file_id);
        //(*merged_offset_index)->vals[(*merged_offset_index)->num].file_id =
            //offset_index->vals[starts[j]].file_id;
            
        (*merged_offset_index)->vals[(*merged_offset_index)->num].file_id =
            *merged_file_id;

        (*merged_offset_index)->vals[(*merged_offset_index)->num].offset =
            offset_index->vals[starts[j]].offset;

        (*merged_starts)[*merged_starts_num] = (*merged_offset_index)->num;
        *merged_starts_num = *merged_starts_num + 1;

        (*merged_offset_index)->num = (*merged_offset_index)->num + 1;
    }
    fprintf(stderr, "\t");

    // loop over ends
    // remove from context tree
    // get merged id
    // - add to new ends list at this position
    fprintf(stderr, "ends:\t");
    for (j = 0; j < ends_size; ++j) {
        uint32_t *merged_id = (uint32_t *)
                indexed_list_get(offset_id_map, ends[j]);
        fprintf(stderr,
                "oid:%u,mid:%u ",
                ends[j], 
                *merged_id);

        (*merged_ends)[*merged_ends_num] = *merged_id;
        *merged_ends_num = *merged_ends_num + 1;

        int r = jsw_avlerase(context_tree, ends + j);
    }
    fprintf(stderr, "\n");

    if (starts != NULL)
        free(starts);
    if (ends != NULL)
        free(ends);
}
//}}}
#endif

