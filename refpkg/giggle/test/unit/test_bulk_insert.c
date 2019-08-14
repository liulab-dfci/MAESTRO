#define _GNU_SOURCE

#include "util.h"

#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <inttypes.h>
#include <string.h>
#include <err.h>
#include <sysexits.h>
#include <glob.h>
#include <string.h>
#include <htslib/kstring.h>

#include "unity.h"
#include "bpt.h"
#include "giggle_index.h"
#include "lists.h"
#include "file_read.h"
#include "wah.h"
#include "ll.h"
#include "jsw_avltree.h"
#include "pq.h"

#define MIN(X, Y) (((X) < (Y)) ? (X) : (Y))
#define MAX(X, Y) (((X) > (Y)) ? (X) : (Y))


void setUp(void) { }
void tearDown(void) { }
//
////{{{ void test_giggle_bulk_insert(void)
//void test_giggle_bulk_insert_base(void)
//{
//    //{{{ Testing constants
//    uint32_t BPT_keys_0[100] = {1, 50481, 50514, 144527, 144577, 404661, 404708,
//        422395, 422436, 562161, 562202, 566581, 567102, 585763, 585804, 816018,
//        816215, 829181, 829221, 935394, 936147, 936482, 936605, 1004540,
//        1004903, 1014910, 1015319, 1243704, 1244171, 1284817, 1284894, 1333543,
//        1334139, 1334987, 1335286, 1406946, 1407324, 1446566, 1446696, 1550487,
//        1551148, 1585271, 1585317, 1624220, 1624282, 1689977, 1690749, 1705115,
//        1705170, 1713720, 1714070, 1777263, 1777312, 1850348, 1850535, 2023080,
//        2023417, 2050152, 2050195, 2111658, 2112294, 2125435, 2125603, 2139348,
//        2139395, 2142838, 2142876, 2144034, 2144521, 2322920, 2323328, 2517907,
//        2518312, 2540473, 2540521, 2574446, 2574807, 2741456, 2741507, 3142094,
//        3142139, 3190381, 3190416, 3252166, 3252207, 3455260, 3455310, 3578853,
//        3579009, 3585079, 3585112, 3612184, 3612525, 3636454, 3636501, 3689319,
//        3691116, 3696692, 3696741, 3774727}; 
//
//    uint32_t BPT_keys_1[100] = {3774986, 3809712, 3809783, 3825025, 3825065,
//        4046152, 4046215, 4324287, 4324327, 4364646, 4364680, 4413593, 4413626,
//        4480292, 4480341, 4569164, 4569196, 4582988, 4583026, 4643102, 4643134,
//        4643898, 4643933, 4739868, 4739913, 4873368, 4873409, 4992134, 4992166,
//        5068291, 5068337, 5068817, 5068850, 5082558, 5082596, 5309091, 5309133,
//        5336759, 5336798, 5748446, 5748479, 5781921, 5782215, 5866231, 5866265,
//        5945881, 5946136, 5976058, 5976417, 6077466, 6077500, 6087279, 6087335,
//        6187730, 6187938, 6192940, 6192971, 6219287, 6219320, 6289843, 6289885,
//        6383797, 6383831, 6423193, 6423390, 6430517, 6430567, 6453604, 6454441,
//        6474367, 6474735, 6483544, 6483800, 6486473, 6486687, 6769020, 6769062,
//        6793250, 6793281, 6800601, 6800646, 6849663, 6849671, 6850435, 6850470,
//        6944351, 6944397, 6959614, 6959655, 6962099, 6962134, 7134888, 7134899,
//        7172076, 7172248, 7189979, 7190016, 7257688, 7257739, 7278281};
//
//    uint32_t bpt_node_tested_0 = 0, bpt_node_tested_1 = 0;
//    //}}}
//
//    char *path_name = "../data/many/*gz";
//
//    glob_t results;
//    int ret = glob(path_name, 0, NULL, &results);
//    if (ret != 0) 
//        errx(1,
//             "Problem with %s (%s), stopping early\n",
//             path_name,
//             (ret == GLOB_ABORTED ? "filesystem problem" :
//             ret == GLOB_NOMATCH ? "no match of pattern" :
//             ret == GLOB_NOSPACE ? "no dynamic memory" :
//             "unknown problem"));
//
//    //Array of open pre-sorted input files
//    struct input_file **i_files = (struct input_file **)
//            malloc(results.gl_pathc * sizeof(struct input_file *));
//
//    // Use these to read intervals from files
//    int chrm_len = 10;
//    char *chrm = (char *)malloc(chrm_len * sizeof(char));
//    uint32_t start, end;
//    long offset;
//    kstring_t line = {0, 0, NULL};
//
//    // Priority queue of starts
//    pri_queue pq_start = priq_new(results.gl_pathc);
//    priority pri_start;
//    // Since we know that each file will have at most one start in the priority
//    // queue, we can reduce mallocs by reusing the array
//    struct pq_data *pqd_starts = (struct pq_data *)
//            malloc(results.gl_pathc * sizeof(struct pq_data));
//
//    // Priority queue of ends
//    pri_queue pq_end = priq_new(results.gl_pathc);
//    priority pri_end;
//    // We cannot assume that there will be some set numberof ends per file
//    // (contained intervals) so we must malloc on each insert
//    struct pq_data *pqd_end;
//
//    // Init other associated indexes
//    struct file_index *file_idx = file_index_init(10,
//                                                  "test_bulk_insert.file.idx");
//    uint32_t interval_id = 0;
//    struct offset_index *offset_idx = 
//            offset_index_init(1000,
//                              "test_bulk_insert.offset.idx");
// 
//    struct chrm_index *chrm_idx = chrm_index_init(24,
//                                                  "test_bulk_insert.chrm.idx");
//   
//    //{{{ add one interval from each file to the priority queue
//
//    uint32_t i;
//    uint32_t num_input_files = results.gl_pathc;
//    for (i = 0; i < results.gl_pathc; i++) {
//        i_files[i] = input_file_init(results.gl_pathv[i]);
//        // register the file with the file index
//        uint32_t file_id = file_index_add(file_idx, results.gl_pathv[i]);
//        TEST_ASSERT_EQUAL(i, file_id);
//        ret = i_files[i]->input_file_get_next_interval(i_files[i],
//                                                       &chrm,
//                                                       &chrm_len,
//                                                       &start,
//                                                       &end,
//                                                       &offset,
//                                                       &line);
//        // register the interval with the offset index
//        interval_id = offset_index_add(offset_idx,
//                                       offset,
//                                       &line,
//                                       file_id);
//
//        //Update the pq data for the start, use the array to reduce mallocs
//        pqd_starts[i].file_id = file_id;
//        pqd_starts[i].interval_id = interval_id;
//        pri_start.pos = start;
//        strcpy(pri_start.chrm, chrm);
//        priq_push(pq_start, &(pqd_starts[i]), pri_start);
//
//        //{{{ debug
//        /*
//        fprintf(stderr,
//                "%s\t%u %u %u\n",
//                results.gl_pathv[i],
//                interval_id,
//                pri.pos, start);
//        */
//        //}}}
//
//        //Update the pq data for the end
//        pqd_end = (struct pq_data *) malloc(sizeof(struct pq_data));
//        pqd_end->file_id = file_id;
//        pqd_end->interval_id = interval_id;
//        pri_end.pos = end + 1; // use end + 1
//        strcpy(pri_end.chrm, chrm);
//        priq_push(pq_end, pqd_end, pri_end);
//    }
//    globfree(&results);
//
//    //}}}
//
//    TEST_ASSERT_EQUAL(22, file_idx->index->num);
//
//    struct pq_data *pqd_start =
//            (struct pq_data *)priq_top(pq_start, &pri_start);
//
//    // curr_pos and curr_chrm track the status of the indexing
//    uint32_t curr_pos = pri_start.pos;
//    char curr_chrm[10];
//    strcpy(curr_chrm, pri_start.chrm);
//
//    // register the chrom with chrom index
//    uint32_t curr_chrm_id = chrm_index_add(chrm_idx, curr_chrm);
//
//    //{{{ init disk store, do this at the start of every chrom
//    char *ds_curr_index_file_name = NULL, *ds_curr_data_file_name = NULL;
//    ret = asprintf(&ds_curr_index_file_name,
//                   "test_bulk_insert.ds_idx.%u",
//                   curr_chrm_id);
//    ret = asprintf(&ds_curr_data_file_name,
//                   "test_bulk_insert.ds_data.%u",
//                   curr_chrm_id);
//    struct disk_store *curr_ds = disk_store_init(10,
//                                                 NULL,
//                                                 ds_curr_index_file_name,
//                                                 NULL,
//                                                 ds_curr_data_file_name);
//    free(ds_curr_index_file_name);
//    free(ds_curr_data_file_name);
//    //}}}
//
//    // Collect the values into this node, then write it and clear 
//    struct bpt_node *bpn = (struct bpt_node *) malloc(sizeof(struct bpt_node));
//    bpn->data = (uint32_t *) malloc(BPT_NODE_NUM_ELEMENTS  * sizeof(uint32_t));
//    memset(bpn->data, 0, BPT_NODE_SIZE);
//
//    //BPT_ID(bpn) =  curr_ds->num;
//    BPT_ID(bpn) = 1;
//    BPT_PARENT(bpn) = 0;
//    BPT_IS_LEAF(bpn) = 1;
//    BPT_LEADING(bpn) = 0;
//    BPT_NEXT(bpn) = 0;
//    BPT_NUM_KEYS(bpn) = 0;
//    BPT_POINTERS_BLOCK(bpn) = 0;
//
//    // These will be used to create the leaf data for each node
//    uint32_t num_leading = 0, num_starts = 0, num_ends = 0;
//    struct uint64_t_array *leading, *starts, *ends;
//    leading = uint64_t_array_init(100);
//    starts = uint64_t_array_init(100);
//    ends = uint64_t_array_init(100);
//
//    struct uint32_t_array *starts_pointers, *ends_pointers;
//    starts_pointers = uint32_t_array_init(100);
//    ends_pointers = uint32_t_array_init(100);
//
//
//    // This tree will track intervals that have begun and not yet ended and
//    // will be used to populate the leading value of nodes
//    jsw_avltree_t *avl = jsw_avlnew(uint_cmp_f, uint_dup_f, uint_rel_f);
//
//    // add the current possition to the node
//    ret = giggle_bulk_insert_append_bpt_key(bpn,
//                                            curr_pos,
//                                            curr_ds,
//                                            avl,
//                                            leading,
//                                            starts,
//                                            ends,
//                                            starts_pointers,
//                                            ends_pointers);
//
//    // Loop over the start queue until it is empty
//    while (priq_top(pq_start, &pri_start) != NULL) {
//        // Grab the top element
//        pqd_start = (struct pq_data *)priq_pop(pq_start, &pri_start);
//        //fprintf(stderr, "%s s:%u\n", pri_start.chrm, pri_start.pos);
//
//        /* The posibilities for this start position are that:
//         * 1) it has been seen before, in which case we will need to add the
//         * interval id associated with that position to the starts leaf data
//         * and leave the bp tree node alone
//         * 2) it has not been seen before, so it will need to be eventually
//         * added it to the tree we need to first let the ends catch up by
//         * popping any end that is less than to the start that was just seen 3)
//         * it is on a new chromosome and we need to do everything that is in 2)
//         * as well as close out the disk store for the current chrom and start
//         * a new one
//         *
//         * - Every start and end must be added to the starts and ends arrays.
//         * - Any time a new node is created, we need to move the leading,
//         *   starts, and ends arrays to a leaf node, and reset the arrays
//         */
//
//        if ((pri_start.pos == curr_pos) && 
//            (strcmp(curr_chrm, pri_start.chrm) == 0)) {
//            // The key didnt' change, so append the current
//            // interval id to the end of the leaf data starts
//
//            uint64_t idx = uint64_t_array_add(starts, pqd_start->interval_id);
//            // bump starts 
//            giggle_bulk_insert_set_starts(bpn, idx);
//
//            // Add interval to tree to track intervals for leading value
//#if DEBUG
//            fprintf(stderr,
//                    "-> %s %u %u\n",
//                    pri_start.chrm,
//                    pri_start.pos,
//                    pqd_start->interval_id);
//#endif
//
//            jsw_avlinsert(avl, &(pqd_start->interval_id));
//        } else {
//            //{{{ #2, we need to go through the ends to catch up
//            pqd_end = (struct pq_data *)priq_top(pq_end, &pri_end);
//
//            // Since the key changed, flush out the ends to this or new keys
//            // up to the value of the next start
//            while ( (pqd_end != NULL) && //not empy
//                    ((strcmp(pri_start.chrm, pri_end.chrm) != 0) || //same chr
//                     (pri_end.pos < pri_start.pos)) ) { // < the start we saw
//
//                pqd_end = (struct pq_data *)priq_pop(pq_end, &pri_end);
//                //fprintf(stderr, "%s e:%u\n", pri_end.chrm, pri_end.pos);
//
//                if (curr_pos == pri_end.pos)  {
//                    // The key didnt' change, so append the current
//                    // interval id to the end of the leaf data ends
//                    uint64_t idx = 
//                            uint64_t_array_add(ends,
//                                               pqd_end->interval_id);
//                    // bump ends
//                    giggle_bulk_insert_set_ends(bpn, idx);
//
//                    // remove end from tree tracking leading values
//#if DEBUG
//                    fprintf(stderr,
//                            "<- %s %u %u\n",
//                            pri_end.chrm,
//                            pri_end.pos,
//                            pqd_end->interval_id);
//#endif
//
//                    ret = jsw_avlerase(avl, &(pqd_end->interval_id));
//                    if (ret == 0)
//                        errx(1, "Error removing element from tree.");
//                } else {
//                    ret = giggle_bulk_insert_append_bpt_key(bpn,
//                                                            pri_end.pos,
//                                                            curr_ds,
//                                                            avl,
//                                                            leading,
//                                                            starts,
//                                                            ends,
//                                                            start_poiners,
//                                                            ends_pointers);
//
//                    uint64_t idx =
//                            uint64_t_array_add(ends,
//                                               pqd_end->interval_id);
//                    // bump ends
//                    //giggle_bulk_insert_set_ends(bpn, idx);
//
//                    // remove end from tree tracking leading values
//#if DEBUG
//                    fprintf(stderr,
//                            "<- %s %u %u\n",
//                            pri_end.chrm,
//                            pri_end.pos,
//                            pqd_end->interval_id);
//#endif
//                    ret = jsw_avlerase(avl, &(pqd_end->interval_id));
//                    if (ret == 0)
//                        errx(1, "Error removing element from tree.");
//
//
//                    curr_pos = pri_end.pos;
//
//                    //{{{ bpt node contents test
//                    if (BPT_NUM_KEYS(bpn) == ORDER) {
//                        if (bpt_node_tested_0 == 0) {
//                            uint32_t i;
//                            for (i = 0; i < ORDER; ++i)
//                                TEST_ASSERT_EQUAL(BPT_keys_0[i],
//                                                  BPT_KEYS(bpn)[i]);
//                            bpt_node_tested_0 = 1;
//                        } else if (bpt_node_tested_1 == 0) {
//                            uint32_t i;
//                            for (i = 0; i < ORDER; ++i)
//                                TEST_ASSERT_EQUAL(BPT_keys_1[i],
//                                                  BPT_KEYS(bpn)[i]);
//                            bpt_node_tested_1 = 1;
//                        }
//                    }
//                    //}}}
//                }
//
//                free(pqd_end);
//                pqd_end = (struct pq_data *)priq_top(pq_end, &pri_end);
//            }
//            //}}}
//
//            // If the chrom did change, we need to sync up the disk store and
//            // open up a new one
//            if (strcmp(curr_chrm, pri_start.chrm) != 0) {
//
//                if (BPT_NUM_KEYS(bpn) > 0) {
//                    BPT_POINTERS_BLOCK(bpn) = (curr_ds->num + 1) + 1;//1-based
//                    BPT_NEXT(bpn) = 0;
//
//                    giggle_bulk_insert_write_leaf_node(bpn,
//                                                       curr_ds,
//                                                       leading,
//                                                       starts,
//                                                       ends,
//                                                       starts_pointers,
//                                                       ends_pointers);
//                    // Reset the bpt node
//                    memset(bpn->data, 0, BPT_NODE_SIZE);
//                    BPT_ID(bpn) =  1;
//                    BPT_PARENT(bpn) = 0;
//                    BPT_IS_LEAF(bpn) = 1;
//                    BPT_LEADING(bpn) = 0;
//                    BPT_NEXT(bpn) = 0;
//                    BPT_NUM_KEYS(bpn) = 0;
//                    BPT_POINTERS_BLOCK(bpn) = 0;
//                }
//
//                strcpy(curr_chrm, pri_start.chrm);
//                //register the new chrom
//                curr_chrm_id = chrm_index_add(chrm_idx,
//                                              curr_chrm);
//                //{{{ fix up the disk store
//                disk_store_sync(curr_ds);
//                disk_store_destroy(&curr_ds);
//                ret = asprintf(&ds_curr_index_file_name,
//                               "test_bulk_insert.ds_idx.%u",
//                               curr_chrm_id);
//                ret = asprintf(&ds_curr_data_file_name,
//                               "test_bulk_insert.ds_data.%u",
//                               curr_chrm_id);
//                curr_ds = disk_store_init(10,
//                                          NULL,
//                                          ds_curr_index_file_name,
//                                          NULL,
//                                          ds_curr_data_file_name);
//                free(ds_curr_index_file_name);
//                free(ds_curr_data_file_name);
//                //}}}
//            }
//
//            curr_pos = pri_start.pos;
//            ret = giggle_bulk_insert_append_bpt_key(bpn,
//                                                    curr_pos,
//                                                    curr_ds,
//                                                    avl,
//                                                    leading,
//                                                    starts,
//                                                    ends);
//            uint64_t idx = uint64_t_array_add(starts, pqd_start->interval_id);
//            // bump starts
//            giggle_bulk_insert_set_starts(bpn, idx);
//
//            // add to tree tracking the leading values
//#if DEBUG
//            fprintf(stderr,
//                    "-> %s %u %u\n",
//                    pri_start.chrm,
//                    pri_start.pos,
//                    pqd_start->interval_id);
//#endif
//
//            jsw_avlinsert(avl, &(pqd_start->interval_id));
//            //{{{ bpt node contents test
//            if (BPT_NUM_KEYS(bpn) == ORDER) {
//                if (bpt_node_tested_0 == 0) {
//                    uint32_t i;
//                    for (i = 0; i < ORDER; ++i)
//                        TEST_ASSERT_EQUAL(BPT_keys_0[i],
//                                          BPT_KEYS(bpn)[i]);
//                    bpt_node_tested_0 = 1;
//                } else if (bpt_node_tested_1 == 0) {
//                    uint32_t i;
//                    for (i = 0; i < ORDER; ++i)
//                        TEST_ASSERT_EQUAL(BPT_keys_1[i],
//                                          BPT_KEYS(bpn)[i]);
//                    bpt_node_tested_1 = 1;
//                }
//            }
//            //}}}
//            //{{{ debug
//            /*
//            fprintf(stderr,
//                    "%s %u\ts: %u(%u)",
//                    curr_chrm,
//                    curr_pos,
//                    pqd_start->interval_id,
//                    pri.pos);
//            */
//            //}}}
//        }
//
//        //{{{ put another interval from the file that just lost one 
//        int ret = i_files[pqd_start->file_id]->
//                    input_file_get_next_interval(i_files[pqd_start->file_id],
//                                                 &chrm,
//                                                 &chrm_len,
//                                                 &start,
//                                                 &end,
//                                                 &offset,
//                                                 &line);
//
//        if (ret >= 0) {
//            interval_id = offset_index_add(offset_idx,
//                                           offset,
//                                           &line,
//                                           pqd_start->file_id);
//
//            pqd_starts[pqd_start->file_id].interval_id = interval_id;
//            pri_start.pos = start;
//            strcpy(pri_start.chrm, chrm);
//            priq_push(pq_start, &(pqd_starts[pqd_start->file_id]), pri_start);
//
//            pqd_end = (struct pq_data *) malloc(sizeof(struct pq_data));
//            pqd_end->file_id = pqd_start->file_id;
//            pqd_end->interval_id = interval_id;
//            pri_end.pos = end + 1;
//            strcpy(pri_end.chrm, chrm);
//            priq_push(pq_end, pqd_end, pri_end);
//        }
//        //}}}
//    }
//
//    if (line.s != NULL)
//        free(line.s);
//
//    // Once the start queue is empty we need to drain the end queue
//    while (priq_top(pq_end, &pri_end) != NULL) {
//        pqd_end = (struct pq_data *)priq_pop(pq_end, &pri_end);
//        //fprintf(stderr, "%s e:%u\n", pri_end.chrm, pri_end.pos);
//
//        if (curr_pos == pri_end.pos)  {
//            uint64_t idx = uint64_t_array_add(ends, pqd_end->interval_id);
//            // bump ends
//            giggle_bulk_insert_set_ends(bpn, idx);
//
//            // remove from tree tracking leading values
//#if DEBUG
//            fprintf(stderr,
//                    "-> %s %u %u\n",
//                    pri_end.chrm,
//                    pri_end.pos,
//                    pqd_end->interval_id);
//#endif
//            ret = jsw_avlerase(avl, &(pqd_end->interval_id));
//
//            if (ret == 0)
//                errx(1, "Error removing element from tree.");
//            //{{{debug
//            /*
//            fprintf(stderr,
//                    " %u(%u)",
//                    pqd_end->interval_id,
//                    pri_end.pos);
//            */
//            //}}}
//        } else {
//            curr_pos = pri_end.pos;
//            ret = giggle_bulk_insert_append_bpt_key(bpn,
//                                                    curr_pos,
//                                                    curr_ds,
//                                                    avl,
//                                                    leading,
//                                                    starts,
//                                                    ends);
//            uint64_t idx = uint64_t_array_add(ends, curr_pos); 
//            // bump ends
//            giggle_bulk_insert_set_ends(bpn, idx);
//
//#if DEBUG
//            fprintf(stderr,
//                    "-> %s %u %u\n",
//                    pri_end.chrm,
//                    pri_end.pos,
//                    pqd_end->interval_id);
//#endif
//            // remove from tree tracking leading values
//            ret = jsw_avlerase(avl, &(pqd_end->interval_id));
//            if (ret == 0)
//                errx(1, "Error removing element from tree.");
//
//            //{{{ test
//            if (BPT_NUM_KEYS(bpn) == ORDER) {
//                if (bpt_node_tested_0 == 0) {
//                    uint32_t i;
//                    for (i = 0; i < ORDER; ++i)
//                        TEST_ASSERT_EQUAL(BPT_keys_0[i],
//                                          BPT_KEYS(bpn)[i]);
//                    bpt_node_tested_0 = 1;
//                } else if (bpt_node_tested_1 == 0) {
//                    uint32_t i;
//                    for (i = 0; i < ORDER; ++i)
//                        TEST_ASSERT_EQUAL(BPT_keys_1[i],
//                                          BPT_KEYS(bpn)[i]);
//                    bpt_node_tested_1 = 1;
//                }
//            }
//            //}}}
//            //{{{ debug
//            /*
//            fprintf(stderr,
//                    "\n%s %u\ts:\te: %u(%u)",
//                    curr_chrm,
//                    curr_pos,
//                    pqd_end->interval_id,
//                    pri_end.pos);
//            */
//            //}}}
//        }
//
//        free(pqd_end);
//    }
//    
//    if (BPT_NUM_KEYS(bpn) > 0) {
//        BPT_POINTERS_BLOCK(bpn) = (curr_ds->num + 1) + 1;//1-based
//        BPT_NEXT(bpn) = 0;
//
//        giggle_bulk_insert_write_leaf_node(bpn,
//                                           curr_ds,
//                                           leading,
//                                           starts,
//                                           ends);
//    }
//
//    TEST_ASSERT_EQUAL(0, jsw_avlsize(avl));
//    TEST_ASSERT_EQUAL(24, chrm_idx->index->num);
//    TEST_ASSERT_EQUAL(21024, offset_idx->index->num);
//
//    disk_store_sync(curr_ds);
//    disk_store_destroy(&curr_ds);
//
//    jsw_avldelete(avl);
//
//    for (i = 0; i < num_input_files; ++i)
//        input_file_destroy(&(i_files[i]));
//    free(i_files);
//
//
//    //{{{ tests
//    for (i = 0; i < chrm_idx->index->num; ++i) {
//        ret = asprintf(&ds_curr_index_file_name,
//                       "test_bulk_insert.ds_idx.%u",
//                       i);
//        ret = asprintf(&ds_curr_data_file_name,
//                       "test_bulk_insert.ds_data.%u",
//                       i);
//        curr_ds = disk_store_load(NULL,
//                                  ds_curr_index_file_name,
//                                  NULL,
//                                  ds_curr_data_file_name);
//
//        free(ds_curr_index_file_name);
//        free(ds_curr_data_file_name);
//
//        uint32_t num_starts = 0, num_ends = 0;
//
//        uint32_t j;
//        for (j = 0; j < curr_ds->num; ++j) {
//            uint64_t size;
//            void *v = disk_store_get(curr_ds, j, &size);
//            struct bpt_node *bpn_in;
//            uint64_t deserialized_size =
//                    bpt_node_deserialize(v,
//                                        size,
//                                        (void **)&bpn_in); 
//        
//            TEST_ASSERT_EQUAL(j + 1, BPT_ID(bpn_in));//1-based
//            TEST_ASSERT_EQUAL(1, BPT_IS_LEAF(bpn_in));
//            if (j < curr_ds->num - 2) {
//                TEST_ASSERT_EQUAL(j + 2, BPT_NEXT(bpn_in) - 1);//1-based
//            } else {
//                TEST_ASSERT_EQUAL(0, BPT_NEXT(bpn_in));
//            }
//            TEST_ASSERT_EQUAL(j + 1, BPT_POINTERS_BLOCK(bpn_in) - 1);//1-based
//
//            free(bpn_in->data);
//            free(bpn_in);
//            free(v);
//            v = NULL;
//
//            j += 1;
//        
//            v = disk_store_get(curr_ds, j, &size);
//            struct leaf_data *ld_in;
//            deserialized_size = leaf_data_deserialize(v,
//                                                      size,
//                                                      (void **)&ld_in); 
//            num_starts += ld_in->num_starts;
//            num_ends += ld_in->num_ends;
//
//            /*
//            fprintf(stderr,
//                    "%u %u %u\n",
//                    ld_in->num_leading,
//                    ld_in->num_starts,
//                    ld_in->num_ends);
//            */
//
//            free(ld_in->data);
//            free(ld_in);
//            free(v);
//            v = NULL;
//        }
//
//        TEST_ASSERT_EQUAL(num_starts, num_ends);
//
//        disk_store_destroy(&curr_ds);
//    }
//    //}}}
//
//    for (i = 0; i < chrm_idx->index->num; ++i) {
//       
//        ret = asprintf(&ds_curr_index_file_name,
//                       "test_bulk_insert.ds_idx.%u",
//                       i);
//        ret = asprintf(&ds_curr_data_file_name,
//                       "test_bulk_insert.ds_data.%u",
//                       i);
//        curr_ds = disk_store_load(NULL,
//                                  ds_curr_index_file_name,
//                                  NULL,
//                                  ds_curr_data_file_name);
//        free(ds_curr_index_file_name);
//        free(ds_curr_data_file_name);
//        //fprintf(stderr, "chrm:%u\n", i);
//
//        // Here we need to loop over each level of the tree until the current
//        // level has just one element in which case that element is the root
//
//        uint32_t num_leaf_node_leaf_data = curr_ds->num;
//        uint32_t curr_level_num_nodes = num_leaf_node_leaf_data / 2;
//        uint32_t curr_level_first_id = 1;
//        uint32_t curr_level_is_leaf = 1;
//        uint32_t new_level_first_id = 0;
//        uint32_t new_level_len = 
//                giggle_bulk_insert_add_tree_level(curr_ds,
//                                                  curr_level_first_id,
//                                                  curr_level_num_nodes,
//                                                  curr_level_is_leaf,
//                                                  &new_level_first_id);
//
//        if (curr_level_num_nodes > 1)
//            TEST_ASSERT_EQUAL((curr_level_num_nodes + ORDER - 1) / ORDER, 
//                              new_level_len);
//        else
//            TEST_ASSERT_EQUAL(0, new_level_len);
//
//        disk_store_destroy(&curr_ds);
//    }
//    
//    //{{{ remove files
//    for (i = 0; i < chrm_idx->index->num; ++i) {
//        ret = asprintf(&ds_curr_index_file_name,
//                       "test_bulk_insert.ds_idx.%u",
//                       i);
//        ret = asprintf(&ds_curr_data_file_name,
//                       "test_bulk_insert.ds_data.%u",
//                       i);
//        remove(ds_curr_index_file_name);
//        remove(ds_curr_data_file_name);
//        free(ds_curr_index_file_name);
//        free(ds_curr_data_file_name);
//    }
//    //}}}
//
//    //frees
//    free(chrm);
//    free(pqd_starts);
//    uint64_t_array_destroy(&leading);
//    uint64_t_array_destroy(&starts);
//    uint64_t_array_destroy(&ends);
//    uint32_t_array_destroy(&starts_pointers);
//    uint32_t_array_destroy(&ends_pointers);
//    free(bpn->data);
//    free(bpn);
//    offset_index_destroy(&offset_idx);
//    chrm_index_destroy(&chrm_idx);
//    file_index_destroy(&file_idx);
//    priq_free(pq_start);
//    priq_free(pq_end);
//    remove("test_bulk_insert.offset.idx");
//}
////}}}
//
//{{{ void test_giggle_bulk_insert_open_files(void)
void test_giggle_bulk_insert_open_files(void)
{
    char *input_path_name = "../data/many/*gz";
    char *output_path_name = "tmp_test_giggle_bulk_insert_open_files";

    /*
    struct stat st = {0};
    if (stat(output_path_name, &st) == -1) {
        mkdir(output_path_name, 0700);
    } else {
        rmrf(output_path_name);
        mkdir(output_path_name, 0700);
    }
    */

    struct file_index *file_idx = NULL;
    struct input_file **i_files = NULL;

    uint32_t num_input_files = giggle_bulk_insert_open_files(input_path_name,
                                                             output_path_name,
                                                             &i_files,
                                                             &file_idx);
    TEST_ASSERT_EQUAL(22, num_input_files);

    uint32_t i;
    for (i = 0; i < num_input_files; ++i)
        input_file_destroy(&(i_files[i]));

    file_index_destroy(&file_idx);
    free(i_files);
}
//}}}

//{{{ void test_giggle_bulk_insert_prime_pq(void)
void test_giggle_bulk_insert_prime_pq(void)
{
    char *input_path_name = "../data/many/*gz";
    char *output_path_name = "tmp_test_giggle_bulk_insert_open_files";

    struct stat st = {0};
    if (stat(output_path_name, &st) == -1) {
        mkdir(output_path_name, 0700);
    } else {
        rmrf(output_path_name);
        mkdir(output_path_name, 0700);
    }

    struct giggle_index *gi = (struct giggle_index *)
            malloc(sizeof(struct giggle_index));
    gi->data_dir = strdup(output_path_name);

    // open files
    gi->file_idx = NULL;
    struct input_file **i_files = NULL;
    uint32_t num_input_files = giggle_bulk_insert_open_files(input_path_name,
                                                             gi->data_dir, 
                                                             &i_files,
                                                             &(gi->file_idx));

    TEST_ASSERT_EQUAL(22, num_input_files);

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

    pri_queue pq_end = priq_new(num_input_files);

    giggle_bulk_insert_prime_pqs(gi,
                                 &pq_start,
                                 pqd_starts,
                                 &pq_end,
                                 i_files,
                                 num_input_files);
    
    TEST_ASSERT_EQUAL(num_input_files, priq_size(pq_start));
    TEST_ASSERT_EQUAL(num_input_files, priq_size(pq_end));

    uint32_t i;
    for (i = 0; i < num_input_files; ++i)
        input_file_destroy(&(i_files[i]));
    free(i_files);

    free(pqd_starts);
    priq_free(pq_start);

    priority pri;
    while (priq_top(pq_end, &pri) != NULL) {
        free(priq_pop(pq_end, &pri));
    }
    priq_free(pq_end);
    free(gi->data_dir);
    file_index_destroy(&(gi->file_idx));
    offset_index_destroy(&(gi->offset_idx));
    chrm_index_destroy(&(gi->chrm_idx));
    free(gi);
    rmrf(output_path_name);
}
//}}}

//{{{ void test_giggle_bulk_insert_build_leaf_levels(void)
void test_giggle_bulk_insert_build_leaf_levels(void)
{
    char *input_path_name = "../data/many/*gz";
    char *output_path_name = "tmp_test_giggle_bulk_insert_build_leaf_levels";

    struct stat st = {0};
    if (stat(output_path_name, &st) == -1) {
        mkdir(output_path_name, 0700);
    } else {
        rmrf(output_path_name);
        mkdir(output_path_name, 0700);
    }

    struct giggle_index *gi = (struct giggle_index *)
            malloc(sizeof(struct giggle_index));
    gi->data_dir = strdup(output_path_name);

    // open files
    gi->file_idx = NULL;
    struct input_file **i_files = NULL;
    uint32_t num_input_files = giggle_bulk_insert_open_files(input_path_name,
                                                             gi->data_dir, 
                                                             &i_files,
                                                             &(gi->file_idx));

    TEST_ASSERT_EQUAL(22, num_input_files);

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

    pri_queue pq_end = priq_new(num_input_files);

    giggle_bulk_insert_prime_pqs(gi,
                                 &pq_start,
                                 pqd_starts,
                                 &pq_end,
                                 i_files,
                                 num_input_files);
    
    TEST_ASSERT_EQUAL(num_input_files, priq_size(pq_start));
    TEST_ASSERT_EQUAL(num_input_files, priq_size(pq_end));

    giggle_bulk_insert_build_leaf_levels(gi,
                                         &pq_start,
                                         pqd_starts,
                                         &pq_end,
                                         i_files,
                                         num_input_files);

    TEST_ASSERT_EQUAL(21024, gi->offset_idx->index->num);
    TEST_ASSERT_EQUAL(24, gi->chrm_idx->index->num);

    TEST_ASSERT_EQUAL(0, priq_size(pq_start));
    TEST_ASSERT_EQUAL(0, priq_size(pq_end));
    
    uint32_t i;
    for (i = 0; i < num_input_files; ++i)
        input_file_destroy(&(i_files[i]));
    free(i_files);

    free(pqd_starts);
    priq_free(pq_start);
    priq_free(pq_end);
    free(gi->data_dir);
    file_index_destroy(&(gi->file_idx));
    offset_index_destroy(&(gi->offset_idx));
    chrm_index_destroy(&(gi->chrm_idx));
    free(gi);
    rmrf(output_path_name);
}
//}}}

//{{{ void test_giggle_bulk_insert_build_tree_on_leaves(void)
void test_giggle_bulk_insert_build_tree_on_leaves(void)
{
    char *input_path_name = "../data/many/*gz";
    char *output_path_name = "tmp_test_giggle_bulk_insert_build_leaf_levels";

    struct stat st = {0};
    if (stat(output_path_name, &st) == -1) {
        mkdir(output_path_name, 0700);
    } else {
        rmrf(output_path_name);
        mkdir(output_path_name, 0700);
    }

    struct giggle_index *gi = (struct giggle_index *)
            malloc(sizeof(struct giggle_index));
    gi->data_dir = strdup(output_path_name);

    // open files
    gi->file_idx = NULL;
    struct input_file **i_files = NULL;
    uint32_t num_input_files = giggle_bulk_insert_open_files(input_path_name,
                                                             gi->data_dir, 
                                                             &i_files,
                                                             &(gi->file_idx));

    TEST_ASSERT_EQUAL(22, num_input_files);

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

    pri_queue pq_end = priq_new(num_input_files);

    giggle_bulk_insert_prime_pqs(gi,
                                 &pq_start,
                                 pqd_starts,
                                 &pq_end,
                                 i_files,
                                 num_input_files);
    
    TEST_ASSERT_EQUAL(num_input_files, priq_size(pq_start));
    TEST_ASSERT_EQUAL(num_input_files, priq_size(pq_end));

    giggle_bulk_insert_build_leaf_levels(gi,
                                         &pq_start,
                                         pqd_starts,
                                         &pq_end,
                                         i_files,
                                         num_input_files);

    TEST_ASSERT_EQUAL(21024, gi->offset_idx->index->num);
    TEST_ASSERT_EQUAL(24, gi->chrm_idx->index->num);

    TEST_ASSERT_EQUAL(0, priq_size(pq_start));
    TEST_ASSERT_EQUAL(0, priq_size(pq_end));
    
    uint32_t i;
    for (i = 0; i < num_input_files; ++i)
        input_file_destroy(&(i_files[i]));
    free(i_files);
    free(pqd_starts);
    priq_free(pq_start);
    priq_free(pq_end);

    giggle_bulk_insert_build_tree_on_leaves(gi);

    free(gi->data_dir);
    file_index_destroy(&(gi->file_idx));
    offset_index_destroy(&(gi->offset_idx));
    chrm_index_destroy(&(gi->chrm_idx));
    free(gi);
    rmrf(output_path_name);
}
//}}}

//{{{ void test_giggle_bulk_insert_build_leaf_levels_few(void)
void test_giggle_bulk_insert_build_leaf_levels_few(void)
{
    ORDER = 5;
    char *input_path_name = "../data/few/*gz";
    char *output_path_name = 
            "tmp_test_giggle_bulk_insert_build_leaf_levels_small";

    struct stat st = {0};
    if (stat(output_path_name, &st) == -1) {
        mkdir(output_path_name, 0700);
    } else {
        rmrf(output_path_name);
        mkdir(output_path_name, 0700);
    }

    struct giggle_index *gi = (struct giggle_index *)
            malloc(sizeof(struct giggle_index));
    gi->data_dir = strdup(output_path_name);

    // open files
    gi->file_idx = NULL;
    struct input_file **i_files = NULL;
    uint32_t num_input_files = giggle_bulk_insert_open_files(input_path_name,
                                                             gi->data_dir, 
                                                             &i_files,
                                                             &(gi->file_idx));

    TEST_ASSERT_EQUAL(3, num_input_files);

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

    pri_queue pq_end = priq_new(num_input_files);

    giggle_bulk_insert_prime_pqs(gi,
                                 &pq_start,
                                 pqd_starts,
                                 &pq_end,
                                 i_files,
                                 num_input_files);
    
    TEST_ASSERT_EQUAL(num_input_files, priq_size(pq_start));
    TEST_ASSERT_EQUAL(num_input_files, priq_size(pq_end));

    giggle_bulk_insert_build_leaf_levels(gi,
                                         &pq_start,
                                         pqd_starts,
                                         &pq_end,
                                         i_files,
                                         num_input_files);

    TEST_ASSERT_EQUAL(16, gi->offset_idx->index->num);
    TEST_ASSERT_EQUAL(1, gi->chrm_idx->index->num);

    TEST_ASSERT_EQUAL(0, priq_size(pq_start));
    TEST_ASSERT_EQUAL(0, priq_size(pq_end));
    
    uint32_t i;
    for (i = 0; i < num_input_files; ++i)
        input_file_destroy(&(i_files[i]));
    free(i_files);

    free(pqd_starts);
    priq_free(pq_start);
    priq_free(pq_end);
    free(gi->data_dir);
    file_index_destroy(&(gi->file_idx));
    offset_index_destroy(&(gi->offset_idx));
    chrm_index_destroy(&(gi->chrm_idx));
    free(gi);


    char *index_file_name = NULL, *data_file_name = NULL;
    ret = asprintf(&index_file_name,
                   "%s/%s0.idx",
                   output_path_name,
                   CACHE_FILE_NAME_PREFIX);
     ret = asprintf(&data_file_name,
                   "%s/%s0.dat",
                   output_path_name,
                   CACHE_FILE_NAME_PREFIX);
    
    struct disk_store *ds = disk_store_load(NULL,
                                            index_file_name,
                                            NULL,
                                            data_file_name);

    TEST_ASSERT_EQUAL(12, ds->num);

    uint32_t A_keys[6][5] =  {{1,2,5,10,16},
                              {21,25,32,36,40},
                              {41,45,51,52,56},
                              {60,61,65,71,72},
                              {76,80,81,85,91},
                              {92,96,101,106,301}};
    uint32_t A_starts_ens[6][5][2] = { {{2,0}, {3,0}, {4,0}, {5,1}, {5,2}},
                                       {{1,0}, {2,0}, {2,1}, {2,2}, {3,2}},
                                       {{1,0}, {2,0}, {2,1}, {2,2}, {2,3}},
                                       {{1,0}, {2,0}, {3,0}, {3,1}, {3,2}},
                                       {{0,1}, {1,1}, {2,1}, {3,1}, {3,2}},
                                       {{0,1}, {0,2}, {0,3}, {0,4}, {0,5}} };

    uint32_t A_leaf_sizes[6][3] = {{0,5,2}, {3,3,2}, {4,2,3}, {3,3,2} ,
                                   {4,3,2,}, {5,0,5}};

    uint64_t A_leading_vals_1[3] = {0,1,2};
    uint64_t A_leading_vals_2[4] = {0,1,2,7};
    uint64_t A_leading_vals_3[3] = {0,1,2};
    uint64_t A_leading_vals_4[4] = {0,1,2,12};
    uint64_t A_leading_vals_5[5] = {0,1,2,14,15};
    uint64_t *A_leading_vals[6] = {NULL,
                                   A_leading_vals_1,
                                   A_leading_vals_2,
                                   A_leading_vals_3,
                                   A_leading_vals_4,
                                   A_leading_vals_5};

    uint32_t id = 1, leaf_i = 0;

    while (id != 0) {
        uint64_t size;
        void *v = disk_store_get(ds, id - 1, &size);
        struct bpt_node *bpn_in;
        uint64_t  deserialized_size = bpt_node_deserialize(v,
                                                          size,
                                                          (void **)&bpn_in); 
        free(v);
        struct leaf_data *ld_in;
        v = disk_store_get(ds, BPT_POINTERS_BLOCK(bpn_in) - 1, &size);
        deserialized_size = leaf_data_deserialize(v,
                                                  size,
                                                  (void **)&ld_in);

        TEST_ASSERT_EQUAL(id, BPT_ID(bpn_in));
        TEST_ASSERT_EQUAL(id+1, BPT_POINTERS_BLOCK(bpn_in));
        if(id < 11)
            TEST_ASSERT_EQUAL(id+2, BPT_NEXT(bpn_in));
        else
            TEST_ASSERT_EQUAL(0, BPT_NEXT(bpn_in));


        TEST_ASSERT_EQUAL(5, BPT_NUM_KEYS(bpn_in));

        TEST_ASSERT_EQUAL(A_leaf_sizes[leaf_i][0], ld_in->num_leading);
        TEST_ASSERT_EQUAL(A_leaf_sizes[leaf_i][1], ld_in->num_starts);
        TEST_ASSERT_EQUAL(A_leaf_sizes[leaf_i][2], ld_in->num_ends);

        for (i = 0; i < ld_in->num_leading; ++i)
            TEST_ASSERT_EQUAL(A_leading_vals[leaf_i][i], ld_in->leading[i]);

        for (i = 0; i < BPT_NUM_KEYS(bpn_in); ++i) {
            TEST_ASSERT_EQUAL(A_keys[leaf_i][i], BPT_KEYS(bpn_in)[i]);

            TEST_ASSERT_EQUAL(A_starts_ens[leaf_i][i][0], 
                              ld_in->starts_pointers[i]);

            TEST_ASSERT_EQUAL(A_starts_ens[leaf_i][i][1],
                              ld_in->ends_pointers[i]);
        }

        /*
        fprintf(stderr," L:");
        for (i = 0; i < ld_in->num_leading; ++i)
            fprintf(stderr,"%u ", ld_in->leading[i]);
        fprintf(stderr," S:");
        for (i = 0; i < ld_in->num_starts; ++i)
            fprintf(stderr,"%u ", ld_in->starts[i]);
        fprintf(stderr," E:");
        for (i = 0; i < ld_in->num_ends; ++i)
            fprintf(stderr,"%u ", ld_in->ends[i]);
        fprintf(stderr,
                "\t%u -> %u\n",
                BPT_POINTERS_BLOCK(bpn_in),
                BPT_NEXT(bpn_in));
        */

        id = BPT_NEXT(bpn_in);
        leaf_i += 1;
        free(v);
        free(bpn_in->data);
        free(bpn_in);
        leaf_data_free_mem((void **)&ld_in);
    }
    rmrf(output_path_name);
}
//}}}

//{{{ void test_giggle_bulk_insert_build_tree_on_leaves_few(void)
void test_giggle_bulk_insert_build_tree_on_leaves_few(void)
{
    ORDER = 5;
    char *input_path_name = "../data/few/*gz";
    char *output_path_name = 
            "tmp_test_giggle_bulk_insert_build_leaf_levels_small";

    struct stat st = {0};
    if (stat(output_path_name, &st) == -1) {
        mkdir(output_path_name, 0700);
    } else {
        rmrf(output_path_name);
        mkdir(output_path_name, 0700);
    }

    struct giggle_index *gi = (struct giggle_index *)
            malloc(sizeof(struct giggle_index));
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

    pri_queue pq_end = priq_new(num_input_files);

    giggle_bulk_insert_prime_pqs(gi,
                                 &pq_start,
                                 pqd_starts,
                                 &pq_end,
                                 i_files,
                                 num_input_files);
    
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

    uint32_t root_id = gi->root_ids[0];

    giggle_index_destroy(&gi);


    char *index_file_name = NULL, *data_file_name = NULL;
    ret = asprintf(&index_file_name,
                   "%s/%s0.idx",
                   output_path_name,
                   CACHE_FILE_NAME_PREFIX);
     ret = asprintf(&data_file_name,
                   "%s/%s0.dat",
                   output_path_name,
                   CACHE_FILE_NAME_PREFIX);
    
    struct disk_store *ds = disk_store_load(NULL,
                                            index_file_name,
                                            NULL,
                                            data_file_name);

    uint64_t size;
    void *v = disk_store_get(ds, root_id - 1, &size);
    struct bpt_node *bpn_in;
    uint64_t  deserialized_size = bpt_node_deserialize(v,
                                                      size,
                                                      (void **)&bpn_in); 
    free(v);

    TEST_ASSERT_EQUAL(5, BPT_NUM_KEYS(bpn_in));
    uint32_t A_keys[5] = {21, 41, 60, 76, 92};
    uint32_t A_pointers[6] = {1, 3, 5, 7, 9, 11};

    for (i = 0; i < 5; ++i)
        TEST_ASSERT_EQUAL(A_keys[i], BPT_KEYS(bpn_in)[i]);

    for (i = 0; i < 6; ++i)
        TEST_ASSERT_EQUAL(A_pointers[i], BPT_POINTERS(bpn_in)[i]);

    free(bpn_in->data);
    free(bpn_in);
    disk_store_destroy(&ds);
    rmrf(output_path_name);
}
//}}}

//{{{void test_giggle_bulk_insert_few(void)
void test_giggle_bulk_insert_few(void)
{
    ORDER = 5;
    char *input_path_name = "../data/few/*gz";
    char *output_path_name = 
            "tmp_test_giggle_bulk_insert_build_leaf_levels_small";

    uint64_t num_intervals = giggle_bulk_insert(input_path_name,
                                                output_path_name,
                                                1);

    struct giggle_index *gi =
                giggle_load(output_path_name,
                            uint64_t_ll_giggle_set_data_handler);
    giggle_data_handler.giggle_collect_intersection =
            giggle_collect_intersection_data_in_block;
    giggle_data_handler.map_intersection_to_offset_list =
            leaf_data_map_intersection_to_offset_list;

    TEST_ASSERT_EQUAL(1, gi->len);
    TEST_ASSERT_EQUAL(1, gi->num);
    TEST_ASSERT_EQUAL(3, gi->file_idx->index->num);
    TEST_ASSERT_EQUAL(1, gi->chrm_idx->index->num);
    TEST_ASSERT_EQUAL(16, gi->offset_idx->index->num);

    struct giggle_query_result *gqr = NULL;
    gqr = giggle_query(gi, "1", 1, 25, gqr);

    TEST_ASSERT_EQUAL(7, gqr->num_hits);

    uint32_t A_query_result_0[3] = {0,3,5};
    uint32_t A_query_result_1[2] = {1,6};
    uint32_t A_query_result_2[2] = {2,4};
    uint32_t *A_query_result[3] = {A_query_result_0,
                                   A_query_result_1,
                                   A_query_result_2};

    uint32_t A_query_len[3] = {3, 2, 2};
    uint32_t i;
    for(i = 0; i < gqr->num_files; i++) {
        TEST_ASSERT_EQUAL(A_query_len[i], giggle_get_query_len(gqr, i));
        struct giggle_query_iter *gqi = giggle_get_query_itr(gqr, i);
        TEST_ASSERT_EQUAL(A_query_len[i], gqi->num);
        giggle_iter_destroy(&gqi);
    }

    giggle_query_result_destroy(&gqr);
    giggle_index_destroy(&gi);
    cache.destroy();
    rmrf(output_path_name);
}
//}}}

//{{{ void test_giggle_bulk_insert_build_tree_on_leaves_few_order_4(void)
void test_giggle_bulk_insert_build_tree_on_leaves_few_order_4(void)
{
    ORDER = 4;
    char *input_path_name = "../data/few/*gz";
    char *output_path_name = 
            "tmp_test_giggle_bulk_insert_build_tree_on_leaves_few_order_4";

    struct stat st = {0};
    if (stat(output_path_name, &st) == -1) {
        mkdir(output_path_name, 0700);
    } else {
        rmrf(output_path_name);
        mkdir(output_path_name, 0700);
    }

    struct giggle_index *gi = (struct giggle_index *)
            malloc(sizeof(struct giggle_index));
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

    pri_queue pq_end = priq_new(num_input_files);

    giggle_bulk_insert_prime_pqs(gi,
                                 &pq_start,
                                 pqd_starts,
                                 &pq_end,
                                 i_files,
                                 num_input_files);
    
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


    uint32_t root_id = gi->root_ids[0];

    giggle_index_destroy(&gi);


    char *index_file_name = NULL, *data_file_name = NULL;
    ret = asprintf(&index_file_name,
                   "%s/%s0.idx",
                   output_path_name,
                   CACHE_FILE_NAME_PREFIX);
     ret = asprintf(&data_file_name,
                   "%s/%s0.dat",
                   output_path_name,
                   CACHE_FILE_NAME_PREFIX);
    
    struct disk_store *ds = disk_store_load(NULL,
                                            index_file_name,
                                            NULL,
                                            data_file_name);

    uint64_t size;
    void *v = disk_store_get(ds, root_id - 1, &size);
    struct bpt_node *bpn_in;
    uint64_t  deserialized_size = bpt_node_deserialize(v,
                                                      size,
                                                      (void **)&bpn_in); 
    free(v);

    TEST_ASSERT_EQUAL(1, BPT_NUM_KEYS(bpn_in));
    TEST_ASSERT_EQUAL(76, BPT_KEYS(bpn_in)[0]);
    TEST_ASSERT_EQUAL(17, BPT_POINTERS(bpn_in)[0]);
    TEST_ASSERT_EQUAL(18, BPT_POINTERS(bpn_in)[1]);

    uint32_t left_id = BPT_POINTERS(bpn_in)[0];
    uint32_t right_id = BPT_POINTERS(bpn_in)[1];

    free(bpn_in->data);
    free(bpn_in);

    v = disk_store_get(ds, left_id - 1, &size);
    deserialized_size = bpt_node_deserialize(v,
                                             size,
                                             (void **)&bpn_in); 
    free(v);

    TEST_ASSERT_EQUAL(4, BPT_NUM_KEYS(bpn_in));
    uint32_t A_left_keys[4] = {16,36,51,61};
    for (i = 0; i < 4; ++i)
        TEST_ASSERT_EQUAL(A_left_keys[i], BPT_KEYS(bpn_in)[i]);

    uint32_t A_left_pointers[5] = {1,3,5,7,9};
    for (i = 0; i < 5; ++i)
        TEST_ASSERT_EQUAL(A_left_pointers[i], BPT_POINTERS(bpn_in)[i]);

    free(bpn_in->data);
    free(bpn_in);

    v = disk_store_get(ds, right_id - 1, &size);
    deserialized_size = bpt_node_deserialize(v,
                                             size,
                                             (void **)&bpn_in); 
    free(v);

    TEST_ASSERT_EQUAL(3, BPT_NUM_KEYS(bpn_in));
    uint32_t A_right_keys[3] = {76,91,106};
    for (i = 0; i < 3; ++i)
        TEST_ASSERT_EQUAL(A_right_keys[i], BPT_KEYS(bpn_in)[i]);

    uint32_t A_right_pointers[4] = {0,11,13,15};
    for (i = 0; i < 4; ++i)
        TEST_ASSERT_EQUAL(A_right_pointers[i], BPT_POINTERS(bpn_in)[i]);

    free(bpn_in->data);
    free(bpn_in);

    disk_store_destroy(&ds);
    rmrf(output_path_name);
}
//}}}

//{{{ void test_giggle_bulk_insert_build_tree_on_leaves_few_order_3(void)
void test_giggle_bulk_insert_build_tree_on_leaves_few_order_3(void)
{
    ORDER = 3;
    char *input_path_name = "../data/few/*gz";
    char *output_path_name = 
            "tmp_test_giggle_bulk_insert_build_tree_on_leaves_few_order_4";

    struct stat st = {0};
    if (stat(output_path_name, &st) == -1) {
        mkdir(output_path_name, 0700);
    } else {
        rmrf(output_path_name);
        mkdir(output_path_name, 0700);
    }

    struct giggle_index *gi = (struct giggle_index *)
            malloc(sizeof(struct giggle_index));
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

    pri_queue pq_end = priq_new(num_input_files);

    giggle_bulk_insert_prime_pqs(gi,
                                 &pq_start,
                                 pqd_starts,
                                 &pq_end,
                                 i_files,
                                 num_input_files);
    
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

    uint32_t root_id = gi->root_ids[0];

    TEST_ASSERT_EQUAL(24, root_id);

    giggle_index_destroy(&gi);

    char *index_file_name = NULL, *data_file_name = NULL;
    ret = asprintf(&index_file_name,
                   "%s/%s0.idx",
                   output_path_name,
                   CACHE_FILE_NAME_PREFIX);
     ret = asprintf(&data_file_name,
                   "%s/%s0.dat",
                   output_path_name,
                   CACHE_FILE_NAME_PREFIX);
    
    struct disk_store *ds = disk_store_load(NULL,
                                            index_file_name,
                                            NULL,
                                            data_file_name);

    uint64_t size;
    void *v = disk_store_get(ds, root_id - 1, &size);
    struct bpt_node *bpn_in;
    uint64_t  deserialized_size = bpt_node_deserialize(v,
                                                      size,
                                                      (void **)&bpn_in); 
    free(v);

    TEST_ASSERT_EQUAL(2, BPT_NUM_KEYS(bpn_in));
    TEST_ASSERT_EQUAL(51, BPT_KEYS(bpn_in)[0]);
    TEST_ASSERT_EQUAL(80, BPT_KEYS(bpn_in)[1]);
    TEST_ASSERT_EQUAL(21, BPT_POINTERS(bpn_in)[0]);
    TEST_ASSERT_EQUAL(22, BPT_POINTERS(bpn_in)[1]);
    TEST_ASSERT_EQUAL(23, BPT_POINTERS(bpn_in)[2]);

    disk_store_destroy(&ds);
    rmrf(output_path_name);
}
//}}}

//{{{ void test_giggle_bulk_insert_build_tree_on_leaves_few_order_4(void)
void test_giggle_bulk_insert_build_tree_on_leaves_few_order_2(void)
{
    ORDER = 2;
    char *input_path_name = "../data/few/*gz";
    char *output_path_name = 
            "tmp_test_giggle_bulk_insert_build_tree_on_leaves_few_order_2";

    struct stat st = {0};
    if (stat(output_path_name, &st) == -1) {
        mkdir(output_path_name, 0700);
    } else {
        rmrf(output_path_name);
        mkdir(output_path_name, 0700);
    }

    struct giggle_index *gi = (struct giggle_index *)
            malloc(sizeof(struct giggle_index));
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

    pri_queue pq_end = priq_new(num_input_files);

    giggle_bulk_insert_prime_pqs(gi,
                                 &pq_start,
                                 pqd_starts,
                                 &pq_end,
                                 i_files,
                                 num_input_files);
    
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

    uint32_t root_id = gi->root_ids[0];

    TEST_ASSERT_EQUAL(41, root_id);

    giggle_index_destroy(&gi);

    char *index_file_name = NULL, *data_file_name = NULL;
    ret = asprintf(&index_file_name,
                   "%s/%s0.idx",
                   output_path_name,
                   CACHE_FILE_NAME_PREFIX);
     ret = asprintf(&data_file_name,
                   "%s/%s0.dat",
                   output_path_name,
                   CACHE_FILE_NAME_PREFIX);
    
    struct disk_store *ds = disk_store_load(NULL,
                                            index_file_name,
                                            NULL,
                                            data_file_name);

    uint64_t size;
    void *v = disk_store_get(ds, root_id - 1, &size);
    struct bpt_node *bpn_in;
    uint64_t  deserialized_size = bpt_node_deserialize(v,
                                                      size,
                                                      (void **)&bpn_in); 
    free(v);

    uint32_t root_keys[2] = {56, 81};
    uint32_t root_values[3] = {38, 39, 40};
    TEST_ASSERT_EQUAL(2, BPT_NUM_KEYS(bpn_in));
    TEST_ASSERT_EQUAL(root_keys[0], BPT_KEYS(bpn_in)[0]);
    TEST_ASSERT_EQUAL(root_keys[1], BPT_KEYS(bpn_in)[1]);
    TEST_ASSERT_EQUAL(root_values[0], BPT_POINTERS(bpn_in)[0]);
    TEST_ASSERT_EQUAL(root_values[1], BPT_POINTERS(bpn_in)[1]);
    TEST_ASSERT_EQUAL(root_values[2], BPT_POINTERS(bpn_in)[2]);

    free(bpn_in->data);
    free(bpn_in);

    v = disk_store_get(ds, root_values[0] - 1, &size);
    deserialized_size = bpt_node_deserialize(v,
                                             size,
                                             (void **)&bpn_in); 
    free(v);

    uint32_t left_keys[2] = {25, 41};
    uint32_t left_values[3] = {31, 32, 33};
    TEST_ASSERT_EQUAL(2, BPT_NUM_KEYS(bpn_in));
    TEST_ASSERT_EQUAL(left_keys[0], BPT_KEYS(bpn_in)[0]);
    TEST_ASSERT_EQUAL(left_keys[1], BPT_KEYS(bpn_in)[1]);
    TEST_ASSERT_EQUAL(left_values[0], BPT_POINTERS(bpn_in)[0]);
    TEST_ASSERT_EQUAL(left_values[1], BPT_POINTERS(bpn_in)[1]);
    TEST_ASSERT_EQUAL(left_values[2], BPT_POINTERS(bpn_in)[2]);

    free(bpn_in->data);
    free(bpn_in);

    v = disk_store_get(ds, root_values[1] - 1, &size);
    deserialized_size = bpt_node_deserialize(v,
                                             size,
                                             (void **)&bpn_in); 
    free(v);


    uint32_t middle_keys[2] = {56, 71};
    uint32_t middle_values[3] = {0, 34, 35};
    TEST_ASSERT_EQUAL(2, BPT_NUM_KEYS(bpn_in));
    TEST_ASSERT_EQUAL(middle_keys[0], BPT_KEYS(bpn_in)[0]);
    TEST_ASSERT_EQUAL(middle_keys[1], BPT_KEYS(bpn_in)[1]);
    TEST_ASSERT_EQUAL(middle_values[0], BPT_POINTERS(bpn_in)[0]);
    TEST_ASSERT_EQUAL(middle_values[1], BPT_POINTERS(bpn_in)[1]);
    TEST_ASSERT_EQUAL(middle_values[2], BPT_POINTERS(bpn_in)[2]);

    free(bpn_in->data);
    free(bpn_in);

    v = disk_store_get(ds, root_values[2] - 1, &size);
    deserialized_size = bpt_node_deserialize(v,
                                             size,
                                             (void **)&bpn_in); 
    free(v);


    uint32_t right_keys[2] = {81, 96};
    uint32_t right_values[3] = {0, 36, 37};
    TEST_ASSERT_EQUAL(2, BPT_NUM_KEYS(bpn_in));
    TEST_ASSERT_EQUAL(right_keys[0], BPT_KEYS(bpn_in)[0]);
    TEST_ASSERT_EQUAL(right_keys[1], BPT_KEYS(bpn_in)[1]);
    TEST_ASSERT_EQUAL(right_values[0], BPT_POINTERS(bpn_in)[0]);
    TEST_ASSERT_EQUAL(right_values[1], BPT_POINTERS(bpn_in)[1]);
    TEST_ASSERT_EQUAL(right_values[2], BPT_POINTERS(bpn_in)[2]);

    free(bpn_in->data);
    free(bpn_in);

    disk_store_destroy(&ds);
    rmrf(output_path_name);
}
//}}}

//{{{void test_giggle_bulk_insert(void)
void test_giggle_bulk_insert_fews(void)
{
    char *input_path_name = "../data/few/*gz";
    char *output_path_name = "test_giggle_bulk_insert_fews";

    uint32_t order;
    for (order = 2; order <= 100; ++order) {
        ORDER = order;
        uint64_t num_intervals = giggle_bulk_insert(input_path_name,
                                                    output_path_name,
                                                    1);
        struct giggle_index *gi =
                    giggle_load(output_path_name,
                                uint64_t_ll_giggle_set_data_handler);
        giggle_data_handler.giggle_collect_intersection =
                giggle_collect_intersection_data_in_block;
        giggle_data_handler.map_intersection_to_offset_list =
                leaf_data_map_intersection_to_offset_list;

        struct giggle_query_result *gqr = NULL;
        gqr = giggle_query(gi, "1", 20, 50, gqr);

        TEST_ASSERT_EQUAL(8, gqr->num_hits);

        giggle_query_result_destroy(&gqr);
        giggle_index_destroy(&gi);
        cache.destroy();
    }
    rmrf(output_path_name);
}
//}}}

//{{{void test_giggle_bulk_insert(void)
void test_giggle_bulk_insert(void)
{
    char *input_path_name = "../data/many/*gz";
    char *output_path_name = "tmp_test_giggle_bulk_insert";

    uint64_t num_intervals = giggle_bulk_insert(input_path_name,
                                                output_path_name,
                                                1);
    TEST_ASSERT_EQUAL(21024, num_intervals);

    struct giggle_index *gi =
                giggle_load(output_path_name,
                            uint64_t_ll_giggle_set_data_handler);
    giggle_data_handler.giggle_collect_intersection =
            giggle_collect_intersection_data_in_block;
    giggle_data_handler.map_intersection_to_offset_list =
            leaf_data_map_intersection_to_offset_list;

    TEST_ASSERT_EQUAL(24, gi->len);
    TEST_ASSERT_EQUAL(24, gi->num);
    TEST_ASSERT_EQUAL(22, gi->file_idx->index->num);
    TEST_ASSERT_EQUAL(24, gi->chrm_idx->index->num);
    TEST_ASSERT_EQUAL(21024, gi->offset_idx->index->num);

    struct giggle_query_result *gqr = NULL;
    gqr = giggle_query(gi, "1", 5000000, 6200000, gqr);

    //fprintf(stderr, "%u\n", gqr->num_hits);

    giggle_index_destroy(&gi);
    rmrf(output_path_name);
}
//}}}
