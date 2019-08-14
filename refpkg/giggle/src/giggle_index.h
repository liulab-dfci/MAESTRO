#ifndef __GIGGLE_INDEX_H__
#define __GIGGLE_INDEX_H__

#include <stdint.h>
#include <htslib/khash.h>
#include "bpt.h"
#include "cache.h"
#include "leaf.h"
#include "jsw_avltree.h"
#include "pq.h"
#include "offset_index.h"

#define PROGRAM_NAME  "giggle"
#define MAJOR_VERSION "0"
#define MINOR_VERSION "6"
#define REVISION_VERSION "3"
#define BUILD_VERSION "0"
#define VERSION MAJOR_VERSION "." MINOR_VERSION "." REVISION_VERSION

char *CHRM_INDEX_FILE_NAME;
char *FILE_INDEX_FILE_NAME;
char *ROOT_IDS_FILE_NAME;
char *CACHE_FILE_NAME_PREFIX;

void c_str_store(void *v, FILE *f, char *file_name);
void *c_str_load(FILE *f, char *file_name);

struct chrm_index
{
    char *file_name;
    struct ordered_set *index;

};
struct chrm_index *chrm_index_init(uint32_t init_size, char *file_name);
struct chrm_index *chrm_index_load(char *file_name);
struct str_uint_pair *chrm_index_get(struct chrm_index *ci, char *chrm);
uint32_t chrm_index_add(struct chrm_index *ci, char *chrm);
void chrm_index_store(struct chrm_index *ci);
void chrm_index_destroy(struct chrm_index **ci);

struct file_index 
{
    char *file_name; //<! file_index file name
    struct unordered_list *index; //<! database file_data elements list
};

struct file_index *file_index_init(uint32_t init_size, char *file_name);
void file_index_destroy(struct file_index **fi);
uint32_t file_index_add(struct file_index *fi, char *file_name);
void file_index_store(struct file_index *fi);
struct file_index *file_index_load(char *file_name);
struct file_data *file_index_get(struct file_index *fi, uint32_t id);

/**
 * @brief The core GIGGLE data structure.
 *
 * Each chrom has its own tree and this struct contains the root node id within
 * each tree. The mapping between chrom strings (e.g., chr1, 1, etc.) is
 * managed by chrm_index. Details about the files contained within the index
 * are stored in file_index. The source file index (within file_index) and
 * relative position of each indexed interval are in offset_index.
 */
struct giggle_index
{
    uint32_t *root_ids; //!< list of root index for each chrom
    uint32_t len; //<! allocated size of root_ids
    uint32_t num; //<! number elements in root_ids
    struct chrm_index *chrm_idx; //<! chrom string/root_id index pairs
    struct file_index *file_idx;
    struct offset_index *offset_idx; //<! file_index/offse pair list
    char *data_dir; //<! database directory
    char *root_ids_file_name; //<! root_ids file name
};

struct giggle_query_result
{
    struct giggle_index *gi;
    uint32_t num_files;
    struct long_uint_ll **offsets;
    uint32_t num_hits;
};

struct giggle_query_result *giggle_query(struct giggle_index *gi,
                                        char *chrm,
                                        uint32_t start,
                                        uint32_t end,
                                        struct giggle_query_result *gqr);

void giggle_query_result_destroy(struct giggle_query_result **gqr);

struct giggle_query_iter
{
    struct giggle_index *gi;
    uint32_t file_id, curr, num;
    struct input_file *ipf;
    long *sorted_offsets;
    struct long_uint_pair *sorted_offset_id_pairs;
};

uint32_t giggle_get_query_len(struct giggle_query_result *gqr,
                              uint32_t file_id);
                      
struct giggle_query_iter *giggle_get_query_itr(struct giggle_query_result *gqr,
                                               uint32_t file_id);

struct giggle_query_iter *giggle_get_query_data_itr(
        struct giggle_query_result *gqr,
        uint32_t file_id);

int giggle_query_next(struct giggle_query_iter *gqi,
                      char **result);

int giggle_query_data(struct giggle_query_iter *gqi,
                      void *result);
int giggle_query_next_data(struct giggle_query_iter *gqi,
                           void **result);

void giggle_iter_destroy(struct giggle_query_iter **gqi);

/**
 * @brief
 *
 * @param domain specifies which cache domain (i.e. which tree) to use,
 * typically chrom id
 * @param root_id index of the root node within the tree
 * @param start interval start
 * @param end interval end
 * @param id interval id, typically an index into the offset_index
 *
 * @retval
 */
uint32_t giggle_insert(uint32_t domain,
                       uint32_t *root_id,
                       uint32_t start,
                       uint32_t end,
                       uint64_t id);

void *giggle_search(uint32_t domain,
                    uint32_t root_id,
                    uint32_t start,
                    uint32_t end);

struct giggle_def 
{
    struct cache_handler non_leading_cache_handler;
    struct cache_handler leading_cache_handler;
    void *(*new_non_leading)(uint32_t domain);
    void *(*new_leading)(uint32_t domain);
    void (*non_leading_SA_add_scalar)(uint32_t domain,
                                      void *non_leading,
                                      void *scalar);
    void (*non_leading_SE_add_scalar)(uint32_t domain, 
                                      void *non_leading,
                                      void *scalar);
    void (*leading_B_add_scalar)(uint32_t domain,
                                 void *leading,
                                 void *scalar);
    void (*leading_union_with_B)(uint32_t domain,
                                 void **result,
                                 void *leading);
    void (*non_leading_union_with_SA)(uint32_t domain,
                                      void **result,
                                      void *non_leading);
    void (*non_leading_union_with_SA_subtract_SE)(uint32_t domain,
                                                  void **result,
                                                  void *non_leading);
    void (*write_tree)(void *arg);
    void *(*giggle_collect_intersection)(uint32_t leaf_start_id,
                                         int pos_start_id,
                                         uint32_t leaf_end_id,
                                         int pos_end_id,
                                         uint32_t domain,
                                         void **r);
    void (*map_intersection_to_offset_list)(struct giggle_index *gi,
                                            struct giggle_query_result *gqr,
                                            void *_R);
};

struct giggle_def giggle_data_handler;

struct giggle_index *giggle_init_index(uint32_t init_size, char *offset_file_name);
void giggle_index_destroy(struct giggle_index **gi);
uint32_t giggle_get_chrm_id(struct giggle_index *gi, char *chrm);
uint32_t giggle_get_file_id(struct giggle_index *gi, char *path);

/**
 * @brief Add the intervals from a file to a GIGGLE index
 *
 * The file will be added to the GIGGLE index file_index, all of the
 * intervals to the offset_index, and giggle_insert will add the interval to
 * the tree where the starnd and end are keys and the value is that interval's
 * index into the offset_index. Stats about the file (mean interval size,
 * number of intervals) are collected for statistical analysis on query
 * results.
 *
 * @param gi an initialized GIGGLE index
 * @param file_name the file to be indexed
 *
 * @retval the total number of intervals indexed
 */
uint32_t giggle_index_file(struct giggle_index *gi,
                           char *file_name); 
/**
 * @brief A wrapper around giggle_index_file, that loops over files in a
 * directory.
 *
 * @param gi the GIGGLE index struct to be populated
 * @param path_name path to a set of interval files to be indexed
 * stored.
 * @param verbose give extra output in the indexing process
 *
 * @retval the total number of intervals indexed
 *
 * Example Usage:
 * @code
 *      struct giggle_index *gi = g
 *              giggle_init(23,
 *                          NULL,
 *                          0,
 *                          uint64_t_ll_giggle_set_data_handler);
 *      char *path_name = "../data/many/\*bed.gz";
 *      uint32_t r = giggle_index_directory(gi, path_name, 0); 
 *      giggle_index_destroy(&gi);
 *      cache.destroy();
 * @endcode
 */ 
uint32_t giggle_index_directory(struct giggle_index *gi,
                                char *path_name,
                                int verbose);
void *giggle_query_region(struct giggle_index *gi,
                          char *chrm,
                          uint32_t start,
                          uint32_t end);

/**
 * @brief Initialize a new GIGGLE index.
 *
 * @param num_chrms an estimate of how many chroms will be considered
 * @param output_dir the directory to store files, can be NULL if the database
 * will not be saved
 * @param force 1 to overwrite any existing index in output_dir, 0 to not
 * @param giggle_set_data_handler defines how to organize the data.
 * uint64_t_ll_giggle_set_data_handler is good for building the index from
 * scratch.
 *
 * @retval an initialize GIGGLE index
 * Example Usage:
 * @code
 *      struct giggle_index *gi = g
 *              giggle_init(23,
 *                          "giggle_i",
 *                          0,
 *                          uint64_t_ll_giggle_set_data_handler);
 *
 *      giggle_index_destroy(&gi);
 *      cache.destroy();
 * @endcode
 */
struct giggle_index *giggle_init(uint32_t num_chrms,
                                 char *output_dir,
                                 uint32_t force,
                                 void (*giggle_set_data_handler)());
uint32_t giggle_store(struct giggle_index *gi);

/**
 * @brief Load a previously stored GIGGLE index.
 *
 * @param data_dir the directory containing all database files
 * @param giggle_set_data_handler defines how data is stored
 *
 * @retval the GIGGLE index
 *
 * Example Usage:
 * @code
 *      
 * @endcode
 */
struct giggle_index *giggle_load(char *data_dir,
                                 void (*giggle_set_data_handler)(void));

struct cache_handler leaf_data_cache_handler;

uint32_t giggle_get_leaf_data(struct giggle_index *gi,
                              uint32_t domain,
                              uint32_t leaf_id,
                              struct leaf_data **lf);

struct cache_handler leaf_data_cache_handler;

uint32_t giggle_leaf_data_get_intersection_size(uint32_t leaf_start_id,
                                                int pos_start_id,
                                                uint32_t leaf_end_id,
                                                int pos_end_id,
                                                uint32_t domain);

uint32_t giggle_leaf_data_get_intersection(uint32_t leaf_start_id,
                                           int pos_start_id,
                                           uint32_t leaf_end_id,
                                           int pos_end_id,
                                           uint32_t domain);

void leaf_data_map_intersection_to_offset_list(struct giggle_index *gi,
                                               struct giggle_query_result *gqr,
                                               void *_R);

//void giggle_write_tree(struct giggle_index *gi);
void giggle_write_tree_cache_dump(void *giggle_index);
void giggle_write_tree_leaf_data(void *giggle_index);


void *giggle_collect_intersection_data_in_pointers(uint32_t leaf_start_id,
                                                   int pos_start_id,
                                                   uint32_t leaf_end_id,
                                                   int pos_end_id,
                                                   uint32_t domain,
                                                   void **r);

void *giggle_collect_intersection_data_in_block(uint32_t leaf_start_id,
                                                int pos_start_id,
                                                uint32_t leaf_end_id,
                                                int pos_end_id,
                                                uint32_t domain,
                                                void **r);
#if 0
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
                           uint32_t *merged_ends_num);

uint32_t giggle_merge_chrom(char *chrm_string,
                            struct giggle_index *gi_0,
                            struct indexed_list *file_index_id_map_0,
                            uint32_t gi_0_cache_name_space,
                            struct giggle_index *gi_1,
                            struct indexed_list *file_index_id_map_1,
                            uint32_t gi_1_cache_name_space,
                            struct disk_store *ds,
                            struct file_id_offset_pairs **merged_offset_index);

uint32_t giggle_merge_chrm_union(struct giggle_index *gi_0,
                                 struct giggle_index *gi_1,
                                 char ***merged_chrm_set);

uint32_t giggle_merge_add_file_index(struct giggle_index *gi,
                                     struct indexed_list *file_index_id_map,
                                     struct unordered_list *merged_file_index);
#endif

int giggle_bulk_insert_append_bpt_key(struct bpt_node *bpn,
                                      uint32_t key_val,
                                      struct disk_store *ds,
                                      jsw_avltree_t *avl,
                                      struct uint64_t_array *leading,
                                      struct uint64_t_array *starts,
                                      struct uint64_t_array *ends,
                                      struct uint32_t_array *starts_pointers,
                                      struct uint32_t_array *ends_pointers);

void giggle_bulk_insert_write_leaf_node(struct bpt_node *bpn,
                                        struct disk_store *ds,
                                        struct uint64_t_array *leading,
                                        struct uint64_t_array *starts,
                                        struct uint64_t_array *ends,
                                        struct uint32_t_array *starts_pointers,
                                        struct uint32_t_array *ends_pointers);

uint32_t giggle_bulk_insert_add_tree_level(struct disk_store *curr_ds,
                                           uint32_t curr_level_first_id,
                                           uint32_t curr_level_num_nodes,
                                           uint32_t curr_level_is_leaf,
                                           uint32_t *new_level_first_id);
uint32_t giggle_bulk_insert_open_files(char *path_name,
                                       char *output_dir_name,
                                       struct input_file ***i_files,
                                       struct file_index **file_idx);

void giggle_bulk_insert_prime_pqs(struct giggle_index *gi,
                                  pri_queue *pq_start,
                                  struct pq_data *pqd_starts,
                                  pri_queue *pq_end,
                                  struct input_file **i_files,
                                  uint32_t num_input_files);
void giggle_bulk_insert_build_leaf_levels(struct giggle_index *gi,
                                          pri_queue *pq_start,
                                          struct pq_data *pqd_starts,
                                          pri_queue *pq_end,
                                          struct input_file **i_files,
                                          uint32_t num_input_files);
void giggle_bulk_insert_build_tree_on_leaves(struct giggle_index *gi);

uint64_t giggle_bulk_insert(char *input_path_name,
                            char *output_path_name,
                            uint32_t force);

uint32_t giggle_get_indexed_files(char *index_dir_name,
                                  char ***names,
                                  uint32_t **num_intervals,
                                  double **mean_interval_sizes);

void block_store_giggle_set_data_handler();
#endif
