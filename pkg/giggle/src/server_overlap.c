#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <sys/types.h>
#include <sys/select.h>
#include <sys/socket.h>
#include <microhttpd.h>
#include <regex.h>
#include <err.h>

#include "util.h"
#include "giggle_index.h"
#include "file_read.h"
#include "ll.h"

struct region
{
    char *chrm;
    uint32_t start, end;
    char *file_patterns_to_be_printed;
    uint32_t full;
};

struct request
{
    int success;
    union {
        struct region *reg;
    } data;
    enum {
        REQUEST_REGION,
        REQUEST_DATA
    } type;
};

struct args {
    struct giggle_index *gi;
    char *data_def;
};

//{{{ struct region *init_region()
struct region *init_region()
{
    struct region *r = (struct region *)malloc(sizeof(struct region));
    r->chrm = NULL;
    r->file_patterns_to_be_printed = NULL;
    r->full = 0;

    return r;
}
//}}}

//{{{ uint32_t parse_file_patterns(char *file_patterns_to_be_printed,
uint32_t parse_file_patterns(char *file_patterns_to_be_printed,
                             uint32_t *num_file_patterns,
                             regex_t **regexs,
                             char ***file_patterns)
{
    *num_file_patterns = 0;
    int s = 0, e = 0;
    while (scan_s(file_patterns_to_be_printed,
                  strlen(file_patterns_to_be_printed),
                  &s,
                  &e,
                  ',') >= 0) {
        *num_file_patterns = (*num_file_patterns) + 1;
        s = e + 1;
    }

    if (*num_file_patterns == 0) {
        fprintf(stderr, "No file patterns detected.\n");
        return 0;
    }

    *regexs = (regex_t *)
            malloc(*num_file_patterns * sizeof(regex_t));

    *file_patterns = (char **)
            malloc(*num_file_patterns * sizeof(char *));
    uint32_t i = 0;
    s = 0;
    e = 0;
    while (scan_s(file_patterns_to_be_printed,
                  strlen(file_patterns_to_be_printed),
                  &s,
                  &e,
                  ',') >= 0) {
        (*file_patterns)[i] = strndup(file_patterns_to_be_printed + s, e-s);
        fprintf(stderr, "%s\n", (*file_patterns)[i]);
        int r = regcomp(&((*regexs)[i]), (*file_patterns)[i], 0);
        if (r != 0) {
            fprintf(stderr, 
                    "Could not compile regex '%s'",
                    (*file_patterns)[i]);
            return 0;
        }
        i += 1;
        s = e + 1;
    }

    return 1;
}
//}}}


//{{{ int scan_url_vals(void *cls,
int scan_url_vals(void *cls,
                  enum MHD_ValueKind kind,
                  const char *key,
                  const char *value)
{
    struct request *r = (struct request *)cls;

    if ((key != NULL)) {
        if (strncmp("data", key, 6) == 0) {
            fprintf(stderr, "scan_url_vals: data\n");
            r->type = REQUEST_DATA;
            r->success = 1;
        } else if (strncmp("region", key, 6) == 0) {
            fprintf(stderr, "scan_url_vals: region\n");
            r->type = REQUEST_REGION;
            if (value != NULL) {
                if (r->data.reg == NULL)
                    r->data.reg = init_region();
                if (parse_region((char *)value,
                                 &(r->data.reg->chrm), 
                                 &(r->data.reg->start),
                                 &(r->data.reg->end)) == 0) {
                    r->success = 1;
                }
            } 
        } else if (strncmp("files", key, 5) == 0) {
            fprintf(stderr, "scan_url_vals: file\n");
            if (r->data.reg == NULL)
                r->data.reg = init_region();
            if (value != NULL) {

                //fprintf(stderr, "files:%s\n", value);
                r->data.reg->file_patterns_to_be_printed = 
                    (char *)malloc((strlen(value)+1) * sizeof(char));
                strcpy(r->data.reg->file_patterns_to_be_printed, value);
            }
        } else if (strncmp("full", key, 5) == 0) {
            fprintf(stderr, "scan_url_vals: full\n");
            if (r->data.reg == NULL)
                r->data.reg = init_region();
            r->data.reg->full = 1;
        }
    }

    return MHD_YES;
}
//}}}

int answer_to_connection(void *cls,
                         struct MHD_Connection *connection, 
                         const char *url, 
                         const char *method, const char *version, 
                         const char *upload_data, 
                         size_t *upload_data_size, void **con_cls)
{
    struct args *arg = (struct args*) cls;
    struct MHD_Response *response;
    int ret = 0;

    struct request r;
    r.success = 0;
    r.data.reg = NULL;

    uint32_t file_patterns_set_is_set = 0;
    uint32_t num_file_patterns = 0;
    regex_t *regexs = NULL;
    char **file_patterns = NULL;

    int num_vals = MHD_get_connection_values(connection,
                                             MHD_GET_ARGUMENT_KIND,
                                             scan_url_vals,
                                             &r);

    if (r.success == 1) {
        if (r.type == REQUEST_DATA) {
            size_t l = strlen(arg->data_def);
            response = MHD_create_response_from_buffer(l,
                                                       (void*) (arg->data_def),
                                                       MHD_RESPMEM_PERSISTENT);
	    MHD_add_response_header(response,
                                    MHD_HTTP_HEADER_ACCESS_CONTROL_ALLOW_ORIGIN,
                                    "*");
            ret = MHD_queue_response (connection, MHD_HTTP_OK, response);
            MHD_destroy_response (response);
        } else if (r.type == REQUEST_REGION) {
            struct region *reg = r.data.reg;
            fprintf(stderr,
                    "c:%s s:%u e:%u\n",
                    reg->chrm,
                    reg->start,
                    reg->end);

            if (reg->file_patterns_to_be_printed != NULL) {
                fprintf(stderr,
                        "files:%s\n",
                        reg->file_patterns_to_be_printed);

                file_patterns_set_is_set = 
                        parse_file_patterns(reg->file_patterns_to_be_printed,
                                            &num_file_patterns,
                                            &regexs,
                                            &file_patterns);
            }

            char *page = NULL, *tmp_page = NULL;;
            struct giggle_index *gi = arg->gi;
            struct giggle_query_result *gqr = giggle_query(gi,
                                                           reg->chrm,
                                                           reg->start,
                                                           reg->end,
                                                           NULL);

            fprintf(stderr,
                    "file_patterns_set_is_set:%u\n",
                    file_patterns_set_is_set);

            uint32_t i, printed_i = 0;
            for(i = 0; i < gqr->num_files; i++) {
                struct file_data *fd = 
                    file_index_get(gi->file_idx, i);
                    //(struct file_data *)unordered_list_get(gi->file_index, i); 

                //fprintf(stderr,
                        //"regexs:%p\t"
                        //"file_patterns:%p\t"
                        //"num_file_patterns:%u\n",
                        //regexs,
                        //file_patterns,
                        //num_file_patterns);

                if (test_pattern_match(gi,
                                       regexs,
                                       file_patterns,
                                       num_file_patterns,
                                       i,
                                       file_patterns_set_is_set)) {
                    if (printed_i == 0) {
                        asprintf(&tmp_page,
                                 "#%s\t"
                                 "%u\t"
                                 "%u\n",
                                 fd->file_name,
                                 fd->num_intervals,
                                 giggle_get_query_len(gqr, i));
                    } else {
                        asprintf(&tmp_page,
                                 "%s"
                                 "#%s\t"
                                 "%u\t"
                                 "%u\n",
                                 page,
                                 fd->file_name,
                                 fd->num_intervals,
                                 giggle_get_query_len(gqr, i));
                    }

                    free(page);
                    page = tmp_page;

                    if ( (reg->full == 1) && 
                         (giggle_get_query_len(gqr, i) > 0 )){

                        char *result = NULL;

                        struct giggle_query_iter *gqi =
                                giggle_get_query_itr(gqr, i);
                        while (giggle_query_next(gqi, &result) == 0) {
                            //printf("%s\n", result);
                            asprintf(&tmp_page, "%s%s\n", page, result);
                            free(page);
                            page = tmp_page;
                        }

                        giggle_iter_destroy(&gqi);

                        /*
                        if (result != NULL) {
                            free(result);
                            result = NULL;
                        }
                        */
                    }

                    printed_i += 1;

                }
            }

	    if (page == NULL) {
	        page = (char *)malloc(sizeof (char));
                page[0] = '\0';
	    }

            response = MHD_create_response_from_buffer(strlen (page),
                                                       (void*) page,
                                                       MHD_RESPMEM_MUST_FREE);
	    MHD_add_response_header(response,
                                    MHD_HTTP_HEADER_ACCESS_CONTROL_ALLOW_ORIGIN,
                                    "*");

            ret = MHD_queue_response (connection, MHD_HTTP_OK, response);
            MHD_destroy_response (response);
            giggle_query_result_destroy(&gqr);
            free(reg->chrm);
            reg->chrm = NULL;
        }
    }

    if (regexs != NULL)
        regfree(regexs);

    return ret;
}

int main(int argc, char **argv)
{
    if (argc != 5) {
        fprintf(stderr,
                "usage:\t%s <num threads> <index dir> "
                "<data definition> <port>\n",
                argv[0]);
        return 0;
    }

    uint32_t NUMBER_OF_THREADS = atoi(argv[1]);
    char *index_dir_name = argv[2]; 
    char *data_def_file_name = argv[3]; 
    uint32_t port = atoi(argv[4]); 

    FILE *fp = fopen(data_def_file_name, "r");
    if (!fp)
        err(1, "Could not open header file '%s'", data_def_file_name);
    fseek(fp, 0, SEEK_END);
    long data_def_size = ftell(fp);
    fseek(fp, 0, SEEK_SET);
    char *data_def = (char *)malloc((data_def_size + 1)*sizeof(char));
    fread(data_def, data_def_size, 1, fp);
    data_def[data_def_size] = '\0';
    fclose(fp);

    struct args *arg = (struct args *)malloc(sizeof(struct args));

    arg->data_def = data_def;

    arg->gi = giggle_load(index_dir_name,
                          uint64_t_ll_giggle_set_data_handler);

#if BLOCK_STORE
    giggle_data_handler.giggle_collect_intersection =
            giggle_collect_intersection_data_in_block;

    giggle_data_handler.map_intersection_to_offset_list =
            leaf_data_map_intersection_to_offset_list;
#endif

    struct MHD_Daemon *daemon;

    daemon = MHD_start_daemon(MHD_USE_SELECT_INTERNALLY,
                              port,
                              NULL,
                              NULL, 
                              &answer_to_connection,
                              arg,
                              MHD_OPTION_THREAD_POOL_SIZE,
                              (unsigned int) NUMBER_OF_THREADS,
                              MHD_OPTION_END);

    if (NULL == daemon) return 1;
    getchar (); 

    free(arg->data_def);
    giggle_index_destroy(&(arg->gi));
    cache.destroy();
    MHD_stop_daemon (daemon);
    return 0;
}
