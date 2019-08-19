#define _GNU_SOURCE
#include <sys/types.h>
#include <sys/select.h>
#include <sys/socket.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <microhttpd.h>
#include <regex.h>
#include <err.h>
#include <sysexits.h>
#include <htslib/kstring.h>
#include <math.h>
#include <getopt.h>
#include <ctype.h>

#include "util.h"
#include "giggle_index.h"
#include "file_read.h"
#include "ll.h"
#include "kfunc.h"


#define MAX(X, Y) (((X) > (Y)) ? (X) : (Y))

//#define PORT            8888
#define POSTBUFFERSIZE  512
#define MAXCLIENTS      2

#define GET             0
#define POST            1

struct region
{
    char *chrm;
    uint32_t start, end;
};

struct query_file
{
    char *file_name;
};

struct request
{
    int success;
    union {
        struct region *reg;
        struct query_file *query;
    } data;
    enum {
        REQUEST_QUERY_FILE,
        REQUEST_REGION,
        REQUEST_DATA
    } type;
    uint32_t full;
    char *file_patterns_to_be_printed;
};

struct args {
    struct giggle_index *gi;
    char *data_def;
    char *upload_dir;
};

static unsigned int nr_of_uploading_clients = 0;

struct connection_info_struct
{
  int connectiontype;
  struct MHD_PostProcessor *postprocessor;
  FILE *fp;
  char *file_name;
  const char *answerstring;
  int answercode;
  struct args *arg;
};

//{{{ static page data
const char *askpage = "<html><body>\n\
                       Upload a file, please!<br>\n\
                       There are %u clients uploading at the moment.<br>\n\
                       <form action=\"/filepost\" method=\"post\" enctype=\"multipart/form-data\">\n\
                       <input name=\"file\" type=\"file\">\n\
                       <input type=\"submit\" value=\" Send \"></form>\n\
                       </body></html>";

const char *busypage =
  "<html><body>This server is busy, please try again later.</body></html>";

const char *completepage =
  "<html><body>The upload has been completed.</body></html>";

const char *errorpage =
  "<html><body>This doesn't seem to be right.</body></html>";
const char *servererrorpage =
  "<html><body>An internal server error has occured.</body></html>";
const char *fileexistspage =
  "<html><body>This file already exists.</body></html>";
const char *queryfile_errorpage = "QueryFileError";
//}}}

int SET_ACCESS_CONTROL_HEADER = 0;

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
    if (*regexs == NULL)
        err(1, "malloc error in parse_file_patterns().");

    *file_patterns = (char **)
            malloc(*num_file_patterns * sizeof(char *));
    if (*file_patterns == NULL)
        err(1, "malloc error in parse_file_patterns().");
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

//{{{ struct query_file *init_query_file()
struct query_file *init_query_file()
{

    struct query_file *q = (struct query_file *)
            malloc(sizeof(struct query_file));
    if (q == NULL)
        err(1, "malloc error in init_query_file().");

    q->file_name = NULL;
    return q;
}
//}}}

//{{{ struct region *init_region()
struct region *init_region()
{
    struct region *r = (struct region *)malloc(sizeof(struct region));
    if (r == NULL)
        err(1, "malloc error in init_region().");

    r->chrm = NULL;
    return r;
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
        } else if (strncmp("query", key, 6) == 0) {
            fprintf(stderr, "scan_url_vals: query\n");
            r->type = REQUEST_QUERY_FILE;
            if (value != NULL) {
                if (r->data.query == NULL)
                    r->data.query = init_query_file();

		r->data.query->file_name = strdup(value);
                r->success = 1;
            } 
        } else if (strncmp("files", key, 5) == 0) {
            fprintf(stderr, "scan_url_vals: file\n");
            r->file_patterns_to_be_printed  =
                    (char *)malloc((strlen(value)+1) * sizeof(char));
            if (r->file_patterns_to_be_printed == NULL)
                err(1, "malloc error in scan_url_vals().");

            strcpy(r->file_patterns_to_be_printed, value);
#if 0
            if (r->data.reg == NULL)
                r->data.reg = init_region();
            if (value != NULL) {

                //fprintf(stderr, "files:%s\n", value);
                r->data.reg->file_patterns_to_be_printed = 
                    (char *)malloc((strlen(value)+1) * sizeof(char));
                strcpy(r->data.reg->file_patterns_to_be_printed, value);
            }
#endif
        } else if (strncmp("full", key, 5) == 0) {
            fprintf(stderr, "scan_url_vals: full\n");
            r->full = 1;
        }
    }

    return MHD_YES;
}
//}}}

//{{{ static int send_page (struct MHD_Connection *connection,
static int send_page (struct MHD_Connection *connection,
                      const char *page,
                      int status_code,
                      const char *type)
{
fprintf(stderr, "%s\n", type);
    int ret;
    struct MHD_Response *response;

    response = MHD_create_response_from_buffer(strlen (page),
                                               (void *) page,
                                               MHD_RESPMEM_MUST_COPY);
    if (!response)
        return MHD_NO;

    if (SET_ACCESS_CONTROL_HEADER != 0 ) 
        MHD_add_response_header(response,
                                MHD_HTTP_HEADER_ACCESS_CONTROL_ALLOW_ORIGIN,
                                "*");

    MHD_add_response_header (response, MHD_HTTP_HEADER_CONTENT_TYPE, type);
    ret = MHD_queue_response (connection, status_code, response);
    MHD_destroy_response (response);

    return ret;
}
//}}}

//{{{static int iterate_post(void *coninfo_cls,
static int iterate_post(void *coninfo_cls,
                        enum MHD_ValueKind kind,
                        const char *key,
                        const char *filename,
                        const char *content_type,
                        const char *transfer_encoding,
                        const char *data,
                        uint64_t off,
                        size_t size)
{
    struct connection_info_struct *con_info = coninfo_cls;
    FILE *fp;

    int ret = asprintf(&(con_info->file_name),
                       "%s/%s",
                       con_info->arg->upload_dir,
                       filename);

    con_info->answerstring = servererrorpage;
    con_info->answercode = MHD_HTTP_INTERNAL_SERVER_ERROR;

    if (strcmp(key, "file") != 0)
        return MHD_NO;

    if (!con_info->fp) {
        // Do not overwrite a file of the same name
        /* 
        if (NULL != (fp = fopen (con_info->file_name, "rb"))) {
            fclose (fp);
            con_info->answerstring = fileexistspage;
            con_info->answercode = MHD_HTTP_FORBIDDEN;
            return MHD_NO;
        }
        */

        con_info->fp = fopen (con_info->file_name, "ab");
        if (!con_info->fp)
            return MHD_NO;
    }

    if (size > 0) {
        if (!fwrite (data, size, sizeof (char), con_info->fp))
            return MHD_NO;
    }

    con_info->answerstring = completepage;
    con_info->answercode = MHD_HTTP_OK;

    return MHD_YES;
}
//}}}

//{{{static void request_completed (void *cls,
static void request_completed (void *cls,
                               struct MHD_Connection *connection,
                               void **con_cls,
                               enum MHD_RequestTerminationCode toe)
{
  struct connection_info_struct *con_info = *con_cls;

  if (NULL == con_info)
    return;

  if (con_info->connectiontype == POST)
    {
      if (NULL != con_info->postprocessor)
        {
          MHD_destroy_post_processor (con_info->postprocessor);
          nr_of_uploading_clients--;
        }

      if (con_info->fp)
        fclose (con_info->fp);
    }

  free (con_info);
  *con_cls = NULL;
}
//}}}

//{{{static int answer_to_connection (void *cls,
static int answer_to_connection (void *cls,
                                 struct MHD_Connection *connection,
                                 const char *url,
                                 const char *method,
                                 const char *version,
                                 const char *upload_data,
                                 size_t *upload_data_size,
                                 void **con_cls)
{
    struct args *arg = (struct args*) cls;

    
    //{{{ set up the connection info if it is null
    if (NULL == *con_cls) {
        struct connection_info_struct *con_info;

        if (nr_of_uploading_clients >= MAXCLIENTS)
            return send_page(connection,
                             busypage,
                             MHD_HTTP_SERVICE_UNAVAILABLE,
                             "text/html");

        con_info = malloc (sizeof (struct connection_info_struct));
        if (con_info == NULL)
            err(1, "malloc error in answer_to_connection().");
        con_info->arg = arg;
        if (NULL == con_info)
            return MHD_NO;

        con_info->fp = NULL;

        if (0 == strcmp(method, "POST")) {
            con_info->postprocessor =
                    MHD_create_post_processor(connection,
                                              POSTBUFFERSIZE,
                                              iterate_post,
                                              (void *)con_info);

            if (NULL == con_info->postprocessor) {
                free (con_info);
                return MHD_NO;
            }

            nr_of_uploading_clients++;

            con_info->connectiontype = POST;
            con_info->answercode = MHD_HTTP_OK;
            con_info->answerstring = completepage;
        } else {
            con_info->connectiontype = GET;
        }

        *con_cls = (void *) con_info;

        return MHD_YES;
    }
    //}}}

    //{{{GET
    if (strcmp(method, "GET") == 0) {
        struct request r;
        r.success = 0;
        r.data.reg = NULL;
        r.file_patterns_to_be_printed = NULL;
        r.full = 0;

        uint32_t file_patterns_set_is_set = 0;
        uint32_t num_file_patterns = 0;
        regex_t *regexs = NULL;
        char **file_patterns = NULL;

        struct MHD_Response *response;
        int ret = 0;

        int num_vals = MHD_get_connection_values(connection,
                                                 MHD_GET_ARGUMENT_KIND,
                                                 scan_url_vals,
                                                 &r);

        if (num_vals == 0) {
            //{{{ empty request
            char buffer[1024];
            snprintf(buffer,
                     sizeof (buffer),
                     askpage,
                     nr_of_uploading_clients);
            return send_page (connection, buffer, MHD_HTTP_OK,"text/html");
            //}}}
        } else if (r.success == 1) {
            if (r.type == REQUEST_DATA) {
                //{{{ asking for data def
                size_t l = strlen(arg->data_def);

                return send_page(connection,
                                 arg->data_def,
                                 MHD_HTTP_OK,
                                 "text/json");
                //}}}
            } else if ( (r.type == REQUEST_REGION) || 
			(r.type == REQUEST_QUERY_FILE)) {

                // check/get file patterns to print
                //char *file_patterns_to_be_printed = NULL;
                uint32_t full = r.full;

                if (r.file_patterns_to_be_printed != NULL) {
                    fprintf(stderr,
                            "files:%s\n",
                            r.file_patterns_to_be_printed);

                    file_patterns_set_is_set = 
                            parse_file_patterns(r.file_patterns_to_be_printed,
                                                &num_file_patterns,
                                                &regexs,
                                                &file_patterns);

                }

                fprintf(stderr,
                        "file_patterns_set_is_set:%u\n",
                        file_patterns_set_is_set);

                struct giggle_index *gi = arg->gi;
                struct giggle_query_result *gqr = NULL;

                // used with query file
                uint32_t num_intervals = 0;
                double mean_interval_size = 0.0;
                double genome_size =  3095677412.0;

                if (r.type == REQUEST_REGION) {
                    struct region *reg = r.data.reg;

                    fprintf(stderr,
                            "c:%s s:%u e:%u\n",
                            reg->chrm,
                            reg->start,
                            reg->end);

                    gqr = giggle_query(gi,
                                       reg->chrm,
                                       reg->start,
                                       reg->end,
                                       NULL);
                } else {
                    struct query_file *q = r.data.query;

                    int chrm_len = 50;
                    char *chrm = (char *)malloc(chrm_len*sizeof(char));
                    if (chrm == NULL)
                        err(1, "malloc error in answer_to_connection().");

                    uint32_t start, end;
                    long offset;

                    char *full_path = NULL;
                    int ret = asprintf(&full_path,
                                       "%s/%s",
                                       arg->upload_dir,
                                       q->file_name);

                    struct input_file *q_f = input_file_init(full_path);
                    if (q_f == NULL) {
                        fprintf(stderr, "Error with file %s.\n", full_path);
                	return send_page (connection, queryfile_errorpage, MHD_HTTP_OK, "text/html");
                    }

                    free(full_path);

                    kstring_t line = {0, 0, NULL};
                    while ( q_f->input_file_get_next_interval(q_f, 
                                                              &chrm,
                                                              &chrm_len,
                                                              &start,
                                                              &end,
                                                              &offset,
                                                              &line) >= 0 ) {
                        gqr = giggle_query(gi,
                                           chrm,
                                           start,
                                           end,
                                           gqr);
                        num_intervals += 1;
                        mean_interval_size += end - start;
                    }

                    if (line.s != NULL)
                        free(line.s);

                    input_file_destroy(&q_f);
                    mean_interval_size = mean_interval_size/num_intervals;
                }

                uint32_t i, printed_i = 0;
                char *page = NULL, *tmp_page = NULL;;
                for(i = 0; i < gqr->num_files; i++) {
                    struct file_data *fd = 
                            file_index_get(gi->file_idx, i);
                        //(struct file_data *)
                        //unordered_list_get(gi->file_index, i); 

                    giggle_get_query_len(gqr, i);
                    if (test_pattern_match(gi,
                                           regexs,
                                           file_patterns,
                                           num_file_patterns,
                                           i,
                                           file_patterns_set_is_set)) {
                        if (printed_i == 0) {
                            int ret = asprintf(&tmp_page,
                                               "#%s\t"
                                               "%u\t"
                                               "%u\n",
                                               fd->file_name,
                                               fd->num_intervals,
                                               giggle_get_query_len(gqr, i));
                        } else {
                            int ret = asprintf(&tmp_page,
                                               "%s"
                                               "#%s\t"
                                               "%u\t"
                                               "%u\n",
                                               page,
                                               fd->file_name,
                                               fd->num_intervals,
                                               giggle_get_query_len(gqr, i));
                        }

                            fprintf(stderr,
                                    "#%s\t"
                                    "%u\t"
                                    "%u\n",
                                    fd->file_name,
                                    fd->num_intervals,
                                    giggle_get_query_len(gqr, i));


                        free(page);
                        page = tmp_page;

                        if ( (full == 1) && 
                             (giggle_get_query_len(gqr, i) > 0 )){

                            char *result = NULL;

                            struct giggle_query_iter *gqi =
                                    giggle_get_query_itr(gqr, i);
                            while (giggle_query_next(gqi, &result) == 0) {
                                int ret = asprintf(&tmp_page,
                                                   "%s%s\n",
                                                   page,
                                                   result);
                                free(page);
                                page = tmp_page;
                            }
                            giggle_iter_destroy(&gqi);
                        }
                        printed_i += 1;
                    }
                }

                if (page == NULL) {
                    page = (char *)malloc(sizeof (char));
                    if (page == NULL)
                        err(1, "malloc error in answer_to_connection().");

                    page[0] = '\0';
                }

                giggle_query_result_destroy(&gqr);

                if (r.type == REQUEST_REGION) {
                    free(r.data.reg->chrm);
                    free(r.data.reg);
                    r.data.reg = NULL;
                } else {
                    free(r.data.query->file_name);
                    free(r.data.query);
                    r.data.query = NULL;
                }

                if (r.file_patterns_to_be_printed != NULL) {
                    free(r.file_patterns_to_be_printed);
                    r.file_patterns_to_be_printed = NULL;
                }
                r.full = 0;


                int ret = send_page(connection,
                                    page,
                                    MHD_HTTP_OK,
                                    "text/txt");
                free(page);
                return ret;
            }
        }
        return 0;
    }
    //}}}

    //{{{ POST 
    if (0 == strcmp (method, "POST")) {
        struct connection_info_struct *con_info = *con_cls;

	fprintf(stderr, "upload_data_size: %zu\n", *upload_data_size);

        if (0 != *upload_data_size) {
            MHD_post_process(con_info->postprocessor,
                             upload_data,
                            *upload_data_size);
            *upload_data_size = 0;

            return MHD_YES;
        } else {
		
            fprintf(stderr, "1\n");
            if (NULL != con_info->fp) {
                fclose (con_info->fp);
		con_info->fp = NULL;
	    }
          /* Now it is safe to open and inspect the file before calling
           * send_page with a response */

            int chrm_len = 50;
            char *chrm = (char *)malloc(chrm_len*sizeof(char));
            if (chrm == NULL)
                err(1, "malloc error in answer_to_connection().");

            uint32_t start, end;
            long offset;
            uint32_t num_intervals = 0;
            double mean_interval_size = 0.0;
            double genome_size =  3095677412.0;

            //uint32_t num_file_patterns = 0;
            //regex_t *regexs = NULL;
            //char **file_patterns = NULL;


            struct input_file *q_f = input_file_init(con_info->file_name);
            if (q_f == NULL) {
                fprintf(stderr, "Error with file %s.\n", con_info->file_name);
                return send_page (connection, queryfile_errorpage, MHD_HTTP_OK, "text/html");
            }

            struct giggle_query_result *gqr = NULL;

            kstring_t line = {0, 0, NULL};
            while ( q_f->input_file_get_next_interval(q_f, 
                                                     &chrm,
                                                     &chrm_len,
                                                     &start,
                                                     &end,
                                                     &offset,
                                                     &line) >= 0 ) {
                gqr = giggle_query(con_info->arg->gi, chrm, start, end, gqr);
                num_intervals += 1;
                mean_interval_size += end - start;
            }

            if (line.s != NULL)
                free(line.s);

            input_file_destroy(&q_f);

            mean_interval_size = mean_interval_size/num_intervals;

            uint32_t i,j;

            char *page = NULL, *tmp_page = NULL;;

            for(i = 0; i < gqr->num_files; i++) {
                struct file_data *fd = 
                        file_index_get(con_info->arg->gi->file_idx, i); 
                        //(struct file_data *)
                        //unordered_list_get(con_info->arg->gi->file_index, i); 

                uint32_t file_counts = giggle_get_query_len(gqr, i);
                long long n11 = (long long)(file_counts);
                long long n12 = (long long)(MAX(0,num_intervals-file_counts));
                long long n21 = (long long)
                        (MAX(0,fd->num_intervals-file_counts));

                double comp_mean = fd->mean_interval_size+mean_interval_size;

                long long n22_full = (long long)
                        MAX(n11 + n12 + n21, genome_size/comp_mean);
                long long n22 = MAX(0, n22_full - (n11 + n12 + n21));

                long double left, right, two;
                long double r = _kt_fisher_exact(n11,
                                                 n12,
                                                 n21,
                                                 n22,
                                                 &left,
                                                 &right,
                                                 &two);
                double ratio = 
                    (((double)n11/(double)n12) / ((double)n21/(double)n22));

                if (i == 0) {
                    int ret = asprintf(&tmp_page,
                                       "#%s\t"
                                       "size:%u\t"
                                       "overlaps:%u\t"
                                       "ratio:%f\t"
                                       "sig:%Lf"
                                       "combo:%Lf"
                                       "\n",
                                       fd->file_name,
                                       fd->num_intervals,
                                       file_counts,
                                       ratio,
                                       right,
                                       log2fc(ratio) * neglog10p(two));
                } else {
                    int ret = asprintf(&tmp_page,
                                       "%s"
                                       "#%s\t"
                                       "size:%u\t"
                                       "overlaps:%u\t"
                                       "ratio:%f\t"
                                       "sig:%Lf\t"
                                       "combo:%Lf"
                                       "\n",
                                       page,
                                       fd->file_name,
                                       fd->num_intervals,
                                       file_counts,
                                       ratio,
                                       right,
                                       log2fc(ratio) * neglog10p(two));
                }

                free(page);
                page = tmp_page;
            }

            int ret = send_page(connection,
                                page,
                                con_info->answercode,
                                "text/txt");
            free(page);
	    //void input_file_destroy(struct input_file **i);/
            fprintf(stderr, "2\n");

            return ret;
        }
    }
    //}}}
    
    return send_page (connection, errorpage, MHD_HTTP_BAD_REQUEST, "text/html");
}
//}}}

//{{{int server_help(int exit_code)
int server_help(int exit_code)
{
    fprintf(stderr,
            "usage: server_enrichment -i <index dir> -u <upload dir> -d <data definition> "
            "-p <port> [options]\n"
            "         options:\n"
            "             -a Set Access-Control-Allow-Origin header to pages\n" );
    return exit_code;
}
//}}}

//{{{int main(int argc, char **argv)
int main(int argc, char **argv)
{
    if (argc < 2) return server_help(EX_OK);

    int i_is_set = 0,
        u_is_set = 0,
        d_is_set = 0,
        p_is_set = 0,
        a_is_set = 0;


    char *index_dir_name = NULL,
         *upload_dir_name = NULL,
         *data_def_file_name = NULL;

    uint32_t port = 0;

    SET_ACCESS_CONTROL_HEADER = 0;

    int c;
    while((c = getopt (argc, argv, "i:u:d:p:ah")) != -1) {
        switch(c) {
            case 'i':
                i_is_set = 1;
                index_dir_name = optarg;
                break;
            case 'u':
                u_is_set = 1;
                upload_dir_name = optarg;
                break;
            case 'd':
                d_is_set = 1;
                data_def_file_name = optarg;
                break;
            case 'p':
                p_is_set = 1;
                port = atoi(optarg);
                break;
            case 'a':
                a_is_set = 1;
                SET_ACCESS_CONTROL_HEADER = 1;
                break;
            case '?':
                 if ( (optopt == 'i') ||
                      (optopt == 'u') ||
                      (optopt == 'd') ||
                      (optopt == 'p') )
                        fprintf (stderr, "Option -%c requires an argument.\n",
                                optopt);
                    else if (isprint (optopt))
                        fprintf (stderr, "Unknown option `-%c'.\n", optopt);
                    else
                    fprintf(stderr,
                            "Unknown option character `\\x%x'.\n",
                            optopt);
                return server_help(EX_USAGE);
            default:
                return server_help(EX_OK);
        }
    }

    if (i_is_set == 0) {
        fprintf(stderr, "Index directory is not set\n");
        return server_help(EX_USAGE);
    }

    if (u_is_set == 0) {
        fprintf(stderr, "Upload directory is not set\n");
        return server_help(EX_USAGE);
    }

    if (d_is_set == 0) {
        fprintf(stderr, "Data def is not set\n");
        return server_help(EX_USAGE);
    }

    if (p_is_set == 0) {
        fprintf(stderr, "Port is not set\n");
        return server_help(EX_USAGE);
    }


    FILE *fp = fopen(data_def_file_name, "r");
    if (!fp)
        err(1, "Could not open header file '%s'", data_def_file_name);
    fseek(fp, 0, SEEK_END);
    long data_def_size = ftell(fp);
    fseek(fp, 0, SEEK_SET);
    char *data_def = (char *)malloc((data_def_size + 1)*sizeof(char));
    if (data_def == NULL)
        err(1, "malloc error in main().");

    int ret = fread(data_def, data_def_size, 1, fp);
    data_def[data_def_size] = '\0';
    fclose(fp);

    struct args *arg = (struct args *)malloc(sizeof(struct args));
    if (arg == NULL)
        err(1, "malloc error in main().");

    arg->data_def = data_def;
    arg->upload_dir = upload_dir_name;

    arg->gi = giggle_load(index_dir_name,
                          block_store_giggle_set_data_handler);

    struct MHD_Daemon *daemon;

    daemon = MHD_start_daemon (MHD_USE_SELECT_INTERNALLY,
                               port,
                               NULL,
                               NULL,
                               &answer_to_connection,
                               arg,
                               MHD_OPTION_NOTIFY_COMPLETED,
                               request_completed,
                               NULL,
                               MHD_OPTION_END);
    if (NULL == daemon) return 1;

    getchar (); 
    
    free(arg->data_def);
    giggle_index_destroy(&(arg->gi));
    cache.destroy();
    MHD_stop_daemon (daemon);
    return 0;
}
//}}}
