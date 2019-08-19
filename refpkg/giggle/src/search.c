#define _GNU_SOURCE
#include <stdlib.h>
#include <stdio.h>
#include <err.h>
#include <string.h>
#include <getopt.h>
#include <ctype.h>
#include <sysexits.h>
#include <regex.h>
#include <htslib/kstring.h>
#include <math.h>

#include "giggle_index.h"
#include "wah.h"
#include "cache.h"
#include "file_read.h"
#include "util.h"
#include "kfunc.h"
#include "ll.h"

#define MAX(X, Y) (((X) > (Y)) ? (X) : (Y))

int search_help(int exit_code);
int print_giggle_query_result(struct giggle_query_result *gqr,
                              struct giggle_index *gi,
                              regex_t *regexs,
                              char **file_patterns,
                              uint32_t num_file_patterns,
                              uint32_t num_intervals,
                              double mean_interval_size,
                              long long genome_size,
                              uint32_t f_is_set,
                              uint32_t v_is_set,
                              uint32_t c_is_set,
                              uint32_t s_is_set,
                              uint32_t o_is_set);

//{{{ int search_help(int exit_code)
int search_help(int exit_code)
{
    fprintf(stderr,
"%s, v%s\n"
"usage:   %s search -i <index directory> [options]\n"
"         options:\n"
"             -i giggle index directory\n"
"             -r <regions (CSV)>\n"
"             -q <query file>\n"
"             -o give reuslts per record in the query file (omits empty results)\n"
"             -c give counts by indexed file\n"
"             -s give significance by indexed file (requires query file)\n"
"             -v give full record results\n"
"             -f print results for files that match a pattern (regex CSV)\n"
"             -g genome size for significance testing (default 3095677412)\n"
"             -l list the files in the index\n",
            PROGRAM_NAME, VERSION, PROGRAM_NAME);
    return exit_code;
}
//}}}

//{{{int print_giggle_query_result(struct giggle_query_result *gqr,
int print_giggle_query_result(struct giggle_query_result *gqr,
                              struct giggle_index *gi,
                              regex_t *regexs,
                              char **file_patterns,
                              uint32_t num_file_patterns,
                              uint32_t num_intervals,
                              double mean_interval_size,
                              long long genome_size,
                              uint32_t f_is_set,
                              uint32_t v_is_set,
                              uint32_t c_is_set,
                              uint32_t s_is_set,
                              uint32_t o_is_set)
{
    if (gqr == NULL)
        return EX_OK;

    uint32_t i,j;

    if (s_is_set == 1) {
        printf("#file\t"
               "file_size\t"
               "overlaps\t"
               "odds_ratio\t"
               "fishers_two_tail\t"
               "fishers_left_tail\t"
               "fishers_right_tail\t"
               "combo_score\n");
    }

    for(i = 0; i < gqr->num_files; i++) {
        struct file_data *fd = file_index_get(gi->file_idx, i);
        if (test_pattern_match(gi,
                               regexs,
                               file_patterns,
                               num_file_patterns,
                               i,
                               f_is_set)) {
            if ( (v_is_set == 1) && (giggle_get_query_len(gqr, i) > 0 )){
                //printf("#%s\n", fd->file_name);
                char *result;
                struct giggle_query_iter *gqi =
                    giggle_get_query_itr(gqr, i);
                while (giggle_query_next(gqi, &result) == 0)
                    printf("%s\t%s\n", result, fd->file_name);
                giggle_iter_destroy(&gqi);
            } else if (c_is_set == 1) {
                printf("#%s\t"
                       "size:%u\t"
                       "overlaps:%u\n",
                       fd->file_name,
                       fd->num_intervals,
                       giggle_get_query_len(gqr, i));
            } else if (s_is_set == 1) {
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

                /*
                double ratio = 
                        (((double)n11/(double)MAX(1,n12)) / 
                         ((double)n21/(double)n22));
                */

                double ratio = 
                        (((double)(n11 + 1)/(double)MAX(1,n12)) / 
                         ((double)(n21 + 1)/(double)(n22 + 1)));

                printf("%s\t"
                       "%u\t"
                       "%u\t"
                       "%.17g\t"
                       "%.17Lg\t"
                       "%.17Lg\t"
                       "%.17Lg\t"
                       "%.17Lg\t"
                       "\n",
                       fd->file_name,
                       fd->num_intervals,
                       file_counts,
                       ratio,
                       two,
                       left,
                       right,
                       log2fc(ratio) * neglog10p(two));
                /*
                printf("#%s\t"
                       "size:%u\t"
                       "overlaps:%u\t"
                       "ratio:%f\t"
                       "combo:%.17Lg\t"
                       "right:%.17Lg\t"
                       "left:%.17Lg\t"
                       "two:%.17Lg\t"
                       "n11:%lld\t"
                       "n12:%lld\t"
                       "n21:%lld\t"
                       "n22:%lld\t"
                       "fd->mean_interval_size:%f\t"
                       "mean_interval_size:%f\t"
                       "comp_mean:%f"
                       "\n",
                       fd->file_name,
                       fd->num_intervals,
                       file_counts,
                       ratio,
                       log2fc(ratio) * neglog10p(two),
                       right,
                       left,
                       two,
                       n11,
                       n12,
                       n21,
                       n22,
                       fd->mean_interval_size,
                       mean_interval_size,
                       comp_mean);
                       */

            }
        }
    }

    return EX_OK;
}
//}}}

int search_main(int argc, char **argv, char *full_cmd)
{
    if (argc < 2) return search_help(EX_OK);

    uint32_t num_chrms = 100;
    int c;
    char *index_dir_name = NULL,
         *regions = NULL,
         *query_file_name = NULL,
         *file_patterns_to_be_printed = NULL;


    char *i_type = "i";

    int i_is_set = 0,
        l_is_set = 0,
        r_is_set = 0,
        q_is_set = 0,
        c_is_set = 0,
        s_is_set = 0,
        v_is_set = 0,
        f_is_set = 0,
        o_is_set = 0;

    double genome_size =  3095677412.0;

    //{{{ cmd line param parsing
    //{{{ while((c = getopt (argc, argv, "i:r:q:cvf:h")) != -1) {
    while((c = getopt (argc, argv, "i:r:q:csvof:g:lh")) != -1) {
        switch (c) {
            case 'i':
                i_is_set = 1;
                index_dir_name = optarg;
                break;
            case 'r':
                r_is_set = 1;
                regions = optarg;
                break;
            case 'q':
                q_is_set = 1;
                query_file_name = optarg;
                break;
            case 'c':
                c_is_set = 1;
                break;
            case 's':
                s_is_set = 1;
                break;
            case 'v':
                v_is_set = 1;
                break;
            case 'o':
                o_is_set = 1;
                break;
            case 'f':
                f_is_set = 1;
                file_patterns_to_be_printed = optarg;
                break;
            case 'g':
                genome_size =  atof(optarg);
                break;
            case 'l':
                l_is_set = 1;
                break;
            case 'h':
                return search_help(EX_OK);
            case '?':
                 if ( (optopt == 'i') ||
                      (optopt == 'r') ||
                      (optopt == 'q') ||
                      (optopt == 'f') )
                        fprintf (stderr, "Option -%c requires an argument.\n",
                                optopt);
                    else if (isprint (optopt))
                        fprintf (stderr, "Unknown option `-%c'.\n", optopt);
                    else
                    fprintf(stderr,
                            "Unknown option character `\\x%x'.\n",
                            optopt);
                return search_help(EX_USAGE);
            default:
                return search_help(EX_OK);
        }
    }
    //}}}
    
    if (i_is_set == 0) {
        fprintf(stderr, "Index directory is not set\n");
        return search_help(EX_USAGE);
    } 

    if (l_is_set == 1) {
        char **names = NULL;
        uint32_t *num_intervals = NULL;
        double *mean_interval_sizes = NULL;
        uint32_t num_files = giggle_get_indexed_files(index_dir_name,
                                                      &names,
                                                      &num_intervals,
                                                      &mean_interval_sizes);
        uint32_t i;
        printf("File name\tNumber of intervals\tMean interval size\n");
        for (i = 0; i < num_files; ++i) {
            printf("%s\t%u\t%lf\n",
                   names[i],
                   num_intervals[i],
                   mean_interval_sizes[i]);
            free(names[i]);
        }
        free(names);
        free(num_intervals);
        free(mean_interval_sizes);
        return EX_OK;
    } 

    if ( (s_is_set == 1) && (q_is_set ==0)) {
        fprintf(stderr, "Significance testing requires a query file input\n");
        return search_help(EX_USAGE);
    }

    // need either a regions or a query file, but not both
    if ((r_is_set == 0) && (q_is_set == 0)) {
        fprintf(stderr, "Neither regions nor query file is set\n");
        return search_help(EX_USAGE);
    } if ((r_is_set == 1) && (q_is_set == 1)) {
        fprintf(stderr, "Both regions and query file is set\n");
        return search_help(EX_USAGE);
    }

    if ((v_is_set == 0) && (s_is_set == 0))
        c_is_set = 1;
    //}}}

    uint32_t num_file_patterns = 0;
    regex_t *regexs = NULL;
    char **file_patterns = NULL;

    //{{{ comiple file name regexs
    if (f_is_set == 1) {
        int s = 0, e = 0;
        while (scan_s(file_patterns_to_be_printed,
                      strlen(file_patterns_to_be_printed),
                      &s,
                      &e,
                      ',') >= 0) {
            num_file_patterns += 1;
            s = e + 1;
        }

        if (num_file_patterns == 0) {
            fprintf(stderr, "No file patterns detected.\n");
            return search_help(EX_USAGE);
        }

        regexs = (regex_t *)
                malloc(num_file_patterns * sizeof(regex_t));
        if (regexs == NULL)
            err(1, "malloc error  in search_main().");

        file_patterns = (char **)
                malloc(num_file_patterns * sizeof(char *));
        if (file_patterns == NULL)
            err(1, "malloc error  in search_main().");

        uint32_t i = 0;
        s = 0;
        e = 0;
        while (scan_s(file_patterns_to_be_printed,
                      strlen(file_patterns_to_be_printed),
                      &s,
                      &e,
                      ',') >= 0) {
            file_patterns[i] = strndup(file_patterns_to_be_printed + s, e-s);
            int r = regcomp(&(regexs[i]), file_patterns[i], 0);
            if (r != 0) {
                errx(EX_USAGE,
                     "Could not compile regex '%s'",
                     file_patterns[i]);
            }
            i += 1;
            s = e + 1;
        }
    }
    //}}}

    struct giggle_index *gi =
                giggle_load(index_dir_name,
                            block_store_giggle_set_data_handler);

    if (gi == NULL)
        errx(1, "Error loading giggle index %s.", index_dir_name);

    struct giggle_query_result *gqr = NULL;

    uint32_t num_intervals = 0;
    double mean_interval_size = 0.0;

    if (r_is_set == 1) {
        // search the list of regions
        uint32_t i, last = 0, len = strlen(regions);
        char *chrm;
        uint32_t start, end;
        int ret;
        for (i = 0; i <= len; ++i) {
            if ((regions[i] == ',') || (regions[i] == '\0') ) {
                regions[i] = '\0';
                char *region;
                ret = asprintf(&region, "%s", regions + last);
                if (parse_region(region, &chrm, &start, &end) == 0) {
                    gqr = giggle_query(gi, chrm, start, end, gqr);
                    free(region);
                } else {
                    errx(EX_USAGE,
                         "Error parsing region '%s'\n",
                         regions + last);
                }
                last = i + 1;
            }
        }
    } else if (q_is_set == 1) {
        // search a file
        int chrm_len = 50;
        char *chrm = (char *)malloc(chrm_len*sizeof(char));
        if (chrm == NULL)
            err(1, "malloc error  in search_main().");
        uint32_t start, end;
        long offset;
        kstring_t line = {0, 0, NULL};

        struct input_file *q_f = input_file_init(query_file_name);
        if (q_f == NULL)
            errx(1, "Error loading query file %s.", query_file_name);

        while ( q_f->input_file_get_next_interval(q_f, 
                                                  &chrm,
                                                  &chrm_len,
                                                  &start,
                                                  &end,
                                                  &offset,
                                                  &line) >= 0 ) {
            gqr = giggle_query(gi, chrm, start, end, gqr);
            if ( (o_is_set == 1) && (gqr->num_hits > 0) ) {
                char *str;
                input_file_get_curr_line_bgzf(q_f, &str);
                printf("##%s",str);
                // Ugh
                if (q_f->type == BED)
                    printf("\n");
                int r = print_giggle_query_result(gqr,
                                                  gi,
                                                  regexs,
                                                  file_patterns,
                                                  num_file_patterns,
                                                  num_intervals,
                                                  mean_interval_size,
                                                  genome_size,
                                                  f_is_set,
                                                  v_is_set,
                                                  c_is_set,
                                                  s_is_set,
                                                  o_is_set);
                giggle_query_result_destroy(&gqr);
            }
            num_intervals += 1;
            mean_interval_size += end - start;
        }

        free(chrm);
        if (line.s != NULL)
            free(line.s);

        mean_interval_size = mean_interval_size/num_intervals;
    }


    int r = print_giggle_query_result(gqr,
                                      gi,
                                      regexs,
                                      file_patterns,
                                      num_file_patterns,
                                      num_intervals,
                                      mean_interval_size,
                                      genome_size,
                                      f_is_set,
                                      v_is_set,
                                      c_is_set,
                                      s_is_set,
                                      o_is_set);

    giggle_query_result_destroy(&gqr);
    giggle_index_destroy(&gi);
    cache.destroy();
    free(full_cmd);
    return r;
}
