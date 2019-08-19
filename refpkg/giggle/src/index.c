#include <stdlib.h>
#include <stdio.h>
#include <err.h>
#include <sysexits.h>
#include <getopt.h>
#include <ctype.h>
#include <inttypes.h>
#include "giggle_index.h"
#include "wah.h"
#include "cache.h"
#include "ll.h"

int index_help(int exit_code);
int index_main(int argc, char **argv, char *full_cmd);


int index_help(int exit_code)
{
    fprintf(stderr,
"%s, v%s\n"
"usage:   %s index -i <input files> -o <output dir> -f\n"
"         options:\n"
"             -s  Files are sorted\n"
"             -i  Files to index (e.g. data/*.gz)\n"
"             -o  Index output directory\n"
"             -f  For reindex if output directory exists\n",
            PROGRAM_NAME, VERSION, PROGRAM_NAME);
    return exit_code;
}

int index_main(int argc, char **argv, char *full_cmd)
{
    if (argc < 2) return index_help(EX_OK);

    uint32_t num_chrms = 100;
    int c;
    char *input_dir_name = NULL,
         *output_dir_name = NULL;
    char *i_type = "i";

    int i_is_set = 0,
        o_is_set = 0,
        s_is_set = 0,
        f_is_set = 0;

    while((c = getopt (argc, argv, "i:o:fsh")) != -1) {
        switch (c) {
            case 'i':
                i_is_set = 1;
                input_dir_name = optarg;
                break;
            case 'o':
                o_is_set = 1;
                output_dir_name = optarg;
                break;
            case 'f':
                f_is_set = 1;
                break;
            case 's':
                s_is_set = 1;
                break;
            case 'h':
                return index_help(EX_OK);
            case '?':
                 if ( (optopt == 'i') ||
                      (optopt == 'o') )
                        fprintf (stderr, "Option -%c requires an argument.\n",
                                optopt);
                    else if (isprint (optopt))
                        fprintf (stderr, "Unknown option `-%c'.\n", optopt);
                    else
                    fprintf(stderr,
                            "Unknown option character `\\x%x'.\n",
                            optopt);
                return index_help(EX_USAGE);
            default:
                return index_help(EX_OK);
        }
    }
    
    if (i_is_set == 0) {
        fprintf(stderr, "Input file is not set\n");
        return index_help(EX_USAGE);
    } else if (o_is_set == 0) {
        fprintf(stderr, "Output file is not set\n");
        return index_help(EX_USAGE);
    }

    struct giggle_index *gi;
        
    if (s_is_set == 1) {
        uint64_t num_intervals = giggle_bulk_insert(input_dir_name,
                                                    output_dir_name,
                                                    f_is_set);
        fprintf(stderr, "Indexed %" PRIu64 " intervals.\n", num_intervals);
    } else {
        gi = giggle_init(num_chrms,
                         output_dir_name,
                         f_is_set,
                         uint64_t_ll_giggle_set_data_handler);
        if (gi == NULL)
            return EX_DATAERR;

        uint32_t r = giggle_index_directory(gi, input_dir_name, 0);
        
        fprintf(stderr, "Indexed %u intervals.\n", r);

#if BLOCK_STORE
        giggle_data_handler.write_tree = giggle_write_tree_leaf_data;
#endif

        r = giggle_store(gi);

        if (r != 0)
            errx(1, "Error storing giggle index.");
        
        giggle_index_destroy(&gi);
        cache.destroy();
    }
    return EX_OK;
}
