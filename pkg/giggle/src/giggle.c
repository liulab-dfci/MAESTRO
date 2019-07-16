#define _GNU_SOURCE
#include <stdlib.h>
#include <stdio.h>
#include <err.h>
#include <string.h>
#include <sysexits.h>
#include "giggle_index.h"
#include "wah.h"
#include "cache.h"

int giggle_help(int argc, char **argv, int exit_code);
int index_main(int argc, char **argv, char *full_cmd);
int search_main(int argc, char **argv, char *full_cmd);

int main(int argc, char **argv)
{
    if (argc < 2) return giggle_help(argc, argv, 0);
    char *cmd = argv[1];
    char *full_cmd = strdup(argv[0]);
    int i, r, quote_next = 0;
    for (i = 1; i < argc; ++i) {
        if (quote_next == 1) {
            r = asprintf(&full_cmd, "%s \"%s\"", full_cmd, argv[i]);
            if (r == -1) err(EX_OSERR, "asprintf error");
            quote_next = 0;
        } else {
            r = asprintf(&full_cmd, "%s %s", full_cmd, argv[i]);
            if (r == -1) err(EX_OSERR, "asprintf error");
        }

        if ( (strcmp(argv[i], "-p") == 0) || (strcmp(argv[i], "-g") == 0))
            quote_next = 1;
    }

    if (strcmp(cmd,"index") == 0)
        return index_main(argc-1, argv+1, full_cmd);
    else if (strcmp(cmd,"search") == 0)
        return search_main(argc-1, argv+1, full_cmd);
    else {
        fprintf(stderr, "Unknown command\n");
        return giggle_help(argc, argv, EX_USAGE);
    }
    free(full_cmd);
}

int giggle_help(int argc, char **argv, int exit_code)
{
    fprintf(stderr,
"%s, v%s\n"
"usage:   %s <command> [options]\n"
"         index     Create an index\n"
"         search    Search an index\n",
            PROGRAM_NAME, VERSION, PROGRAM_NAME);
    return exit_code;
}
