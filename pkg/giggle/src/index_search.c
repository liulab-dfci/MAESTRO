#include <stdlib.h>
#include <stdio.h>

#include "giggle_index.h"
#include "ll.h"


int main(int argc, char **argv)
{
    ORDER = 10;
    struct simple_cache *sc = simple_cache_init(1000, 30, NULL);
    uint64_t_ll_giggle_set_data_handler();
    struct giggle_index *gi = giggle_init_index(30);
    char *path_name = argv[1];
    uint32_t r = giggle_index_directory(gi, path_name, 0);
    giggle_query_region(gi, argv[2], atoi(argv[3]), atoi(argv[4]));
}
