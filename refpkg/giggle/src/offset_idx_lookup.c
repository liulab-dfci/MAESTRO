#include <stdlib.h>
#include <stdio.h>
#include <err.h>
#include <string.h>
#include <htslib/kstring.h>

#include "giggle_index.h"
#include "offset_index.h"
#include "file_read.h"


int main(int argc, char **argv)
{

    if ((argc != 4)) {
        errx(1,
             "usage:\t%s <offest idx> <file idx> <idx>",
             argv[0]);
    }

    char *offset_idx_path = argv[1];
    char *file_idx_path = argv[2];
    uint64_t offset_id = atoll(argv[3]);

    struct offset_index *offset_idx = offset_index_load(offset_idx_path);
    struct file_index *file_idx = file_index_load(file_idx_path);

    struct file_id_offset_pair fid_off =
        offset_index_get(offset_idx, offset_id);

    struct file_data *fd = file_index_get(file_idx, fid_off.file_id);

    struct input_file *ipf = input_file_init(fd->file_name);

    ipf->input_file_seek(ipf, fid_off.offset);

    char *result = NULL;
    ipf->input_file_get_next_line(ipf, &result);

    printf("%s\n", result);

    return 0;
}
