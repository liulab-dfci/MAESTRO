#define _GNU_SOURCE

#include <float.h>
#include <limits.h>
#include <err.h>
#include <sysexits.h>
#include <stdlib.h>
#include <stdio.h>
#include <stdint.h>
#include <string.h>
#include <math.h>

#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>

#include <ftw.h>

#include "util.h"

int long_uint_pair_cmp(const void *_a, const void *_b)
{
    struct long_uint_pair *a = (struct long_uint_pair *)_a;
    struct long_uint_pair *b = (struct long_uint_pair *)_b;

    return a->long_val - b->long_val;
}


int doubles_uint32_t_tuple_cmp(const void *_a, const void *_b)
{
    struct doubles_uint32_t_tuple *a = (struct doubles_uint32_t_tuple *)_a;
    struct doubles_uint32_t_tuple *b = (struct doubles_uint32_t_tuple *)_b;

    if (a->d1 < b->d1)
        return -1;
    else if (a->d1 > b->d1)
        return 1;
    else {
        if (a->d2 < b->d2)
            return -1;
        else if (a->d2 > b->d2)
            return 1;

        return 0;
    }
}


void check_file_read(char *file_name, FILE *fp, size_t exp, size_t obs)
{
    if (exp != obs) {
        if (feof(fp))
            errx(EX_IOERR,
                 "Error reading file \"%s\": End of file",
                 file_name);
        err(EX_IOERR, "Error reading file \"%s\"", file_name);
    }
}

int unlink_cb(const char *fpath,
              const struct stat *sb,
              int typeflag,
              struct FTW *ftwbuf)
{
    int rv = remove(fpath);
    if (rv)
        perror(fpath);
    return rv;
}

int rmrf(char *path)
{
    return nftw(path, unlink_cb, 64, FTW_DEPTH | FTW_PHYS);
}

int uint32_t_cmp(const void *_a, const void *_b)
{
    uint32_t *a = (uint32_t *)_a;
    uint32_t *b = (uint32_t *)_b;

    if (*a < *b)
        return -1;
    else if (*a > *b)
        return 1;
    else
        return 0;
}

int uint64_t_cmp(const void *_a, const void *_b)
{
    uint64_t *a = (uint64_t *)_a;
    uint64_t *b = (uint64_t *)_b;

    if (*a < *b)
        return -1;
    else if (*a > *b)
        return 1;
    else
        return 0;
}

int char_p_cmp(const void *_a, const void *_b)
{
    return strcmp(*(char **)_a, *(char **)_b);
}

uint32_t bin_char_to_int(char *bin)
{
    uint32_t i = 0;
    int j = 0;

    while (bin[j] != '\0') {
        i = i << 1;
        if (bin[j] == '1')
            i += 1;
        j+=1;
    }

    return i;
}

int long_cmp(const void *_a, const void *_b)
{
    long *a = (long *)_a;
    long *b = (long *)_b;

    if (*a < *b)
        return -1;
    else if (*a > *b)
        return 1;
    else
        return 0;
}


//{{{int parse_region(char *region_s, char **chrm, uint32_t *start, uint32_t
int parse_region(char *region_s, char **chrm, uint32_t *start, uint32_t *end)
{
    *chrm = NULL;
    *start = 0;
    *end = 0;
    uint32_t i, len = strlen(region_s);
    
    for (i = 0; i < len; ++i) {
        if (region_s[i] == ':') {
            region_s[i] = '\0';
            *chrm = strndup(region_s, strlen(region_s) + 1);
            *start = atoi(region_s + i + 1);
        } else if (region_s[i] == '-') {
            region_s[i] = '\0';
            *end = atoi(region_s + i + 1);
            if (*chrm != NULL)
                return 0;
            else 
                return 1;
        }
    }

    return 1;
}
//}}}

//{{{int test_pattern_match(struct giggle_index *gi,
int test_pattern_match(struct giggle_index *gi,
                       regex_t *regexs,
                       char **file_patterns,
                       uint32_t num_file_patterns,
                       uint32_t file_id,
                       uint32_t f_is_set) {
    if (f_is_set == 0)
        return 1;

    struct file_data *fd = file_index_get(gi->file_idx, file_id);

    int match = 0;
    uint32_t j;
    for(j = 0; j < num_file_patterns; j++) {
        int r = regexec(&(regexs[j]), fd->file_name, 0, NULL, 0);
        if (r == 0) {
            match = 1;
            break;
        } else if (r != REG_NOMATCH) {
            char msgbuf[100];
            regerror(r, &regexs[j], msgbuf, sizeof(msgbuf));
            errx(EX_USAGE,
                 "Regex '%s' match failed: %s\n",
                 file_patterns[file_id],
                 msgbuf);
        }
    }

    return match;
}
//}}}

//{{{double log2fc(double ratio)
double log2fc(double ratio)
{
    if (fabs(ratio) < 0.0001)
        return 0.0;

    if (ratio < 1) {
        ratio = 1.0/ratio;
        return -1.0 * log2(ratio);
    }

    return log2(ratio);
}
//}}}

//{{{double neglog10p(double sig)
long double neglog10p(long double sig)
{
    if (fabsl(sig) < -DBL_MAX)
        return 10.0;
    return -1.0 * log10l(sig);
}
//}}}

