#ifndef __UTIL_H__
#define __UTIL_H__

#define _GNU_SOURCE

#include <stdio.h>
#include <sys/stat.h>
#include <ftw.h>
#include <stdint.h>
#include <regex.h>
#include <stdint.h>

#include "giggle_index.h"
#include "file_read.h"

struct doubles_uint32_t_tuple
{
    double d1,d2;
    uint32_t u1,u2,u3;
};
int doubles_uint32_t_tuple_cmp(const void *_a, const void *_b);

struct long_uint_pair
{
    long long_val;
    uint32_t uint_val;
};
int long_uint_pair_cmp(const void *_a, const void *_b);

extern struct FTW *ftwbuf;
void check_file_read(char *file_name, FILE *fp, size_t exp, size_t obs);
int unlink_cb(const char *fpath,
              const struct stat *sb,
              int typeflag,
              struct FTW *ftwbuf);
int rmrf(char *path);
int uint32_t_cmp(const void *_a, const void *_b);
int uint64_t_cmp(const void *_a, const void *_b);
int char_p_cmp(const void *_a, const void *_b);
uint32_t bin_char_to_int(char *bin);
int long_cmp(const void *_a, const void *_b);
int parse_region(char *region_s, char **chrm, uint32_t *start, uint32_t *end);
int test_pattern_match(struct giggle_index *gi,
                       regex_t *regexs,
                       char **file_patterns,
                       uint32_t num_file_patterns,
                       uint32_t file_id,
                       uint32_t f_is_set);


double log2fc(double ratio);
long double neglog10p(long double sig);

#endif
