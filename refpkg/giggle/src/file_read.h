#ifndef __FILE_READ_H__
#define __FILE_READ_H__

#include <htslib/bgzf.h>
#include <htslib/vcf.h>
#include <htslib/kstring.h>
#include <stdint.h>

struct file_data
{
    char *file_name;
    uint32_t num_intervals;
    double mean_interval_size;
};

void file_data_free(void **v);
void file_data_store(void *v, FILE *f, char *file_name);
void *file_data_load(FILE *f, char *file_name);

struct input_file
{
    char *file_name;
    bcf1_t *line;
    BGZF *bed_fp;
    htsFile *bcf_fp;
    bcf_hdr_t *hdr;
    kstring_t *kstr;
    long last_offset;
    enum {BED,VCF,BCF} type;
    int (*input_file_get_next_interval)(struct input_file *i,
                                        char **chrm,
                                        int *chrm_len,
                                        uint32_t *start,
                                        uint32_t *end,
                                        long *offset,
                                        kstring_t *line);
    int (*input_file_get_next_line)(struct input_file *i,
                                    char **str);
    void (*input_file_seek)(struct input_file *i, long offset);
};

int scan_s(char *str, int str_len, int *s, int *e, const char delim);
struct input_file *input_file_init(char *file_name);
void input_file_destroy(struct input_file **i);
void input_file_seek_bgzf(struct input_file *i, long offset);
int input_file_get_next_interval_bed(struct input_file *i,
                                     char **chrm,
                                     int *chrm_len,
                                     uint32_t *start,
                                     uint32_t *end,
                                     long *offset,
                                     kstring_t *line);
int input_file_get_next_interval_vcf(struct input_file *i,
                                     char **chrm,
                                     int *chrm_len,
                                     uint32_t *start,
                                     uint32_t *end,
                                     long *offset,
                                     kstring_t *line);
int input_file_get_next_line_bgzf(struct input_file *i,
                                  char **str);

void input_file_get_curr_line_bgzf(struct input_file *i,
                                   char **str);
#endif
