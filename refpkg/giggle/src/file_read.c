#define _GNU_SOURCE

#include <htslib/bgzf.h>
#include <htslib/kstring.h>
#include <htslib/hts.h>
#include <stdlib.h>
#include <stdio.h>
#include <stdint.h>
#include <string.h>
#include <err.h>
#include <sysexits.h>

#include "util.h"
#include "file_read.h"

//{{{int scan_s(char *str, int str_len, int *s, int *e, const char delim)
int scan_s(char *str, int str_len, int *s, int *e, const char delim) {
    if (*e == str_len)
        return -1;
    for (*e = *s; *e <= str_len; *e+=1) {
        if ((str[*e] == delim) || (*e == str_len)) {
            return *e - *s;
        }
    }
    return -2;
}
//}}}

//{{{struct input_file *input_file_init(char *file_name)
struct input_file *input_file_init(char *file_name)
{
    struct input_file *i = (struct input_file *)
            malloc(sizeof(struct input_file));
    if (i == NULL)
        err(1, "malloc error in input_file_init().");

    i->file_name = strdup(file_name);

    // .bed.gz
    if ( (strlen(i->file_name) > 7) &&
            strcmp(".bed.gz", file_name + strlen(i->file_name) - 7) == 0) {
        i->type = BED;
    // .vcf.gz
    } else if ( (strlen(i->file_name) > 7) &&
                strcmp(".vcf.gz", file_name + strlen(i->file_name) - 7) == 0) {
        i->type = VCF;
    } else if ( (strlen(i->file_name) > 4) &&
                strcmp(".bcf", file_name + strlen(i->file_name) - 4) == 0) {
        i->type = BCF;
    } else {
        fprintf(stderr, "File type not supported '%s'.\n", i->file_name);
        free(i->file_name);
        free(i);
        return NULL;
    }

    if (i->type == BED) {
        i->kstr = (kstring_t*)calloc(1, sizeof(kstring_t));
        if (i->kstr == NULL) {
            fprintf(stderr, "calloc error in input_file_init().");
            return NULL;
        }

        if ((i->bed_fp = bgzf_open(file_name, "r")) == 0) {
            fprintf(stderr ,"Could not open file '%s'\n", file_name);
            return NULL;
        }

        if ( bgzf_compression(i->bed_fp) != bgzf ) {
            fprintf(stderr, "Not a BGZF file '%s'\n", file_name);
            return NULL;
        }

        if ( !i->bed_fp->is_compressed ) {
            fprintf(stderr, "Not a bgzip compressed file '%s'\n", file_name);
            return NULL;
        }

        i->last_offset = 0;

        i->input_file_get_next_interval = 
            input_file_get_next_interval_bed;
        i->input_file_get_next_line = 
            input_file_get_next_line_bgzf;
        i->input_file_seek = 
            input_file_seek_bgzf;
    } else if ( (i->type == VCF) || (i->type == BCF) ) {
        i->bcf_fp = hts_open(i->file_name,"rb");
        if (!i->bcf_fp) {
            fprintf(stderr, "Could not read file '%s'", i->file_name);
            return NULL;
        }

        if (i->bcf_fp->format.compression != bgzf) {
            fprintf(stderr, "Not a BGZF file '%s'\n", i->file_name);
            return NULL;
        }

        i->line = bcf_init1();
        i->hdr = bcf_hdr_read(i->bcf_fp);

        if ( !i->hdr ) {
            fprintf(stderr, "Could not read the header '%s'", i->file_name);
            return NULL;
        }

        htsFormat type = *hts_get_format(i->bcf_fp);
        if (type.format == bcf) {
            fprintf(stderr, "File type not supported '%s'.\n", i->file_name);
            free(i->file_name);
            free(i);
            return NULL;
        } else {
            hts_close(i->bcf_fp);

            i->kstr = (kstring_t*)calloc(1, sizeof(kstring_t));

            if (i->kstr == NULL) {
                fprintf(stderr, "calloc error in input_file_init().");
                return NULL;
            }

            if ((i->bed_fp = bgzf_open(file_name, "r")) == 0) {
                fprintf(stderr, "Could not open file '%s'\n", i->file_name);
                return NULL;
            }

            if ( !i->bed_fp->is_compressed ) {
                fprintf(stderr,
                        "Not a bgzip compressed file '%s'\n",
                        i->file_name);
                return NULL;
            }

            // move past the header
            while (bgzf_getline(i->bed_fp, '\n', i->kstr) >= 0) {
                if (i->kstr->s[0] == '#')
                    i->last_offset = bgzf_tell(i->bed_fp);
                else
                    break;
            }

            if (bgzf_seek(i->bed_fp, i->last_offset, SEEK_SET) != 0) {
                fprintf(stderr, "Error moving past header '%s'", i->file_name);
                return NULL;
            }

            i->input_file_get_next_interval = 
                input_file_get_next_interval_vcf;
            i->input_file_get_next_line = 
                input_file_get_next_line_bgzf;
            i->input_file_seek = 
                input_file_seek_bgzf;
        }
    }

    return i;
}
//}}}

//{{{void input_file_destroy(struct input_file **i)
void input_file_destroy(struct input_file **i)
{
    if ((*i)->type == VCF) {
        bcf_hdr_destroy((*i)->hdr); 
    }
    free((*i)->file_name);
    bgzf_close((*i)->bed_fp);
    free((*i)->kstr->s);
    free((*i)->kstr);
    free(*i);
    *i = NULL;
}
//}}}

//{{{int input_file_get_next_interval_bed(struct input_file *i,
int input_file_get_next_interval_bed(struct input_file *i,
                                     char **chrm,
                                     int *chrm_len,
                                     uint32_t *start,
                                     uint32_t *end,
                                     long *offset,
                                     kstring_t *line)
{
    *offset = i->last_offset;
    int ret = bgzf_getline(i->bed_fp, '\n', i->kstr);
    i->last_offset = bgzf_tell(i->bed_fp);

    if (ret < 0 )
        return ret;

    line->l = 0;
    kputs(i->kstr->s, line);
 
    int s=0, e=0;
    int col = 0;
    while ( scan_s(i->kstr->s, i->kstr->l, &s, &e, '\t') >= 0 ) {
        col += 1;
        if (col == 1) {
            if ((e - s + 1) > *chrm_len) {
                *chrm_len = *chrm_len *2;
                *chrm = realloc(*chrm, *chrm_len * sizeof(char));
                if (*chrm == NULL)
                    err(1,
                        "realloc error in input_file_get_next_interval_bed().");
            }

            if( strncmp("chr", i->kstr->s + s, 3) == 0) {
                memcpy(*chrm, i->kstr->s + s + 3, e - s - 3);
                (*chrm)[e - s - 3] = '\0';
            } else {
                memcpy(*chrm, i->kstr->s + s, e - s);
                (*chrm)[e - s] = '\0';
            }

        } else if (col == 2) {
            *start = strtol(i->kstr->s + s, NULL, 0);
        } else if (col == 3) {
            *end = strtol(i->kstr->s + s, NULL, 0);
        } else {
            break;
        }
        s = e + 1;
    }

    return ret;
}
//}}}

//{{{int input_file_get_next_interval_vcf(struct input_file *i,
int input_file_get_next_interval_vcf(struct input_file *i,
                                     char **chrm,
                                     int *chrm_len,
                                     uint32_t *start,
                                     uint32_t *end,
                                     long *offset,
                                     kstring_t *line)
{
    *offset = i->last_offset;
    int ret = bgzf_getline(i->bed_fp, '\n', i->kstr);
    i->last_offset = bgzf_tell(i->bed_fp);

    if (ret < 0 )
        return ret;

    vcf_parse(i->kstr, i->hdr, i->line);

    line->l = 0;
    kputs(i->kstr->s, line);

    int lensize = sizeof(uint32_t);
    uint32_t *len = (uint32_t *) calloc(1,lensize);
    if (len == NULL)
        err(1, "calloc error in input_file_get_next_interval_vcf()");

    int size = sizeof(uint32_t);
    if (bcf_get_info_int32(i->hdr, i->line, "END", &end, &size) < 0) {
        size = sizeof(uint32_t);
        if (bcf_get_info_int32(i->hdr, i->line, "SVLEN", &len, &lensize) >= 0)
            *end = i->line->pos + *len;
        else
            *end = i->line->pos + i->line->rlen;
    }

    *end = *end + 1;

    const char *_chrm = bcf_hdr_id2name(i->hdr, i->line->rid);

    while (strlen(_chrm) + 1 > *chrm_len) {
        *chrm_len = *chrm_len *2;
        *chrm = realloc(*chrm, *chrm_len * sizeof(char));
        if (*chrm == NULL)
            err(1, "realloc error in input_file_get_next_interval_vcf()");
    }

    if( strncmp("chr", _chrm, 3) == 0)
        memcpy(*chrm, _chrm + 3, (strlen(_chrm) + 1) - 3);
    else
        memcpy(*chrm, _chrm, strlen(_chrm) + 1);

    *start = i->line->pos + 1;

    return ret;
}
//}}}

//{{{int input_file_get_next_line_bgzf(struct input_file *i,
int input_file_get_next_line_bgzf(struct input_file *i,
                                  char **str)
{
    int ret = bgzf_getline(i->bed_fp, '\n', i->kstr);
    i->last_offset = bgzf_tell(i->bed_fp);
    *str = i->kstr->s;
    return ret;
}
//}}}

//{{{void input_file_get_curr_line_bgzf(struct input_file *i,
void input_file_get_curr_line_bgzf(struct input_file *i,
                                  char **str)
{
    //fprintf(stdout, "****%zu\n", i->kstr->l);

    if (i->type == VCF) {
        i->kstr->l = 0;
        vcf_format1(i->hdr, i->line, i->kstr);
    }
    *str = i->kstr->s;
}
//}}}

//{{{void input_file_seek_bgzf(struct input_file *i, long offset)
void input_file_seek_bgzf(struct input_file *i, long offset)
{
    if (bgzf_seek(i->bed_fp, offset, SEEK_SET) != 0)
        err(EX_DATAERR, "Error seeking in bgzf file.");
    i->last_offset = offset;
}
//}}}

//{{{void file_data_store(void *v, FILE *f, char *file_name)
void file_data_store(void *v, FILE *f, char *file_name)
{
    struct file_data *fd = (struct file_data *)v;
    uint32_t size = strlen(fd->file_name) +
                    sizeof(uint32_t) + 
                    sizeof(double);

    if (fwrite(&size, sizeof(uint32_t), 1, f) != 1)
        err(EX_IOERR,
            "Error writing file_data size to '%s'.", file_name);

    if (fwrite(fd->file_name,
               sizeof(char),
               strlen(fd->file_name), f) != strlen(fd->file_name))
        err(EX_IOERR,
            "Error writing file_data file_name to '%s'.", file_name);

    if (fwrite(&(fd->num_intervals),
               sizeof(uint32_t),
               1, f) != 1)
        err(EX_IOERR,
            "Error writing file_data num_intervals to '%s'.", file_name);

    if (fwrite(&(fd->mean_interval_size),
               sizeof(double),
               1, f) != 1)
        err(EX_IOERR,
            "Error writing file_data mean_interval_size to '%s'.", file_name);

}
//}}}

//{{{ void *file_data_load(FILE *f, char *file_name)
void *file_data_load(FILE *f, char *file_name)
{
    uint32_t size;
    size_t fr = fread(&size, sizeof(uint32_t), 1, f);
    check_file_read(file_name, f, 1, fr);

    uint32_t str_len = size - sizeof(uint32_t) - sizeof(double);

    struct file_data *fd = (struct file_data *)
        calloc(1, sizeof(struct file_data));
    if (fd == NULL)
        err(1, "calloc error in file_data_load()");

    fd->file_name = (char *)calloc(str_len + 1, sizeof(char));
    if (fd->file_name == NULL)
        err(1, "calloc error in file_data_load()");

    fr = fread(fd->file_name, sizeof(char), str_len, f);
    check_file_read(file_name, f, str_len, fr);

    fd->file_name[str_len] = '\0';

    uint32_t v;
    fr = fread(&v, sizeof(uint32_t), 1, f);
    check_file_read(file_name, f, 1, fr);

    fd->num_intervals = v;

    fr = fread(&(fd->mean_interval_size), sizeof(double), 1, f);
    check_file_read(file_name, f, 1, fr);

    return (void *)fd;
}
//}}}

//{{{void file_data_free(void **v)
void file_data_free(void **v)
{
    struct file_data **fd = (struct file_data **)v;
    free((*fd)->file_name);
    free(*fd);
    *fd = NULL;
}
//}}}
