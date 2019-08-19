#define _GNU_SOURCE

#include "util.h"

#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <inttypes.h>
#include <string.h>
#include <err.h>
#include <sysexits.h>
#include <glob.h>
#include <string.h>

#include "unity.h"
#include "bpt.h"
#include "giggle_index.h"
#include "lists.h"
#include "file_read.h"
#include "wah.h"
#include "ll.h"
#include "jsw_avltree.h"
#include "pq.h"

#define MIN(X, Y) (((X) < (Y)) ? (X) : (Y))
#define MAX(X, Y) (((X) > (Y)) ? (X) : (Y))


void setUp(void) { }
void tearDown(void) { }

//{{{ void test_giggle_bulk_insert(void)
void test_generic_offset_index_basic(void)
{
}
