BIN=../bin
OBJ=../obj
LIBD=../lib
HTS_ROOT=$(LIBD)/htslib
LIBMHD_INCLUDES=$(HOME)/usr/local/include
LIBMHD_LIBS=$(HOME)/usr/local/lib
CFLAGS=-O2 -D_FILE_OFFSET_BITS=64 -Werror -Wuninitialized
#CFLAGS=-g -D_FILE_OFFSET_BITS=64 -Werror -Wuninitialized
CFLAGS+=-DBLOCK_STORE 
INCLUDES=-I$(HTS_ROOT)
LIBS=-ldl -lz -lm -pthread -lcurl
LIBS+=-lcrypto
CC=gcc


LIB=bpt.o \
    disk_store.o \
    cache.o \
    data_reg.o \
    file_read.o \
    giggle_index.o \
    lists.o \
    ll.o \
    util.o \
    timer.o \
    wah.o \
    kfunc.o \
    index.o \
    offset_index.o \
    search.o \
    leaf.o \
    pq.o \
    jsw_avltree.o \
    fastlz.o


.SUFFIXES:

OBJS=$(addprefix $(OBJ)/, $(LIB))

.SECONDARY: $(OBJS)


PROG=giggle \
     api_test \
     offset_idx_lookup

LIST=$(addprefix $(BIN)/, $(PROG))


all: check-env $(LIST) library

server: INCLUDES+=-I$(LIBMHD_INCLUDES) -L$(LIBMHD_LIBS)
server: LIBS+=-lmicrohttpd 
server: $(BIN)/server_enrichment

library: $(OBJS)
	ar -cvq $(LIBD)/libgiggle.a $(OBJS)

$(OBJ)/%.o: %.c 
	$(CC) -c $(CFLAGS) -o $@ $< \
		-I$(HTS_ROOT)


$(BIN)/%: %.c $(OBJS)
	$(CC) $(CFLAGS) -o $@ $^ \
		$(INCLUDES) \
		-DSAMTOOLS=1 \
		$(HTS_ROOT)/libhts.a \
		$(LIBS)


check-env:
ifdef TRAVIS_BUILD_DIR
HTS_ROOT=$(TRAVIS_BUILD_DIR)/htslib
endif


