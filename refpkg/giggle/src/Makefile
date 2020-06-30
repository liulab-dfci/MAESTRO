BIN=../bin
OBJ=../obj
LIBD=../lib
HTS_ROOT=$(LIBD)/htslib
CFLAGS+=-DBLOCK_STORE  -O2 -D_FILE_OFFSET_BITS=64 -Werror -Wuninitialized
INCLUDES=-I$(HTS_ROOT)
LIBS=-ldl -lz -lm -pthread -lcurl
LIBS+=-lcrypto

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

PROG=giggle

LIST=$(addprefix $(BIN)/, $(PROG))

all: check-env $(LIST) library

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
		$(LIBS) $(LDFLAGS)

check-env:
ifdef TRAVIS_BUILD_DIR
HTS_ROOT=$(TRAVIS_BUILD_DIR)/htslib
endif

