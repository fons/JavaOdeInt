LIBNAME = cdvode


##########################
OS = $(shell uname)
ifeq ($(OS), linux)
endif

ifeq ($(OS), Darwin)
OS_ARCH = darwin_universal
CPPFLAGS  += -D_DARWIN_C_SOURCE
CPPFLAGS  += -mmacosx-version-min=10.11
LDFLAGS   += -mmacosx-version-min=10.11
SHARED_EXT = dylib
else
OS_ARCH = linux_x64
SHARED_EXT = so
LDFLAGS   += -Wl,--hash-style=both
endif


############################
OUT_BASE_DIR = $(PWD)/../../../target/
LIB_DIR      = $(OUT_BASE_DIR)/lib/$(OS_ARCH)
OBJ_DIR      = $(OUT_BASE_DIR)/obj/$(OS_ARCH)

SOURCES      = $(wildcard *.c)
INCLUDES     = $(wildcard ../include/*.h)
OBJECTS      = $(addprefix $(OBJ_DIR)/, $(SOURCES:.c=.o))

LINFLAGS =   -ldvode 

CC     = gcc 
CFLAGS = -Wall -g -fPIC -O3 -ftree-vectorize -msse2  -ffast-math -std=gnu99
LIB    = $(LIB_DIR)/lib$(LIBNAME).$(SHARED_EXT)

all: $(LIB) 

$(LIB) : $(OBJECTS)
	@mkdir -p $(LIB_DIR)
	$(CC) -L$(LIB_DIR)  $(LINFLAGS) -dynamiclib $(LDFLAGS) -o $(LIB) $(OBJECTS)

$(OBJ_DIR)/%.o: %.c $(INCLUDES)
	@mkdir -p $(OBJ_DIR)
	$(CC) -c $(CFLAGS)  $< -o $@

clean:
	- rm $(OBJECTS)
	- rm $(LIB)
