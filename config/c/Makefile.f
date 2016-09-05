PICFLAG   = -fPIC
DEBUGFLAG = -g
OPTFLAG   = -O3

##########################
OS = $(shell uname)


ifeq ($(OS), Darwin)
  OS_ARCH = darwin_universal
  LDFLAGS   += -mmacosx-version-min=10.11
  SHARED_EXT = dylib
  SHARED_FLAG=-dynamiclib
 else
     OS_ARCH = linux_x64
     SHARED_EXT = so
     LDFLAGS   += -Wl,--hash-style=both
     SHARED_FLAG=-shared
endif

OUT_BASE_DIR = $(PWD)/../../../target/

LIB_DIR      = $(OUT_BASE_DIR)/lib/$(OS_ARCH)
OBJ_DIR      = $(OUT_BASE_DIR)/obj/$(LIBNAME)/$(OS_ARCH)

SOURCES      = $(wildcard *.c)
INCLUDES     = $(wildcard ../include/*.h)
OBJECTS      = $(addprefix $(OBJ_DIR)/, $(SOURCES:.c=.o))

CC     = gcc

CFLAGS = -c -Wall $(INCFLAG) $(DEBUGFLAG) $(PICFLAG) $(OPTFLAG) -ftree-vectorize -msse2  -ffast-math -std=gnu99

LIB    = $(LIB_DIR)/lib$(LIBNAME).$(SHARED_EXT)


############################

all: $(LIB) 

$(LIB) : $(OBJECTS)
	@mkdir -p $(LIB_DIR)
	$(CC) -L$(LIB_DIR)  $(LINFLAGS) $(SHARED_FLAG) $(LDFLAGS) -o $(LIB) $(OBJECTS)

$(OBJ_DIR)/%.o: %.c $(INCLUDES)
	@mkdir -p $(OBJ_DIR)
	$(CC) $(CFLAGS)  $< -o $@

clean:
	- rm $(OBJECTS)
	- rm $(LIB)

