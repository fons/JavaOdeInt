PICFLAG   = -fPIC
DEBUGFLAG = -g
OPTFLAG   = -O3

OUT_BASE_DIR = $(PWD)/../../../target/
LIB_DIR      = $(OUT_BASE_DIR)/lib/$(OS_ARCH)
OBJ_DIR      = $(OUT_BASE_DIR)/obj/$(LIBNAME)/$(OS_ARCH)

##########################
OS = $(shell uname)

ifeq ($(OS), Darwin)
  OS_ARCH = darwin_universal
  LDFLAGS   += -mmacosx-version-min=10.11
  SHARED_EXT = dylib
  SHARED_FLAG=-dynamiclib
  LD=gcc
 else
     OS_ARCH = linux_x64
     SHARED_EXT = so
#     LDFLAGS   += -Wl,--hash-style=both,-rpath-link,$(LIB_DIR)
     LDFLAGS   = --hash-style=both -rpath-link=$(LIB_DIR)
     SHARED_FLAG=-shared
     LD=ld
     F77_LIB_DIR_OPTION=-L /usr/lib/gcc/x86_64-linux-gnu/4.8  
     F77_LIB_OPTION=-lgfortran
endif

SOURCES   = $(wildcard *.c)
INCLUDES  = $(wildcard ../include/*.h)
OBJECTS   = $(addprefix $(OBJ_DIR)/, $(SOURCES:.c=.o))

CC     = gcc

CFLAGS = -c -Wall $(INCFLAG) $(DEBUGFLAG) $(PICFLAG) $(OPTFLAG) -ftree-vectorize -msse2  -ffast-math -std=gnu99

LIB    = $(LIB_DIR)/lib$(LIBNAME).$(SHARED_EXT)

############################

all: $(LIB) 

$(LIB) : $(OBJECTS)
	@mkdir -p $(LIB_DIR)
	$(LD) $(LDFLAGS2) -L$(LIB_DIR) $(F77_LIB_DIR_OPTION) $(F77_LIB_OPTION) $(LINFLAGS) $(SHARED_FLAG) -o $(LIB) $(OBJECTS)
#	ld $(LDFLAGS2) $(LINFLAGS) -L$(LIB_DIR) $(SHARED_FLAG) -o $(LIB) $(OBJECTS)
#	$(CC) $(LDFLAGS) -L$(LIB_DIR)  $(LINFLAGS) $(SHARED_FLAG) $(LDFLAGS) -o $(LIB) $(OBJECTS)
#	$(CC) -L$(LIB_DIR)  $(LINFLAGS) $(SHARED_FLAG) $(LDFLAGS) -o $(LIB) $(OBJECTS)

$(OBJ_DIR)/%.o: %.c $(INCLUDES)
	@mkdir -p $(OBJ_DIR)
	$(CC) $(CFLAGS)  $< -o $@

clean:
	- rm $(OBJECTS)
	- rm $(LIB)

