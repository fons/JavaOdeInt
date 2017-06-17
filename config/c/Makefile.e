##########################
CC = gcc 

OUT_BASE_DIR = $(PWD)/../../../target/

LIB_DIR      = $(OUT_BASE_DIR)/lib/$(OS_ARCH)
OBJ_DIR      = $(OUT_BASE_DIR)/test/$(OS_ARCH)/obj
EXE_DIR      = $(OUT_BASE_DIR)/test/$(OS_ARCH)/bin


OS = $(shell uname)

ifeq ($(OS), Darwin)
  OS_ARCH = darwin_universal
  LDFLAGS   += -mmacosx-version-min=10.12
  SHARED_EXT = dylib
  SHARED_FLAG=-dynamiclib
 else
     OS_ARCH = linux_x64
     SHARED_EXT = so
     LDFLAGS   += -Wl,--hash-style=both,-rpath=$(LIB_DIR)
     SHARED_FLAG=-shared
endif


#
# directory name is executable name
#
EXE          = $(shell basename `pwd`)
SOURCES      = $(wildcard *.c)
INCLUDES     = $(wildcard ../include/*.h)
OBJECTS      = $(addprefix $(OBJ_DIR)/, $(SOURCES:.c=.o))
EXECUTABLE   = $(addprefix $(EXE_DIR)/, $(EXE))

#######################

all : $(EXECUTABLE)

$(EXECUTABLE): $(OBJECTS) 
	@mkdir -p $(EXE_DIR)
	$(CC) $(LDFLAGS)  $(OBJECTS)  -o $@ -L $(LIB_DIR) $(LIBS) -lm 

$(OBJ_DIR)/%.o: %.c
	@mkdir -p $(OBJ_DIR)
	$(CC) -c -Wall $(INCFLAG) $(DEBUGFLAG) -O3 -ftree-vectorize -msse2 -ffast-math  $< -o $@

clean:
	- rm -f $(OBJECTS) *~ $(EXECUTABLE)
