LIBNAME =  gamd
CC=gcc
##########################
OS = $(shell uname)
ifeq ($(OS), Darwin)
OS_ARCH = darwin_universal
SHARED_EXT = dylib
LIBARGS=-Wl
F77_LIB_DIR= -L /usr/local/gfortran/lib/
SHARED_FLAG=-dynamiclib
else
OS_ARCH = linux_x64
SHARED_EXT = so
SHARED_FLAG=-shared
LIBARGS=-Wl,-soname,lib$(LIBNAME).$(SHARED_EXT) 
endif
#-Wl,-soname,libopks.so

OUT_BASE_DIR = $(PWD)/../../../target/
OBJ_DIR      = $(OUT_BASE_DIR)/obj/$(OS_ARCH)
MOD_DIR      = $(OUT_BASE_DIR)/mod/$(OS_ARCH)
LIB_DIR      = $(OUT_BASE_DIR)/lib/$(OS_ARCH)

FC      = gfortran 
LDFLAGS =  -lm -lgfortran $(LIBS)

SOURCES  = $(wildcard *.f90)
OBJECTS  = $(addprefix $(OBJ_DIR)/, $(SOURCES:.f90=.o))
LIB      = $(LIB_DIR)/lib$(LIBNAME).$(SHARED_EXT) 

all: $(LIB)
#
# order is important; this forces compile order
#
$(LIB) : $(OBJ_DIR)/gamda.o $(OBJ_DIR)/gamd.o 
	@mkdir -p $(LIB_DIR)
	$(CC)  $(F77_LIB_DIR) -L $(LIB_DIR) $(LDFLAGS) $(SHARED_FLAG) $(LIBARGS)  -o $(LIB) $(OBJECTS)

$(OBJ_DIR)/%.o: %.f90
	@mkdir -p $(OBJ_DIR)
	@mkdir -p $(MOD_DIR)
	$(FC) $(FCFLAGS) -J$(MOD_DIR) -c $< -o $@

clean:
	- rm $(LIB)
	- rm $(OBJECTS)
