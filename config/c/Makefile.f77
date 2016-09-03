##########################
OS = $(shell uname)
ifeq ($(OS), Darwin)
OS_ARCH = darwin_universal
SHARED_EXT = dylib
else
OS_ARCH = linux_x64
SHARED_EXT = so
endif

OUT_BASE_DIR = $(PWD)/../../../target/
OBJ_DIR      = $(OUT_BASE_DIR)/obj/$(OS_ARCH)
LIB_DIR      = $(OUT_BASE_DIR)/lib/$(OS_ARCH)

FC      = gfortran 
FCFLAGS = -fPIC -std=legacy -ftree-vectorize -msse2  -ffast-math
LDFLAGS =  -lm -lgfortran $(LIBS)

SOURCES  = $(wildcard *.f)
OBJECTS  = $(addprefix $(OBJ_DIR)/, $(SOURCES:.f=.o))
LIB      =  $(LIB_DIR)/lib$(LIBNAME).$(SHARED_EXT) 


all: $(LIB)

$(LIB) : $(OBJECTS)
	@mkdir -p $(LIB_DIR)
	$(CC)  -L /usr/local/gfortran/lib/ -L $(LIB_DIR) $(LDFLAGS) -dynamiclib -Wl  -o $(LIB) $(OBJECTS)

$(OBJ_DIR)/%.o: %.f
	@mkdir -p $(OBJ_DIR)
	$(FC) $(FCFLAGS) -c $< -o $@

clean:
	- rm $(LIB)
	- rm $(OBJECTS)
        
