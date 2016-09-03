CC = gcc 

CFLAGS    = -c -Wall $(INCFLAG) $(DEBUGFLAG) -O3 -ftree-vectorize -msse2 -ffast-math
LDFLAGS   = -L $(PWD)/../../../target/lib/darwin_universal -lm -lgfortran $(LIBS)


OBJ_DIR      = $(PWD)/../../../target/test/obj
EXE_DIR      = $(PWD)/../../../target/test/bin


#
# directory name is executable name
#
EXE          = $(shell basename `pwd`)
SOURCES      = $(wildcard *.c)
OBJECTS      = $(addprefix $(OBJ_DIR)/, $(SOURCES:.c=.o))
EXECUTABLE   = $(addprefix $(EXE_DIR)/, $(EXE))

all : $(EXECUTABLE)

$(EXECUTABLE): $(OBJECTS) 
	@mkdir -p $(EXE_DIR)
	$(CC)  $(OBJECTS) -o $@ $(LDFLAGS)

$(OBJ_DIR)/%.o: %.c
	@mkdir -p $(OBJ_DIR)
	$(CC) $(CFLAGS) $< -o $@

clean:
	- rm -f $(OBJECTS) *~ $(EXECUTABLE)
