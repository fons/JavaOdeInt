CC = gcc 
DEBUG    = -g
INCFLAGS = -I $(PWD)/../../codepack/include
CFLAGS   = -c -Wall $(INCFLAGS) $(DEBUG) -O3 -ftree-vectorize -msse2 -ffast-math
LDFLAGS  = -L $(PWD)/../../../target/lib/darwin_universal 

EXE          = f77split
OBJDIR       = $(PWD)/../../../target/util/obj/$(EXE)
EXEDIR       = $(PWD)/../../../target/util/bin

#
# directory name is executable name
#

SOURCES      = $(wildcard *.c)
OBJECTS      = $(addprefix $(OBJDIR)/, $(SOURCES:.c=.o))
EXECUTABLE   = $(addprefix $(EXEDIR)/, $(EXE))

all : $(EXECUTABLE)

$(EXECUTABLE): $(OBJECTS) 
	@mkdir -p $(EXEDIR)
	$(CC)  $(OBJECTS) -o $@ $(LDFLAGS)

$(OBJDIR)/%.o: %.c
	@mkdir -p $(OBJDIR)
	$(CC) $(CFLAGS) $< -o $@

clean:
	- rm -f $(OBJECTS) *~ $(EXECUTABLE)
