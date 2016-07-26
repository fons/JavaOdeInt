CC = gcc 
DEBUG    = -g
INCFLAGS = -I $(PWD)/../../codepack/include
CFLAGS   = -c -Wall $(INCFLAGS) $(DEBUG) -O3 -ftree-vectorize -msse2 -ffast-math
LDFLAGS  = -L $(PWD)/../../../target/lib/darwin_universal -lcodepack -lodepack -lm -lgfortran


OBJDIR       = $(PWD)/../../../target/test/obj
EXEDIR       = $(PWD)/../../../target/test/bin


#
# directory name is executable name
#
EXE          = $(shell basename `pwd`)
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
