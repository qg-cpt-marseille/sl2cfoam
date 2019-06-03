QUIET ?= @

CC = gcc
AR = ar

default: lib test
all: default

# folder names
EXTDIR = ext
SRCDIR = src
OBJDIR = obj
INCDIR = inc
LIBDIR = lib
BINDIR = bin
TESTDIR = test

# two sets of CFLAGS
#CFLAGS = -std=c99 -fopenmp -fPIC -g -O0 -I$(EXTDIR)/wigxjpf-1.9/inc/ -I$(EXTDIR)/ -Iinc/ -Isrc/ # debug
CFLAGS = -std=c99 -fopenmp -fPIC -O2 -I$(EXTDIR)/wigxjpf-1.9/inc/ -I$(EXTDIR)/ -Iinc/ -Isrc/ # optimized

LDFLAGS = -L$(EXTDIR)/wigxjpf-1.9/lib/
LDLIBS =  -lm -lgsl -lgslcblas -lgmp -lmpfr -lmpc -lwigxjpf
ARFLAGS = rcs

INCS = $(EXTDIR)/khash.h src/utilities.h src/jsymbols.h src/common.h src/b4function.h         \
       src/config.h src/coherentstates.h src/cgamma.h src/b3function.h src/recouplingSL2C.h   \
       src/hashing.h src/recursion.h

_OBJS = utilities.o khash.o jsymbols.o coherentstates.o b4function.o b3function.o cgamma.o    \
        recouplingSL2C.o init.o ampl_3simplex.o ampl_4simplex.o config.o recursion.o          \
	hash_3simplex.o	hash_4simplex.o

OBJS = $(patsubst %,$(OBJDIR)/%,$(_OBJS))

_TESTS = 3simplex_test 4simplex_test b3_test b4_test recoupling_test coherent_test bf_test
TESTS = $(patsubst %,$(BINDIR)/%,$(_TESTS))

# khash object file
$(OBJDIR)/khash.o: $(EXTDIR)/khash.h $(INCS)
	@echo "   CC    $@"
	$(QUIET)mkdir -p $(dir $@)
	$(QUIET)$(CC) $(CFLAGS) -c -o $@ $< 

# library/src object files
$(OBJDIR)/%.o: $(SRCDIR)/%.c $(INCS)
	@echo "   CC    $@"
	$(QUIET)mkdir -p $(dir $@)
	$(QUIET)$(CC) $(CFLAGS) -c -o $@ $< 

# build library
$(LIBDIR)/libsl2cfoam.a: $(OBJS)
	@echo "   AR    $@"
	$(QUIET)mkdir -p $(dir $@)
	$(QUIET)$(AR) $(ARFLAGS) $@ $(OBJS)

# compile test programs
$(BINDIR)/%: $(TESTDIR)/%.c $(LIBDIR)/libsl2cfoam.a 
	@echo "   CC    $@"
	$(QUIET)mkdir -p $(dir $@)
	$(QUIET)$(CC) $(CFLAGS) -o $@ -Iinc/ $< -Llib/ $(LDFLAGS) -lsl2cfoam $(LDLIBS)

.PHONY: default all clean

lib: $(LIBDIR)/libsl2cfoam.a

test: lib $(TESTS)
	
clean:
	rm -rf $(OBJDIR)
	rm -rf $(LIBDIR)
	rm -rf $(BINDIR)
