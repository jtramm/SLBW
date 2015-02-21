COMPILER    = gnu

# Standard Flags
CFLAGS := -std=gnu99

# Linker Flags
LDFLAGS = -lm

# Regular gcc Compiler
ifeq ($(COMPILER),gnu)
  CC = gcc
  LDFLAGS += -fopenmp
  CFLAGS += -Ofast -ffast-math
endif

# intel Compiler
ifeq ($(COMPILER),intel)
  CC = icc
  LDFLAGS += -openmp
  CFLAGS += -O3 -xhost -ansi-alias -no-prec-div -DINTEL -vec-report6
endif

all:
	$(CC) $(CFLAGS) main.c Faddeeva.c -o SLBW $(LDFLAGS)

clean:
	rm -f SLBW data.dat
run:
	./SLBW
graph:
	gnuplot graph.gp
edit:
	vim main.c
