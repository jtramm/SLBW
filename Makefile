COMPILER    = gnu

program = SLBW

source = \
main.c \
Faddeeva.c

obj = $(source:.c=.o)

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

$(program): $(obj)
	$(CC) $(CFLAGS) $(obj) -o $@ $(LDFLAGS)

%.o: %.c
	$(CC) $(CFLAGS) -c $< -o $@

clean:
	rm -f SLBW data.dat main.o Faddeeva.o
run:
	./SLBW
graph:
	gnuplot graph.gp
edit:
	vim main.c
