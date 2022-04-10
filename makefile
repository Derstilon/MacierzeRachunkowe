OLD  := MMultBlock
NEW  := MMultBlock
#
# sample makefile
#

CC         := gcc
LINKER     := $(CC)
CFLAGS     := -O2 -Wall -msse3 -g
LDFLAGS    := -lm

UTIL       := copy_matrix.o \
              compare_matrices.o \
              random_matrix.o \
              dclock.o \
              REF_MMult.o \
              print_matrix.o

TEST_OBJS  := test_MMult.o $(NEW).o 

%.o: %.c
	$(CC) $(CFLAGS) -c $< -o $@
%.o: %.c
	$(CC) $(CFLAGS) -c $< -o $@

all: 
	make clean;
	make test_MMult.x

test_MMult.x: $(TEST_OBJS) $(UTIL) parameters.h
	$(LINKER) $(TEST_OBJS) $(UTIL) $(LDFLAGS) \
        $(BLAS_LIB) -o $(TEST_BIN) $@ 

run:	
	make all
	export OMP_NUM_THREADS=1
	export GOTO_NUM_THREADS=1
	echo "version = '$(NEW)';" > output_$(NEW).m
	./test_MMult.x >> output_$(NEW).m
	cp output_$(OLD).m output_old.m
	cp output_$(NEW).m output_new.m

test_invers_matrix.x: $(UTIL) parameters.h invers_matrix.o test_invers_matrix.o MMultBlock.o
	$(LINKER) $(UTIL) invers_matrix.o test_invers_matrix.o MMultBlock.o $(LDFLAGS) \
        $(BLAS_LIB) -o $(TEST_BIN) $@ 

test_factorization.x: $(UTIL) parameters.h lu_factorization.o test_factorization.o invers_matrix.o MMultBlock.o
	$(LINKER) $(UTIL) lu_factorization.o test_factorization.o invers_matrix.o MMultBlock.o $(LDFLAGS) \
		$(BLAS_LIB) -o $(TEST_BIN) $@

run_invers:
	make clean;
	make test_invers_matrix.x
	./test_invers_matrix.x > inverse_matrix.m

run_factor:
	make clean;
	make test_factorization.x
	./test_factorization.x > factorization.m

clean:
	rm -f *.o *~ core *.x

cleanall:
	rm -f *.o *~ core *.x output*.m *.eps *.png
