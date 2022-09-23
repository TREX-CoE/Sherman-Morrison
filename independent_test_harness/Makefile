#FC = ifx
#CC = nvc

CFLAGS=-std=c99 -O3 -Wall -g -mp -target=gpu

INCLUDE=-I$(NVHPC_ROOT)/math_libs/include

LDFLAGS=-L/usr/lib/x86_64-linux-gnu/hdf5/serial -lhdf5 -lhdf5_hl
LDFLAGS+=-L$(MKLROOT)/lib/intel64 -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lpthread -lm -ldl
LDFLAGS+=-L$(NVHPC_ROOT)/math_libs/lib64 -lcublas -lcusolver -mp -target=gpu

all: test

## Link with icc
# test: sm.o test.o detupdate21.o meuk.o
# 	$(CC) $(LDFLAGS) -o test sm.o detupdate21.o test.o meuk.o
test: sm.o test.o meuk.o
	$(CC) $(LDFLAGS) -o test sm.o test.o meuk.o

%.o: %.f90
	$(FC) $(FFLAGS) -c -o $@ $< 

%.o : %.c
	$(CC) $(CFLAGS) $(INCLUDE) -c -o $@ $<

clean:  
	rm -rf *.o *genmod* test test

debug_n_2: test
	gdb --silent --args ./test n 2

debug_a_2: test
	gdb --silent --args ./test a 2
