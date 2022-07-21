# FC = gfortran
# CC = gcc
# FFLAGS=-O0 -finline -g -lm -Wall -pedantic
# CFLAGS=-std=c99 -O0 -finline -g -lm -Wall -pedantic
FC = ifx
CC = icx
# FFLAGS=-O0 -warn all -g -pedantic
# CFLAGS=-std=c99 -O0 -Wall -g -pedantic
FFLAGS=-O3 -warn all -finline -xCORE-AVX2 -g -qopenmp -fopenmp-targets=spir64
CFLAGS=-std=c99 -O3 -Wall -finline -xCORE-AVX2 -g -qopenmp -fopenmp-targets=spir64
INCLUDE=-I/usr/include/hdf5/serial -I/usr/local/cuda/include
LFLAGS=-L/usr/lib/x86_64-linux-gnu/hdf5/serial -lhdf5 -lhdf5_hl -qmkl=sequential -L/usr/local/cuda-11.7/targets/x86_64-linux/lib -lcublas
#FC = verificarlo-f
#CC = verificarlo-c
#FFLAGS=-O3 -finline -g
#CFLAGS=-O3 -finline -g

## Link with icc
# test: sm.o test.o detupdate21.o meuk.o
# 	$(CC) $(LFLAGS) -o test sm.o detupdate21.o test.o meuk.o
test: sm.o test.o meuk.o
	$(CC) $(LFLAGS) -o test sm.o test.o meuk.o



## Link with ifort
# test: sm.o test.o detupdate21.o meuk.o
# 	$(FC) $(LFLAGS) -nofor-main -o test sm.o detupdate21.o test.o meuk.o

## Link with gfortran
# test: sm.o test.o detupdate21.o meuk.o
# 	$(FC) $(LFLAGS) -Wno-main -o test sm.o detupdate21.o test.o meuk.o

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
