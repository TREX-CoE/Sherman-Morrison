icpc -c worker.cpp && ifort -c main.f90 && ifort -lstdc++ worker.o main.o -o test
