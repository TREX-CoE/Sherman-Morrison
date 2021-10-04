#!/usr/bin/env bash

## THIS SCRIPT ONLY WORKS IF QMCKL HAS BEEN BUILT CORRECTLY

gfortran -c qmckl_f.f90
gfortran -ffree-line-length-none -c qmckl_test_f.f90
gfortran -o qmckl_test_f qmckl_test_f.o -L. -lqmckl
export LD_LIBRARY_PATH=../../qmckl/src/.libs:/opt/intel/oneapi/compiler/2021.3.0/linux/compiler/lib/intel64_lin
#export LD_LIBRARY_PATH=$HOME/code/qmcchem/lib

echo
echo "Current '\$LD_LIBRARY_PATH':" $LD_LIBRARY_PATH
echo
./qmckl_test_f && octave qmckl_test_f.m
