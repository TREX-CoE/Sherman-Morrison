#!/bin/bash

if [[ $1 == "cpu"  ]]
then
	TEST=./test_icc_mkl
	MAKEFILE=Makefile.icc_mkl.cpu
elif [[ $1 == "gpu" ]]
then
	TEST=./test_nvc_ompol
	MAKEFILE=Makefile.nvc_ompol.gpu
else
	echo -e "Need either 'cpu' or gpu'\n"
	exit 1
fi

rm -vf kay cu em
make -f $MAKEFILE clean all
time $TEST m > em
time $TEST k > kay
[[ $1 == "gpu" ]] && time $TEST c > cu
head -n 4 em
awk 'NR==5' em
awk 'NR==5' kay
[[ $1 == "gpu" ]] && awk 'NR==5' cu

