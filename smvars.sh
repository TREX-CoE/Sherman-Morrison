#!/usr/bin/env bash
unset THRESHOLD
unset MKL
ENV=$1

## Set Sherman-Morrison root dir
PWD=$(pwd)
SRCDIR=$(dirname $BASH_SOURCE)
case $SRCDIR in
	/*) SMROOT=$SRCDIR ;;
	*)
		if [[ $SRCDIR = . ]]
		then
			SMROOT=$PWD
		else
			SMROOT=$PWD/$SRCDIR
		fi
		;;
esac
export SMROOT

## Set environment for hdf5-tools and Makefile
case $ENV in
  intel)
	echo "* SM build environment set to 'intel'"
	export HDF5_CXX=icpc
	export HDF5_CXXLINKER=icpc
	export ENV=INTEL
    ;;
  llvm)
  	echo "* SM build environment set to 'llvm'"
	export HDF5_CXX=clang++
	export HDF5_CXXLINKER=clang++
	export ENV=LLVM
    ;;
  vfc)
  	echo "* SM build environment set to 'vfc'"
	export HDF5_CXX=clang++
	export HDF5_CXXLINKER=clang++
	export ENV=LLVM
	export VFC_BACKENDS="libinterflop_ieee.so --count-op"
    ;;
  gnu)
  	echo "* SM build environment set to 'gnu'"
	export HDF5_CXX=g++
	export HDF5_CXXLINKER=g++
	export ENV=GNU
	;;
  *)
    echo "Unknown environment descriptor given."
	echo "Usage: source smvars.sh {intel | llvm | vfc | gnu}"
	echo "Usage: source smvars.sh {intel | llvm | vfc | gnu} [mkl]"
	echo "Usage: source smvars.sh {intel | llvm | vfc | gnu} [threshold]"
	echo "Usage: source smvars.sh {intel | llvm | vfc | gnu} [mkl] [threshold]"
	echo "Usage: source smvars.sh {intel | llvm | vfc | gnu} [threshold] [mkl]"
	return 1
	;;
esac

## Export path, but only once
if [[ -z $SMVARS ]]
then
	export PATH=$SMROOT/bin:$PATH
	export SMVARS=true
fi

## Argument parsing
if [[ $# -eq 1 ]]
then
	echo "* Default threshold of '1e-3' will be used in next build."
fi
if [[ $# -eq 2 ]]
then
	if [[ $2 == mkl ]]
	then
		echo "* oneMKL-dependent parts of the code are enable in next build."
		echo "* Default threshold of '1e-3' will be used in next build."
		export MKL="-DMKL"
	else
		echo "* Threshold of '$2' will be used in next build."
		export THRESHOLD="-DTHRESHOLD=$2"
	fi
fi
if [[ $# -eq 3 ]]
then
	echo "* oneMKL-dependent parts of the code are enable in next build."
	export MKL="-DMKL"
	if [[ $2 == mkl ]]
	then
		echo "* Threshold of '$3' will be used in next build."
		export THRESHOLD="-DTHRESHOLD=$3"
	else
		echo "* Threshold of '$2' will be used in next build."
		export THRESHOLD="-DTHRESHOLD=$2"	
	fi
fi
