#!/usr/bin/env bash

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
	export HDF5_CXX=icpc
	export HDF5_CXXLINKER=icpc
	export ENV=INTEL
    ;;
  llvm)
	export HDF5_CXX=clang++
	export HDF5_CXXLINKER=clang++
	export ENV=LLVM
    ;;
  vfc)
	export HDF5_CXX=clang++
	export HDF5_CXXLINKER=clang++
	export ENV=LLVM
	export VFC_BACKENDS="libinterflop_ieee.so --count-op"
    ;;
  gnu)
	export HDF5_CXX=g++
	export HDF5_CXXLINKER=g++
	export ENV=GNU
	;;
  *)
    echo "Unknown environment descriptor given."
	echo "Usage: source smvars.sh {intel | llvm | vfc | gnu}"
	return 1
	;;
esac

## Export path, but only once
if [[ -z $SMVARS ]]
then
	export PATH=$SMROOT/bin:$PATH
	export SMVARS=true
fi
