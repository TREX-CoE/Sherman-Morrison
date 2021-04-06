#!/usr/bin/env bash

ENV=$1

PWD=$(pwd)
SRCDIR=$(dirname $BASH_SOURCE)
case $SRCDIR in
	/*) SMROOT=$SRCDIR ;; ## sourced from absolute path
	*) ## sourced from absolute path
		if [[ $SRCDIR = . ]] ## check if already in root
		then
			SMROOT=$PWD
		else
			SMROOT=$PWD/$SRCDIR
		fi
		;;
esac
export SMROOT

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
  gnu)
	export HDF5_CXX=g++
	export HDF5_CXXLINKER=g++
	export ENV=GNU
	;;
  *)
    echo "Unknown environment descriptor given."
	echo "Usage: source smvars.sh {intel | llvm | gnu}"
	return 1
	;;
esac

if [[ -z $SMVARS ]]
then
	export PATH=$SMROOT/bin:$PATH
	export SMVARS=true
fi
