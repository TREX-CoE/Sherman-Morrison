#!/bin/bash

INDIR="cycles_329_dets"
# INDIR="cycles_15784_dets"
OUTDIR=$1

mkdir -v ${OUTDIR}
cp -av get_stats.m ${OUTDIR}

## All cycles
ln -svf ${INDIR}/all_cycles.h cycles.h
make clean && make
./test a > ${OUTDIR}/ANTHONY.dat
./test n > ${OUTDIR}/NAIVE.dat
./test l > ${OUTDIR}/LATER.dat
./test s > ${OUTDIR}/SPLITTING.dat
./test b > ${OUTDIR}/BLOCKED.dat
./test m > ${OUTDIR}/MKL_LAPACK.dat

## Cycles w/ 2 upds excl. w/ WB2
ln -svf ${INDIR}/2_cycles.h cycles.h
make clean && make
./test 2 > ${OUTDIR}/WB2.dat

## Cycles w/ 3 upds excl. w/ WB3
ln -svf ${INDIR}/3_cycles.h cycles.h
make clean && make
./test 3 > ${OUTDIR}/WB3.dat

make clean
rm cycles.h

(cd ${OUTDIR} && ./get_stats.m)
