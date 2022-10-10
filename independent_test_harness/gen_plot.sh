#!/bin/bash

for SIZE in 32 64 128 256 512 1024 2048 4096 8192 16384
do
    echo $SIZE >> SIZES
    for LOAD in 25 50 75 100
    do
        NUPDS=$((SIZE*LOAD/100))
        ln -svf ../../random_matrix_generator/ds_${SIZE}_${NUPDS}.hdf5 dataset
        for KERNEL in MKL WBK_CPU WBK_GPU
        do
            case $KERNEL in
                MKL)
                    ./test_nvc_ompol m | awk 'NR==5 {print $11}' >> ${KERNEL}_${LOAD}.dat
                    ;;
                WBK_CPU)
                    ./test_nvc_ompol k | awk 'NR==5 {print $11}' >> ${KERNEL}_${LOAD}.dat
                    ;;
                WBK_GPU)
                    ./test_nvc_ompol c | awk 'NR==5 {print $11}' >> ${KERNEL}_${LOAD}.dat
                    ;;
            esac
        done
    done
done

for KERNEL in MKL WBK_CPU WBK_GPU
do
    paste SIZES ${KERNEL}_25.dat ${KERNEL}_50.dat ${KERNEL}_75.dat ${KERNEL}_100.dat > ${KERNEL}.dat
    rm -v ${KERNEL}_25.dat ${KERNEL}_50.dat ${KERNEL}_75.dat ${KERNEL}_100.dat
done
rm -v SIZES
