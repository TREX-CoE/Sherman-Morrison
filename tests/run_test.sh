#!/usr/bin/env bash

## IMPORTANT: start cycle will always be 1 because of the intention of this
## script and the way it is written now. In fact, var. $START is not used at
## the moment.

## Call: $ ./run_test.sh <{sm1|sm2|sm3|maponia3}> <start cycle> <stop cycle>
## Example: ./run_test.sh sm2 1 500

TEST=test_h5
ALGO=$1
START=$2
STOP=$3
BLOCKSIZE=329
TOLERANCE=1e-3

BLOCKS=$((STOP / BLOCKSIZE))
LEFT=$((STOP % BLOCKSIZE))

if [[ $LEFT -ne 0 ]]
then
    BSTART=0
    BSTOP=0
    LAST=0
    for ((i=1; i<=$BLOCKS; i++))
    do
        BSTART=$(((i-1)*BLOCKSIZE+1))
        BSTOP=$((i*BLOCKSIZE-1))
        $TEST $ALGO $BSTART $BSTOP $TOLERANCE
        LAST=$i
    done
    LSTART=$((LAST*BLOCKSIZE+1))
    LSTOP=$((LSTART+LEFT-1))
    $TEST $ALGO $LSTART $LSTOP $TOLERANCE
else
    for ((i=1; i<=$BLOCKS; i++))
    do
        BSTART=$(((i-1)*BLOCKSIZE+1))
        BSTOP=$((i*BLOCKSIZE-1))
        $TEST $ALGO $BSTART $BSTOP $TOLERANCE
    done
fi
