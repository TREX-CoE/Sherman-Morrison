#!/usr/bin/env bash

start_cycle=$1
stop_cycle=$2

if [ "$#" -ne 2 ]
then
    echo "usage: ./run-tests.sh <start cycle> <stop cycle>"
    exit 1
fi

if [ ! -f "bin/test_external_h5" ]
then
    make bin/test_external_h5
fi

cd datasets/

for ((cycle = start_cycle; cycle < stop_cycle+1; cycle++))
do
    ../bin/test_external_h5 $cycle
done