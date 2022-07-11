#!/bin/bash
IN=$1 # Input dataset (hdf5)
NU=$2 # Number of updates

CYCLE_LIST=$(h5ls ${IN} | awk '{print$1}' | sed 's/cycle_//g' | sort -n)
SELECTION=()

# Filter CYCLE_LIST and add to SELECTION
for CYCLE in ${CYCLE_LIST}
do
    NUPDS=$(h5ls -d ${IN}/cycle_${CYCLE}/nupdates | awk 'FNR == 3 {print $2}')
    if (( NUPDS == NU ))
    then
        SELECTION+=($CYCLE)
    fi
    # SELECTION+=($CYCLE)
done

# Generate C-header file
NELEMENTS=${#SELECTION[@]}
echo "const uint32_t n_cycles = $NELEMENTS;" > ${NU}_cycles.h
echo -n "uint32_t cycles[n_cycles] = {" >> ${NU}_cycles.h
for VAL in "${SELECTION[@]}"
do
    echo -n "$VAL," >> ${NU}_cycles.h
done
truncate -s-1 ${NU}_cycles.h # remove last ','
echo "};" >> ${NU}_cycles.h
