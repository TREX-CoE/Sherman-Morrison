#!/bin/bash
make
rm -v blocked kay cu
./test b > blocked && ./test k > kay && ./test c > cu && head -n 4 cu && awk 'NR==5' blocked && awk 'NR==5' kay && awk 'NR==5' cu 

