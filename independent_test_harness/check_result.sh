#!/bin/bash

rm -v blocked kay cu em

make && 
	./test m > em && \
	./test b > blocked && \
	./test k > kay && \
	./test c > cu && \
	head -n 4 em && \
	awk 'NR==5' em && \
	awk 'NR==5' blocked && \
	awk 'NR==5' kay && \
	awk 'NR==5' cu

