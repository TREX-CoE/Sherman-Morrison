#!/usr/bin/env bash

STYLE='--style=LLVM'
FORMATER='clang-format -i'

if [[ -z $SMVARS ]]
then
	echo '$SMVARS is not set. Please source '/path/to/Sherman-Morrison/smvars.sh''
	exit 1
fi

for ext in c cc cpp h hpp
do
    find $SMROOT -type f -iname "*.${ext}" -exec echo "$FORMATER $STYLE" {} \;
done
