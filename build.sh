

mkdir build/ bin/

## compile tests/test.cpp
icpc -O0 -fPIC -c src/SM_MaponiA3.cpp -o build/SM_MaponiA3.o -I include/

## compile tests/fMaponiA3_test_3x3_3.f90
ifort -c src/SM_MaponiA3_mod.f90 -o build/SM_MaponiA3_mod.o -module build
ifort -c tests/fMaponiA3_test_3x3_3.f90 -o build/fMaponiA3_test_3x3_3.o -I build/
## link bin/QMCChem_dataset_test
ifort -o bin/fMaponiA3_test_3x3_3 build/fMaponiA3_test_3x3_3.o build/SM_MaponiA3_mod.o build/SM_MaponiA3.o -lstdc++

## compile tests/fMaponiA3_test_4x4_2.f90
ifort -c src/SM_MaponiA3_mod.f90 -o build/SM_MaponiA3_mod.o -module build
ifort -c src/Helpers_mod.f90 -o build/Helpers_mod.o -module build
ifort -c tests/fMaponiA3_test_4x4_2.f90 -o build/fMaponiA3_test_4x4_2.o -I build/
## link bin/QMCChem_dataset_test
ifort -o bin/fMaponiA3_test_4x4_2 build/fMaponiA3_test_4x4_2.o build/SM_MaponiA3_mod.o build/Helpers_mod.o build/SM_MaponiA3.o -lstdc++

## compile tests/QMCChem_dataset_test.f90
ifort -c src/SM_MaponiA3_mod.f90 -o build/SM_MaponiA3_mod.o -module build
ifort -c src/Helpers_mod.f90 -o build/Helpers_mod.o -module build
ifort -c tests/QMCChem_dataset_test.f90 -o build/QMCChem_dataset_test.o -I build/
## link bin/QMCChem_dataset_test
ifort -o bin/QMCChem_dataset_test build/QMCChem_dataset_test.o build/SM_MaponiA3_mod.o build/Helpers_mod.o build/SM_MaponiA3.o  -lstdc++

## compile tests/test.cpp
h5c++ -c tests/test.cpp -o build/test.o -I include/
## link bin/test
h5c++ -o bin/test build/test.o build/SM_MaponiA3.o

## compile tests/cMaponiA3_test_3x3_3.cpp
icpc -c tests/cMaponiA3_test_3x3_3.cpp -o build/cMaponiA3_test_3x3_3.o -I include/
## link bin/cMaponiA3_test_3x3_3
icpc -o bin/cMaponiA3_test_3x3_3 build/cMaponiA3_test_3x3_3.o build/SM_MaponiA3.o
