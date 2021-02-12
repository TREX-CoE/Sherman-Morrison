## Used compilers
CXX = icpc
FC = ifort

## Compiler flags
CXXFLAGS = -O0 #-debug full -traceback
FFLAGS = -O0 #-debug full -traceback
# ARCH = -xCORE-AVX2

## Deps & objs for C++ cMaponiA3_test
cMaponiA3_testDEP = cMaponiA3_test.cpp SM_MaponiA3.cpp SM_MaponiA3.hpp Helpers.hpp
cMaponiA3_testOBJ = cMaponiA3_test.o SM_MaponiA3.o

## Deps & objs for Fortran fMaponiA3_test
fMaponiA3_testDEP = fMaponiA3_test.f90 SM_MaponiA3_mod.f90
fMaponiA3_testOBJ = SM_MaponiA3.o SM_MaponiA3_mod.o fMaponiA3_test.o
fMaponiA3_testLIB = -lstdc++

## Deps & objs for Fortran QMCChem_dataset_test
QMCChem_dataset_testDEP = QMCChem_dataset_test.f90 SM_MaponiA3_mod.f90 Utils_mod.f90
QMCChem_dataset_testOBJ = SM_MaponiA3.o Utils_mod.o SM_MaponiA3_mod.o QMCChem_dataset_test.o
QMCChem_dataset_testLIB = -lstdc++

## Compile recipes for C++ cMaponiA3_test
%.o: %.cpp $(cMaponiA3_testDEP)
	$(CXX) $(ARCH) $(CXXFLAGS) -c -o $@ $<

## Compile recepies for Fortran fMaponiA3_test
%.o: %.f90 $(fMaponiA3_testDEP)
	$(FC) $(ARCH) $(FFLAGS) -c -o $@ $<

## Build tagets
.PHONY: all clean distclean

all: cMaponiA3_test fMaponiA3_test QMCChem_dataset_test

clean:
	@rm -vf *.o *.mod

distclean: clean
	@rm -vf cMaponiA3_test fMaponiA3_test QMCChem_dataset_test

## Linking the C++ example program
cMaponiA3_test: $(cMaponiA3_testOBJ)
	$(CXX) $(ARCH) $(CXXFLAGS) -o $@ $^

## Linking Fortran example program calling the C++ function 'Sherman_Morrison()'
fMaponiA3_test: $(fMaponiA3_testOBJ)
	$(FC) $(ARCH) $(FFLAGS) $(fMaponiA3_testLIB) -o $@ $^

## Linking Fortran example program calling the C++ function 'Sherman_Morrison()'
QMCChem_dataset_test: $(QMCChem_dataset_testOBJ)
	$(FC) $(ARCH) $(FFLAGS) $(QMCChem_dataset_testLIB) -o $@ $^
