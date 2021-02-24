# ARCH = -xCORE-AVX2

## Used compilers
H5CXX = h5c++
CXX = icpc
FC = ifort

## Compiler flags & common obs & libs
H5CXXFLAGS = -O0 -g
CXXFLAGS = -O0 -g -traceback
FFLAGS = -O0 -g -traceback
FLIBS = -lstdc++
OBJS = SM_MaponiA3.o

## Deps & objs for C++ cMaponiA3_test_3x3_3
cMaponiA3_test_3x3_3OBJ = cMaponiA3_test_3x3_3.o
fMaponiA3_test_3x3_3OBJ = SM_MaponiA3_mod.o fMaponiA3_test_3x3_3.o
fMaponiA3_test_4x4_2OBJ = Helpers_mod.o SM_MaponiA3_mod.o fMaponiA3_test_4x4_2.o
QMCChem_dataset_testOBJ = Helpers_mod.o SM_MaponiA3_mod.o QMCChem_dataset_test.o


## Default build target: build everything
all: cMaponiA3_test_3x3_3 fMaponiA3_test_3x3_3 fMaponiA3_test_4x4_2 QMCChem_dataset_test tests/test


## Compile recipes for C++
%.o: %.cpp
	$(CXX) $(ARCH) $(CXXFLAGS) -c -o $@ $<

## Compile recepies for Fortran
%.o: %.f90
	$(FC) $(ARCH) $(FFLAGS) -c -o $@ $<

## Explicit recipe to trigger rebuild and relinking when headerfile is changed
SM_MaponiA3.o: SM_MaponiA3.cpp Helpers.hpp
	$(CXX) $(ARCH) $(CXXFLAGS) -fPIC -c -o $@ $<


## Build tagets
.PHONY: all clean distclean

clean:
	@rm -vf *.o *.mod

distclean: clean
	@rm -vf cMaponiA3_test_3x3_3 \
	fMaponiA3_test_3x3_3 fMaponiA3_test_4x4_2 \
	QMCChem_dataset_test \
	Slater* Updates.dat \
	tests/test


## Linking the C++ example program
cMaponiA3_test_3x3_3: $(cMaponiA3_test_3x3_3OBJ) $(OBJS)
	$(CXX) $(ARCH) $(CXXFLAGS) -o $@ $^

## Linking Fortran example program calling the C++ function
fMaponiA3_test_3x3_3: $(fMaponiA3_test_3x3_3OBJ) $(OBJS)
	$(FC) $(ARCH) $(FFLAGS) $(FLIBS) -o $@ $^

fMaponiA3_test_4x4_2: $(fMaponiA3_test_4x4_2OBJ) $(OBJS)
	$(FC) $(ARCH) $(FFLAGS) $(FLIBS) -o $@ $^

QMCChem_dataset_test: $(QMCChem_dataset_testOBJ) $(OBJS)
	$(FC) $(ARCH) $(FFLAGS) $(FLIBS) -o $@ $^

tests/test: tests/test.cpp SM_MaponiA3.o
	$(H5CXX) $(ARCH) $(H5CXXFLAGS) -o $@ $^
