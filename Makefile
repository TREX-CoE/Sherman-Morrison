## Used compilers
CXX = icpc
FC = ifort

## Compiler flags
CXXFLAGS = -O0 -debug full -traceback
FFLAGS = -O0 -debug full -traceback
# ARCH = -xCORE-AVX2 

## Deps & objs for the C++ stuff
cppDEPS = cppmain.cpp SM_MaponiA3.cpp SM_MaponiA3.hpp Helpers.hpp
cppOBJ = cppmain.o SM_MaponiA3.o

## Deps & objs for the Fortran stuff
fDEPS = fmain.f90 SM_MaponiA3_mod.f90
fOBJ = SM_MaponiA3_f.o SM_MaponiA3_mod.o fmain.o
fLIBS = -lstdc++

## Compile recipes for C++ stuff
%.o: %.cpp $(cppDEPS)
	$(CXX) $(ARCH) $(CXXFLAGS) -c -o $@ $<

## Compile recepies for Fortran stuff
%.o: %.f90 $(fDEPS)
	$(FC) $(ARCH) $(FFLAGS) -c -o $@ $<

## Build tagets
.PHONY: all clean distclean

all: cppSherman-Morrison fSherman-Morrison

clean:
	@rm -vf *.o *.mod
	
distclean: clean
	@rm -vf cppSherman-Morrison fSherman-Morrison

## Linking the C++ example program
cppSherman-Morrison: $(cppOBJ)
	$(CXX) $(ARCH) $(CXXFLAGS) -o $@ $^

## Linking Fortran example program calling the C++ function 'Sherman_Morrison()'
fSherman-Morrison: $(fOBJ)
	$(FC) $(ARCH) $(FFLAGS) $(fLIBS) -o $@ $^
