CXX=icpc
CXXFLAGS=-O0 -debug full -traceback
FC=ifort
FFLAGS=-O0 -debug full -traceback
# ARCH=-xCORE-AVX2 

## Deps & objs for the C++ stuff
cppDEPS = cppmain.cpp SM_MaponiA3.cpp SM_MaponiA3.hpp Helpers.hpp
cppOBJ = cppmain.o SM_MaponiA3.o
## Deps & objs for the Fortran stuff
fDEPS = fmain.f90 SM_MaponiA3_mod.f90
fOBJ = SM_MaponiA3_mod.o fmain.o

%.o: %.cpp $(cppDEPS)
	$(CXX) $(ARCH) $(CXXFLAGS) -c -o $@ $<

%.o: %.f90 $(fDEPS)
	$(FC) $(ARCH) $(FFLAGS) -c -o $@ $<

all: cppSherman-Morrison fSherman-Morrison

cppSherman-Morrison: $(cppOBJ)
	$(CXX) $(ARCH) $(CXXFLAGS) -o $@ $^

fSherman-Morrison: $(fOBJ)
	$(FC) $(ARCH) $(FFLAGS) -o $@ $^

clean:
	@rm -vf *.o *.mod
	