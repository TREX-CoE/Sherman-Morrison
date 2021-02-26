SRC_DIR := src
INC_DIR := include
OBJ_DIR := build
BIN_DIR := bin

## Used compilers
ARCH = -xCORE-AVX2
H5CXX = h5c++
CXX = icpc
FC = ifort

## Compiler flags & common obs & libs
H5CXXFLAGS = -O0 -g
CXXFLAGS = -O0 -g -traceback
FFLAGS = -O0 -g -traceback
INCLUDE = -I$(INC_DIR)
FLIBS = -lstdc++

## Deps & objs for C++ cMaponiA3_test_3x3_3
OBJS = $(OBJ_DIR)/SM_MaponiA3.o
cMaponiA3_test_3x3_3OBJ = cMaponiA3_test_3x3_3.o
fMaponiA3_test_3x3_3OBJ = SM_MaponiA3_mod.o fMaponiA3_test_3x3_3.o
fMaponiA3_test_4x4_2OBJ = Helpers_mod.o SM_MaponiA3_mod.o fMaponiA3_test_4x4_2.o
QMCChem_dataset_testOBJ = Helpers_mod.o SM_MaponiA3_mod.o QMCChem_dataset_test.o


## Default build target: build everything
all: cMaponiA3_test_3x3_3 fMaponiA3_test_3x3_3 fMaponiA3_test_4x4_2 QMCChem_dataset_test tests/test


## Compile recipes for C++
$(OBJ_DIR)/%.o: $(SRC_DIR)/%.cpp | $(OBJ_DIR)
	$(CXX) $(ARCH) $(CXXFLAGS) $(INCLUDE) -c -o $@ $<

## Compile recepies for Fortran
$(OBJ_DIR)/%.o: $(SRC_DIR)/%.f90 | $(OBJ_DIR)
	$(FC) $(ARCH) $(FFLAGS) -c -o $@ $<

## Explicit recipe to trigger rebuild and relinking when headerfile is changed
$(OBJ_DIR)/SM_MaponiA3.o: $(SRC_DIR)/SM_MaponiA3.cpp $(INC_DIR)/Helpers.hpp| $(OBJ_DIR)
	$(CXX) $(ARCH) $(CXXFLAGS) $(INCLUDE) -fPIC -c -o $@ $<


## Build tagets
.PHONY: all clean distclean

clean:
	@rm -vrf $(OBJ_DIR)

distclean: clean
	@rm -vrf $(BIN_DIR) \
	Slater* Updates.dat


## Linking the C++ example program
$(BIN_DIR)/cMaponiA3_test_3x3_3: $(OBJ_DIR)/$(cMaponiA3_test_3x3_3OBJ) $(OBJS) | $(BIN_DIR)
	$(CXX) $(ARCH) $(CXXFLAGS) $(INCLUDE) -o $@ $^

## Linking Fortran example program calling the C++ function
$(BIN_DIR)/fMaponiA3_test_3x3_3: $(OBJ_DIR)/$(fMaponiA3_test_3x3_3OBJ) $(OBJS) | $(BIN_DIR)
	$(FC) $(ARCH) $(FFLAGS) $(FLIBS) -o $@ $^

$(BIN_DIR)/fMaponiA3_test_4x4_2: $(OBJ_DIR)/$(fMaponiA3_test_4x4_2OBJ) $(OBJS) | $(BIN_DIR)
	$(FC) $(ARCH) $(FFLAGS) $(FLIBS) -o $@ $^

$(BIN_DIR)/QMCChem_dataset_test: $(OBJ_DIR)/$(QMCChem_dataset_testOBJ) $(OBJS) | $(BIN_DIR)
	$(FC) $(ARCH) $(FFLAGS) $(FLIBS) -o $@ $^

$(BIN_DIR)/test: $(SRC_DIR)/test.cpp $(OBJS) | $(BIN_DIR)
	$(H5CXX) $(H5CXXFLAGS) $(INCLUDE) -o $@ $^

$(BIN_DIR) $(OBJ_DIR):
	mkdir -p $@
