## Compilers
ARCH =
CXX = clang++-7
FC = flang-7
H5CXX = h5c++

## Compiler flags
CXXFLAGS = -O0 
FFLAGS = $(CXXFLAGS)
H5CXXFLAGS = $(CXXFLAGS) -fPIC
FLIBS = -lstdc++

INCLUDE = -I $(INC_DIR)/
DEPS_CXX = $(OBJ_DIR)/SM_MaponiA3.o $(OBJ_DIR)/SM_Standard.o
DEPS_F = $(DEPS_CXX) $(OBJ_DIR)/SM_MaponiA3_mod.o $(OBJ_DIR)/Helpers_mod.o

SRC_DIR := src
TST_DIR := tests
INC_DIR := include
OBJ_DIR := build
BIN_DIR := bin

EXEC := $(BIN_DIR)/cMaponiA3_test_3x3_3 \
		$(BIN_DIR)/test_internal_h5 \
		$(BIN_DIR)/fMaponiA3_test_3x3_3 \
		$(BIN_DIR)/fMaponiA3_test_4x4_2 \
		$(BIN_DIR)/QMCChem_dataset_test

## Build tagets
.PHONY: all clean distclean

all: $(EXEC)

clean:
	@rm -vrf $(OBJ_DIR)

distclean: clean
	@rm -vrf $(BIN_DIR) \
	Slater* Updates.dat


#### COMPILING
$(BIN_DIR) $(OBJ_DIR):
	mkdir -p $@

### IMPLICIT BUILD RULES
## C++ objects
$(OBJ_DIR)/%.o: $(TST_DIR)/%.cpp $(INC_DIR)/* | $(OBJ_DIR)
	$(CXX) $(CXXFLAGS) $(ARCH) $(INCLUDE) -c -o $@ $<

## HDF5/C++ objects
$(OBJ_DIR)/%_h5.o: $(TST_DIR)/%_h5.cpp $(INC_DIR)/* | $(OBJ_DIR)
	$(H5CXX) $(H5CXXFLAGS) $(INCLUDE) -c -o $@ $<

## Fortran modules
$(OBJ_DIR)/%_mod.o: $(SRC_DIR)/%_mod.f90 | $(OBJ_DIR)
	$(FC) $(FFLAGS) $(ARCH) -J$(OBJ_DIR)/ -c -o $@ $<

## Fortran objects
$(OBJ_DIR)/%.o: $(TST_DIR)/%.f90 | $(OBJ_DIR)
	$(FC) $(FFLAGS) $(ARCH) -I $(OBJ_DIR)/ -c -o $@ $<

### EXPLICIT BUILD RULES
## special compiler flag -fPIC otherwise h5c++ builds fail
$(OBJ_DIR)/SM_MaponiA3.o: $(SRC_DIR)/SM_MaponiA3.cpp $(INC_DIR)/* | $(OBJ_DIR)
	$(CXX) $(CXXFLAGS) -fPIC $(ARCH) $(INCLUDE) -c -o $@ $<

$(OBJ_DIR)/SM_Standard.o: $(SRC_DIR)/SM_Standard.cpp $(INC_DIR)/* | $(OBJ_DIR)
	$(CXX) $(CXXFLAGS) -fPIC $(ARCH) $(INCLUDE) -c -o $@ $<


#### LINKING
$(BIN_DIR)/cMaponiA3_test_3x3_3: $(OBJ_DIR)/cMaponiA3_test_3x3_3.o $(DEPS_CXX) | $(BIN_DIR)
	$(CXX) -o $@ $^

$(BIN_DIR)/test_internal_h5: $(OBJ_DIR)/test_internal_h5.o $(DEPS_CXX) | $(BIN_DIR)
	$(H5CXX) -o $@ $^

$(BIN_DIR)/fMaponiA3_test_3x3_3: $(DEPS_F) $(OBJ_DIR)/fMaponiA3_test_3x3_3.o  | $(BIN_DIR)
	$(FC) $(FLIBS) -o $@ $^

$(BIN_DIR)/fMaponiA3_test_4x4_2: $(DEPS_F) $(OBJ_DIR)/fMaponiA3_test_4x4_2.o  | $(BIN_DIR)
	$(FC) $(FLIBS) -o $@ $^

$(BIN_DIR)/QMCChem_dataset_test: $(DEPS_F) $(OBJ_DIR)/QMCChem_dataset_test.o  | $(BIN_DIR)
	$(FC) $(FLIBS) -o $@ $^
