## Compilers, compiler flags & external libs
ifeq ($(ENV),INTEL)
	CXX = icpx
	FC = ifort
	ARCH = -xCORE-AVX2
	OPT = -O3
	DEBUG = -g -debug full
else ifeq ($(ENV),LLVM)
	CXX = clang++
	FC = flang
	ARCH = -march=native
	OPT = -O3
	DEBUG = -g
else ifeq ($(ENV),GNU)
	CXX = g++
	FC = gfortran
	# ARCH = -mavx
	ARCH =
	OPT = -O0
	DEBUG = -g
else
    $(error No valid compiler environment set in $$ENV. \
	First run: $$ source smvars.sh {intel | llvm | gnu})
endif
HDF5_CXX = $(CXX)
H5CXX = h5c++
FLIBS = -lstdc++
CXXFLAGS = $(OPT) $(ARCH) $(DEBUG) $(THRESHOLD) -fPIC

## MKL linker flags
ifeq ($(MKL),-DMKL)
	CXXFLAGS += $(MKL)
	H5LFLAGS = -L$(MKLROOT)/lib/intel64 -Wl,--no-as-needed -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lpthread -lm -ldl
	ifeq ($(ENV),INTEL)
		LFLAGS = -mkl=sequential # implicit
	else
		LFLAGS = $(H5LFLAGS)
	endif
endif
H5CXXFLAGS = $(CXXFLAGS)
FFLAGS = $(CXXFLAGS)

## Includes and dependencies
INCLUDE = -I $(INC_DIR)/
DEPS_CXX = $(OBJ_DIR)/SM_Maponi.o \
		   $(OBJ_DIR)/SM_Standard.o \
		   $(OBJ_DIR)/Woodbury.o \
		   $(OBJ_DIR)/SMWB.o \
		   $(OBJ_DIR)/Helpers.o
DEPS_F = $(DEPS_CXX) \
		 $(OBJ_DIR)/finterface_mod.o \
		 $(OBJ_DIR)/helpers_mod.o

## QMCkl includes and linking
QMCKL_INCLUDE = -I $(SMROOT)/qmckl/build/include
QMCKLLFLAGS = -L$(SMROOT)/qmckl/build/lib -lqmckl


## Directory structure
SRC_DIR := src
TST_DIR := tests
INC_DIR := include
OBJ_DIR := build
BIN_DIR := bin

EXEC := $(BIN_DIR)/cMaponiA3_test_3x3_3 \
		$(BIN_DIR)/test_h5 \
		$(BIN_DIR)/fnu_test_h5 \
		$(BIN_DIR)/fMaponiA3_test_3x3_3 \
		$(BIN_DIR)/fMaponiA3_test_4x4_2 \
		$(BIN_DIR)/QMCChem_dataset_test

## Build tagets
.PHONY: all clean distclean

all: $(EXEC)

clean:
	@rm -vrf $(OBJ_DIR) *.dbg *.cmdx *.cmod *.ilm *.stb

distclean: clean
	@rm -vrf $(BIN_DIR) \
	Slater* Updates.dat


#### COMPILING
$(BIN_DIR) $(OBJ_DIR):
	mkdir -p $@

### IMPLICIT BUILD RULES
## C++ objects
$(OBJ_DIR)/%.o: $(TST_DIR)/%.cpp $(INC_DIR)/* | $(OBJ_DIR)
	$(CXX) $(CXXFLAGS) $(INCLUDE) -c -o $@ $<

$(OBJ_DIR)/%.o: $(SRC_DIR)/%.cpp $(INC_DIR)/* | $(OBJ_DIR)
	$(CXX) $(CXXFLAGS) $(INCLUDE) -c -o $@ $<

## HDF5/C++ objects
$(OBJ_DIR)/%_h5.o: $(TST_DIR)/%_h5.cpp $(INC_DIR)/* | $(OBJ_DIR)
	$(H5CXX) $(H5CXXFLAGS) $(INCLUDE) $(QMCKL_INCLUDE) -c -o $@ $<

## Fortran modules
$(OBJ_DIR)/%_mod.o: $(SRC_DIR)/%_mod.f90 | $(OBJ_DIR)
ifeq ($(ENV),$(filter $(ENV),LLVM GNU))
	$(FC) $(FFLAGS) -J $(OBJ_DIR)/ -c -o $@ $<
else
	$(FC) $(FFLAGS) -module $(OBJ_DIR)/ -c -o $@ $<
endif

## Fortran objects
$(OBJ_DIR)/%.o: $(TST_DIR)/%.f90 | $(OBJ_DIR)
	$(FC) $(FFLAGS) -I $(OBJ_DIR)/ -c -o $@ $<


#### LINKING
$(BIN_DIR)/cMaponiA3_test_3x3_3: $(OBJ_DIR)/cMaponiA3_test_3x3_3.o $(DEPS_CXX) | $(BIN_DIR)
	$(CXX) $(LFLAGS) -o $@ $^

$(BIN_DIR)/test_h5: $(OBJ_DIR)/test_h5.o $(DEPS_CXX) | $(BIN_DIR)
	$(H5CXX) $(H5LFLAGS) -o $@ $^

$(BIN_DIR)/fnu_test_h5: $(OBJ_DIR)/fnu_test_h5.o $(DEPS_CXX) | $(BIN_DIR)
	$(H5CXX) $(H5LFLAGS) -o $@ $^

$(BIN_DIR)/qmckl_test_h5: $(OBJ_DIR)/qmckl_test_h5.o | $(BIN_DIR)
	$(H5CXX) $(H5LFLAGS) $(QMCKLLFLAGS) -o $@ $^

$(BIN_DIR)/fMaponiA3_test_3x3_3: $(DEPS_F) $(OBJ_DIR)/fMaponiA3_test_3x3_3.o  | $(BIN_DIR)
	$(FC) $(LFLAGS) $(FLIBS) -o $@ $^

$(BIN_DIR)/fMaponiA3_test_4x4_2: $(DEPS_F) $(OBJ_DIR)/fMaponiA3_test_4x4_2.o  | $(BIN_DIR)
	$(FC) $(LFLAGS) $(FLIBS) -o $@ $^

$(BIN_DIR)/QMCChem_dataset_test: $(DEPS_F) $(OBJ_DIR)/QMCChem_dataset_test.o  | $(BIN_DIR)
	$(FC) $(LFLAGS) $(FLIBS) -o $@ $^
