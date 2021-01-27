CC=icc
CXX=icpc
CFLAGS=-O0 -debug full
CXXFLAGS=-O0 -debug full -traceback
# ARCH=-xCORE-AVX2 

DEPS = SM-MaponiA3.cpp
OBJ = SM-MaponiA3.o 

%.o: %.c $(DEPS)
	$(CXX) $(ARCH) -c -o $@ $< $(CFLAGS)

SM-MaponiA3: $(OBJ)
	$(CXX) $(ARCH) -o $@ $^ $(CFLAGS)

clean:
	@rm -vf *.o
