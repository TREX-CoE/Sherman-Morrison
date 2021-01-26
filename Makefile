CC=icc
CXX=icpc
CFLAGS=-xCORE-AVX2 -O0 -debug full
CXXFLAGS=-xCORE-AVX2 -O0 -debug full

DEPS = SM-MaponiA3.cpp
OBJ = SM-MaponiA3.o 

%.o: %.c $(DEPS)
	$(CXX) -c -o $@ $< $(CFLAGS)

SM-MaponiA3: $(OBJ)
	$(CXX) -o $@ $^ $(CFLAGS)

clean:
	@rm -vf *.o
