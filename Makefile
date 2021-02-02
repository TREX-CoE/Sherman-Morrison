CC=icc
CXX=icpc
CFLAGS=-O0 -debug full
CXXFLAGS=-O0 -debug full -traceback
# ARCH=-xCORE-AVX2 

DEPS = SM-MaponiA3.cpp
OBJ = SM-MaponiA3.o main.o

%.o: %.cpp $(DEPS)
	$(CXX) $(ARCH) -c -o $@ $< $(CFLAGS)

Sherman-Morrison: $(OBJ)
	$(CXX) $(ARCH) -o $@ $^ $(CFLAGS)

clean:
	@rm -vf *.o
