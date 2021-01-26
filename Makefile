CC=icc
CXX=icpc
CFLAGS=

DEPS = SM-MaponiA3.cpp
OBJ = SM-MaponiA3.o 

%.o: %.c $(DEPS)
	$(CXX) -c -o $@ $< $(CFLAGS)

SM-MaponiA3: $(OBJ)
	$(CXX) -o $@ $^ $(CFLAGS)
