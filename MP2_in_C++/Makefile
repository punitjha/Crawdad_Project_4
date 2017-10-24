CC=g++-5
CFLAGS= -std=c++11 -g -O3 -larmadillo -llapack -lblas 
OBJ = mp2.o
DEPS = mass.h hartree.h mmult.h
all: mp2.exe

mp2.exe: $(OBJ)
	$(CC) -o mp2.exe $^ $(CFLAGS)

%.o: %.cpp $(DEP)
	$(CC) -c $< $(CFLAGS)

clean:
	rm *.o *.exe
