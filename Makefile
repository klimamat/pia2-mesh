all: mesh

CPP = g++
CPPFLAGS = -std=c++11 -g #-DDEBUG

%.o: %.cpp
	$(CPP) $(CPPFLAGS) -c -o $@ $<

mesh: main.o Mesh.o output.o Vector2D.o Compressible.o init.o
	$(CPP) $(CPPFLAGS) -o mesh main.o Mesh.o output.o Vector2D.o Compressible.o init.o

clean:
	rm -f *.o mesh
