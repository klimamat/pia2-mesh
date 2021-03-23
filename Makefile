all: mesh

CPP = g++
CPPFLAGS = -std=c++11 #-DDEBUG

%.o: %.cpp
	$(CPP) $(CPPFLAGS) -c -o $@ $<

mesh: main.o Mesh.o output.o Vector2D.o Compressible.o
	$(CPP) $(CPPFLAGS) -o mesh main.o Mesh.o output.o Vector2D.o Compressible.o

clean:
	rm -f *.o mesh
