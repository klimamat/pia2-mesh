#include "init.h"

void initSod(Mesh *& m, Field<Compressible> *& W) {
	m = new Mesh(0.,1.,0.,0.05,20,1);
	W = new Field<Compressible>(*m); 
	
	
}