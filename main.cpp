#include "Mesh.h"
#include "Field.h"
#include "Compressible.h"
#include "init.h"
#include "output.h"
#include <iostream>
#include <cstdlib>

int main(int iargc, char* iargv[]) {
    
    Mesh *m;
    Field<Compressible> *W;
    
	initSod(m,W);
	
	double dt, t = 0.0;
	const double t_max = 0.1;
	int n = 0;
	
	while (t < t_max) {
		dt = timestep(*m,*W);
		
		FVMstep(*m,*W,dt);
		
		applyBC(*m,*W);
		
		n++;
		std::cout << "Step " << n << ", dt = " << dt << "\n";
		t += dt;
	}

	outputVTK("output.vtk",*m,*W);

	delete m; delete W;
    	
	return 0;
}

