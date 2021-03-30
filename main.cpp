#include "Mesh.h"
#include "Field.h"
#include "Compressible.h"
#include "init.h"
#include "output.h"
#include <iostream>
#include <cstdlib>

int main(int iargc, char* iargv[]) {

    int nx = 20;
    int ny = 20;
    
    Mesh *m;
    Field<Compressible> *W;
    
	initSod(m,W);
	
	double dt, t = 0.0;
	const double t_max = 1.0;
	int n = 0;
	
	while (t < t_max) {
		dt = timestep(*m,*W);
		
		FVMstep(*m,*W,dt);
		
		n++;
		std::cout << "Step " << n << ", dt = " << dt << "\n";
		t += dt;
	}
	
	outputVTK("output.vtk",*m,*W);

	delete m; delete W;
    	
	return 0;
}

