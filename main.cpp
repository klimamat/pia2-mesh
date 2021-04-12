#include "Mesh.h"
#include "Field.h"
#include "Compressible.h"
#include "init.h"
#include "BC.h"
#include "output.h"
#include <fenv.h>
#include <iostream>
#include <cstdlib>

int main(int iargc, char* iargv[]) {
    
    Mesh *m;
    Field<Compressible> *W;
    std::vector<BC<Compressible>*> boundary_conds; // Array of boundary conditions
    
    //feenableexcept(FE_DIVBYZERO || FE_INVALID || FE_OVERFLOW);
    
	//initSod(m,W,boundary_conds);
	initJet(m,W,boundary_conds);
	
	outputVTK("output_init.vtk",*m,*W);
	
	double dt, t = 0.0;
	const double t_max = 0.2;
	int n = 0;
	
	while (t < t_max) {
		dt = timestep(*m,*W);
		
		FVMstep(*m,*W,dt);
		
		for (auto bc : boundary_conds) bc->apply(*m,*W);
		
		n++;
		std::cout << "Step " << n << ", dt = " << dt << "\n";
		t += dt;
	}

	outputVTK("output.vtk",*m,*W);

	delete m; delete W;
	for (auto bc : boundary_conds) delete bc;
    	
	return 0;
}

