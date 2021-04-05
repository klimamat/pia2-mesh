#include "init.h"

void initSod(Mesh *& m, Field<Compressible> *& W, std::vector<BC<Compressible>*>& boundary_conds) {
	int nx = 1000;
	m = new Mesh(0.,1.,0.,1./double(nx),nx,1);
	W = new Field<Compressible>(*m); 
	boundary_conds.push_back(new SlipWallBC({1}));
			
	for (int i=0; i<m->nc; ++i) {
		Polygon const& T_x = m->cell[i];
		double rho, e, p;
		Vector2D u;
		if(T_x.centroid().x < 0.5){
			rho = 1.0;
			u = {0.0,0.0};
			p = 1.0;
		}
		else{
			rho = 0.125;
			u = {0.0,0.0};
			p = 0.1;
		}
			(*W)[i].rho = rho;
			(*W)[i].rhoU = rho*u;
			(*W)[i].e = p / (Compressible::kappa - 1.0) - 0.5 * rho * dot(u,u);
	}			
	
	for (auto bc : boundary_conds) bc->apply(*m,*W);				
}
