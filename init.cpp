#include "init.h"
#include "MeshGmsh.h"

void initSod(Mesh *& m, Field<Compressible> *& W, std::vector<BC<Compressible>*>& boundary_conds) {
	int nx = 1000;
	m = new Mesh(0.,1.,0.,1./double(nx),nx,1);
	//m = new MeshGmsh("square.msh");
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
			(*W)[i].e = (*W)[i].eos_e_from_p(p);
	}			
	
	for (auto bc : boundary_conds) bc->apply(*m,*W);				
}

void initJet(Mesh *& m, Field<Compressible> *& W, std::vector<BC<Compressible>*>& boundary_conds) {
	m = new MeshGmsh("jet.msh");
	W = new Field<Compressible>(*m); 
	boundary_conds.push_back(new SlipWallBC({1}));
	boundary_conds.push_back(new ReservoirBC({2}));
			
	for (int i=0; i<m->nc; ++i) {
		Polygon const& T_x = m->cell[i];
			
			(*W)[i].rho = 1.0;
			(*W)[i].rhoU = Vector2D(0.0,0.0);
			(*W)[i].e = (*W)[i].eos_e_from_p(0.1);
	}			
	
	for (auto bc : boundary_conds) bc->apply(*m,*W);				
}

void initKH(Mesh *& m, Field<Compressible> *& W, std::vector<BC<Compressible>*>& boundary_conds) {
	m = new MeshGmsh("KelvinHelmholtz_3x05.msh");
	W = new Field<Compressible>(*m); 
	boundary_conds.push_back(new SlipWallBC({9}));  //horni
	boundary_conds.push_back(new SlipWallBC({10})); //dolni
	boundary_conds.push_back(new puBC({11})); //levo dole
	boundary_conds.push_back(new muBC({12})); //pravo hore
	boundary_conds.push_back(new puBC({13})); //pravo dole
	boundary_conds.push_back(new muBC({14}));	//levo hore	


	for (int i=0; i<m->nc; ++i) {
		Polygon const& T = m->cell[i];
		double rho, e, p;
		Vector2D u;
		if(T.centroid().y < 0.25){
			rho = 1.0;
			u = {0.5,0.0};
			p = 2.5;
		}
		else{
			rho = 2.0;
			u = {-0.5,0.0};
			p = 2.5;
		}
			(*W)[i].rho = rho;
			(*W)[i].rhoU = rho*u;
			(*W)[i].e = (*W)[i].eos_e_from_p(p);
	}			
	
	for (auto bc : boundary_conds) bc->apply(*m,*W);				
}
