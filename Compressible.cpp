#include "Compressible.h"
#include "Vector2D.h"

double Compressible::epsilon() const{
	Vector2D u=rhoU/rho;
	double uSq=dot(u,u);
	return e/rho+0.5*uSq;
}
double Compressible::p() const{
	double kappa = 1.4;
	return (kappa-1)*rho*epsilon();
}

Compressible fluxUpwind(Compressible Wl, Compressible Wr, Vector2D ne){
	Compressible F=Wl;
	
	Vector2D ul=Wl.rhoU / Wl.rho;
	Vector2D ur=Wr.rhoU / Wr.rho;
	
	Vector2D ue = 0.5*(ul+ur);
	
	if(dot(ue,ne)<0){
		F=Wr;
	}
	
	F.rho=F.rho*dot(ue,ne);
	F.rhoU=F.rhoU*dot(ue,ne);
	F.e=F.e*dot(ue,ne);
	
	double pe=0.5*(Wl.p()+Wr.p());
	
	F.rhoU = F.rhoU + pe*ne;
	F.e = F.e + pe * dot(ue,ne);
		
	return F;
}


double timestep(Mesh const& m, Field<Compressible> const& W) {
	
}

void FVMstep(Mesh const& m, Field<Compressible> & W, double dt) {
	
	Field<Compressible> res(m);
	
	for(int j=0;i < m.edge.size(); ++j){
		int l = m.edge[j].left();  // Index of the cell on the left
		int r = m.edge[j].right(); // Index of the cell on the right
		Compressible F = fluxUpwind(W[l],W[r],m.edge[j].normal());
		if (r != -1) {
			res[l] = res[l] + F;
			res[r] = res[r] - F;
		}
	}
	
	for(int j=0;j < m.cell.size(); ++j){
		W[j] = W[j]-(dt/m.cell[j].area())*res[j];
	}
}
