#include "Compressible.h"
#include "Vector2D.h"
#include <iostream>

double Compressible::epsilon() const{
	Vector2D u=rhoU/rho;
	double uSq=dot(u,u);
	return e/rho+0.5*uSq;
}
double Compressible::p() const{
	return (kappa-1)*rho*epsilon();
}
double Compressible::c() const {
	return std::sqrt(kappa*p()/rho);
}
Vector2D Compressible::u() const {
	return rhoU/rho;
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
	const double cfl = 0.9;
	double dt = 1e12;
	
	for (int i=0;i<m.cell.size();++i) {
		Polygon const& p = m.cell[i];
		Vector2D u = W[i].u();
		double c = W[i].c();
		double lambda = 0;
		
		for (int j=0;j<p.node_id.size();++j) {
			int j2; 
			if (j == p.node_id.size() - 1) j2 = 0;
			else j2 = j + 1;
			Point const& n1 = m.node[p.node_id[j]];
			Point const& n2 = m.node[p.node_id[j2]];
			Vector2D e = Vector2D(n1,n2);
			lambda += std::fabs(dot(u,e.normal())) + c*e.norm();
		}
		
		double dt_i = cfl * p.area() / lambda;
		if (dt_i < dt) dt = dt_i;
	}
	
	return dt;
}

void FVMstep(Mesh const& m, Field<Compressible> & W, double dt) {
	
	Field<Compressible> res(m);
	
	for(int j=0;j< m.edge.size(); ++j){
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
