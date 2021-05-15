#include "BC.h"
#include <iostream>
#include <cmath>

void SlipWallBC::apply(Mesh const& m, Field<Compressible> & W) {
	for (auto const& e : m.edge) {
		if (isCorrectLocation(e.location)) {
			int cl = e.left();  // Internal cell index
			int cr = e.right(); // Ghost cell index
			
			Vector2D en = e.unitNormal();
			Vector2D un = en*dot(W[cl].u(),en); // Normal component of the fluid velocity relative to the edge
			Vector2D ur = W[cl].u() - 2.0*un; // Velocity vector mirrored along the edge
			
			W[cr].rho = W[cl].rho;
			W[cr].rhoU = W[cl].rho * ur;
			W[cr].e = W[cl].e;
		}
	}
}

void ReservoirBC::apply(Mesh const& m, Field<Compressible> & W) {
	for (auto const& e : m.edge) {
		if (isCorrectLocation(e.location)) {
			int cl = e.left();  // Internal cell index
			int cr = e.right(); // Ghost cell index
			
			W[cr].rho = 1.0;
			W[cr].rhoU = Vector2D(0.0,0.0);
			W[cr].e = W[cr].eos_e_from_p(1.0);
			
			
		}
	}
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////
void InletBC::apply(Mesh const& m, Field<Compressible> & W) {
	double pl, p_in, rho_in;
	double p0=101391.8555;
	double gama=1.4;

	
	for (auto const& e : m.edge) {
		if (isCorrectLocation(e.location)) {
			int cl = e.left();  // Internal cell index
			int cr = e.right(); // Ghost cell index
			
		//pl=p0;
		pl=W[cl].p();
			
		if (pl<p0){
			p_in=pl;
		}
		else{
		p_in=p0;
		}
			rho_in=pow((p_in/p0),(1/gama));
			double c_in=gama*(p_in/rho_in);
			double M_in=(2.0/(gama-1.0))*((pow((p_in/p0),-(gama-1.0)/gama))-1.0);
			double RhoU=rho_in*(sqrt(M_in*c_in));
			
			W[cr].rho = rho_in;
			W[cr].rhoU = Vector2D(RhoU,0.0);
			W[cr].e = W[cr].eos_e_from_p(p_in);
			
			
		}
	}
}

void OutletBC::apply(Mesh const& m, Field<Compressible> & W) {
	double pout;
	double p0=101391.8555;
	double Mout =0.675;
	double gama=1.4;
	
	for (auto const& e : m.edge) {
		if (isCorrectLocation(e.location)) {
			int cl = e.left();  // Internal cell index
			int cr = e.right(); // Ghost cell index
			
			pout=p0*(pow(1.0+(gama-1.0)*(pow(Mout,2.0)/2.0),(gama/(1.0-gama))));
			
			W[cr].e = W[cl].eos_e_from_p(pout);
			W[cr].rho = W[cl].rho;
			W[cr].rhoU = W[cl].rhoU;			
		}
	}
}