#include <omp.h>
#include "Compressible.h"
#include "Vector2D.h"
#include <iostream>
#include <iostream>
#include <fstream>
#include <string>

double Compressible::g = 0.0;

double Compressible::epsilon() const{
	Vector2D u=rhoU/rho;
	double uSq=dot(u,u);
	return e/rho-0.5*uSq;
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
double Compressible::eos_e_from_p(double p) const {
	return p / (kappa - 1.0) + 0.5 * rho * dot(u(),u());
}

Compressible fluxUpwind(Compressible Wl, Compressible Wr, Vector2D ne){
	double g=0.1;
	Compressible F=Wl;

	Vector2D ul=Wl.rhoU / Wl.rho;
	Vector2D ur=Wr.rhoU / Wr.rho;

	Vector2D ue = 0.5*(ul+ur);

	if(dot(ue,ne)<0){
		F=Wr;
	}

	F = F * dot(ue,ne);

	double pe=0.5*(Wl.p()+Wr.p());

	F.rhoU = F.rhoU + pe*ne;
	F.e = F.e + pe * dot(ue,ne);

	//F.rhoU.y+=F.rho*g;

	return F;
}

Compressible fluxNS(Compressible Wl, Compressible Wr, Vector2D ne, Point rl, Point rr){

	double mu =5.0e-5;
	Compressible D;
	Vector2D tau;
	Vector2D taux;
	Vector2D tauy;
	double tauU;

	Vector2D ul=Wl.rhoU / Wl.rho;
	Vector2D ur=Wr.rhoU / Wr.rho;
	Vector2D ue = 0.5*(ul+ur);

	double rd = pow((rr.x-rl.x),2)+pow((rr.y-rl.y),2);
	double a = (double)2/3;
	double pe=0.5*(Wl.p()+Wr.p());

	double uxx=(ur.x-ul.x)*(rr.x-rl.x)/rd;  // dux/dx
	double uxy=(ur.x-ul.x)*(rr.y-rl.y)/rd;	// dux/dy
	double uyy=(ur.y-ul.y)*(rr.y-rl.y)/rd;	// duy/dy
	double uyx=(ur.y-ul.y)*(rr.x-rl.x)/rd;	// duy/dx

	 tau.x=mu*((2*uxx-a*(uxx+uyy))*ne.x+(uxy+uyx)*ne.y);
	 tau.y=mu*((uxy+uyx)*ne.x+(2*uyy-a*(uxx+uyy))*ne.y);
	 tauU=mu*( (2*uxx-a*(uxx+uyy))*ue.x*ne.x+(uxy+uyx)*ue.y*ne.x+(uxy+uyx)*ue.x*ne.y+(2*uyy-a*(uxx+uyy))*ue.y*ne.y );

	D.rhoU = tau;
	D.e = tauU;

	return D;
}


Compressible fluxHLL(Compressible Wl, Compressible Wr, Vector2D ne){
	Vector2D nu = ne/ne.norm();
	double ul = dot(Wl.u(),ne)/ne.norm(); double ur = dot(Wr.u(),ne)/ne.norm();

	double rhols = std::sqrt(Wl.rho); double rhors = std::sqrt(Wr.rho);
	double ub = (rhols*ul + rhors*ur)/(rhols + rhors);
	double Hl = (Wl.e + Wl.p())/Wl.rho; double Hr = (Wr.e + Wr.p())/Wr.rho;
	double Hb = (rhols*Hl + rhors*Hr)/(rhols + rhors);
	double ab = std::sqrt((Compressible::kappa - 1.0)*(Hb - 0.5*ub*ub));
	double Sl = ub - ab; double Sr = ub + ab;

	Vector2D vl = Wl.u(); Vector2D vr = Wr.u();
	Vector2D ve = 0.5*(vl+vr);
	double pe=0.5*(Wl.p()+Wr.p());


	Compressible F;
	if (Sl >= 0.0) {
		F = Wl * dot(vl,ne) + Compressible(0.0,Wl.p()*ne,Wl.p()*dot(vl,ne));
    	//std::cout << "Fl = " << F.rho << " " << F.rhoU.norm() << " " << F.e << "\n";
	}
	else if (Sr <= 0.0) {
		F = Wr * dot(vr,ne) + Compressible(0.0,Wr.p()*ne,Wr.p()*dot(vr,ne));
    	//std::cout << "Fr = " << F.rho << " " << F.rhoU.norm() << " " << F.e << "\n";
	}
	else {
		Compressible Fl = Wl * dot(vl,nu) + Compressible(0.0,Wl.p()*nu,Wl.p()*dot(vl,nu));
		Compressible Fr = Wr * dot(vr,nu) + Compressible(0.0,Wr.p()*nu,Wr.p()*dot(vr,nu));
		F = ne.norm()*(Sr*Fl - Sl*Fr + Sl*Sr*(Wr - Wl))/(Sr-Sl);
    	/*if (Fl.rhoU.x != 0.0) {
    		std::cout << "Sr = " << Sr << ", Sl = " << Sl << "\n";
    		std::cout << "Fl = " << Fl.rho << " " << Fl.rhoU.x << " " << Fl.rhoU.y << " " << Fl.e << "\n";
    		std::cout << "Fr = " << Fr.rho << " " << Fr.rhoU.x << " " << Fr.rhoU.y << " " << Fr.e << "\n";
    		std::cout << "Fhll = " << F.rho << " " << F.rhoU.x << " " << F.rhoU.y << " " << F.e << "\n";
        }*/
	}
	return F;
}

double timestep(Mesh const& m, Field<Compressible> const& W) {
	const double cfl = 0.3;
	double dt = 1e12;
	double dt_local_min;

	#pragma omp parallel private(dt_local_min) shared(dt)
	{
	dt_local_min = 1e12;

	#pragma omp for
	for (int i=0;i<m.nc;++i) {
		Polygon const& p = m.cell[i];
		Vector2D u = W[i].u(); // Fluid velocity
		double c = W[i].c();   // Sound speed
		double lambda = 0.0;

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
		if (dt_i < dt_local_min) dt_local_min = dt_i;
	}

	#pragma omp critical
	{
		if (dt_local_min < dt) dt = dt_local_min;
	}

	}

	return dt;
}

void FVMstep(Mesh const& m, Field<Compressible> & W, double dt) {

	Field<Compressible> res(m);

	double g = Compressible::g;

	#pragma omp parallel for schedule(dynamic)
	for (int i=0; i<m.edge.size(); ++i) {
		auto const& e = m.edge[i];
		int l = e.left();  // Index of the cell on the left
		int r = e.right(); // Index of the cell on the right
		//Compressible F = fluxUpwind(W[l],W[r],e.normal());
		Compressible F = fluxHLL(W[l],W[r],e.normal());

//		if (!e.boundary) {
//		F = F - fluxNS(W[l],W[r],e.normal(),m.cell[l].centroid(),m.cell[r].centroid());
//		}
		#pragma omp critical
		{
			res[l] = res[l] + F;
			if (!e.boundary) {
				res[r] = res[r] - F;
			}
		}
	}

	#pragma omp parallel for
	for(int j=0;j<m.nc; ++j){
//		double const rho = W[j].rho;
//		Compressible Fg = -Compressible::g*Compressible(0.0,0.0,W[j].rho,W[j].rhoU.y);
		W[j] = W[j]-(dt/m.cell[j].area())*res[j]; // + dt*Fg;
	}
}
