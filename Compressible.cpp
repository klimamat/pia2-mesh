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
	
	double pe=0.5*(Wl.p()+Wr.p());
	
	F.rhoU = pe*ne;
	F.e = pe * dot(ue,ne);
		
	return F;
}
