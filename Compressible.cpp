#include "Compressible.h"
#include "Vector2D.h"


double Compressible::epsilon() const{
	return e/rho;
}


double Compressible::p() const{
	float kappa = 1.4;
	
	double p = (kappa - 1) * rho * (epsilon() - 0.5 * sqrt( pow(rhoU.x,2) + pow(rhoU.y,2) ));
		
	return p;
}