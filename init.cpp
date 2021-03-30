#include "init.h"

void initSod(Mesh *& m, Field<Compressible> *& W) {
	m = new Mesh(0.,1.,0.,0.05,20,1);
	W = new Field<Compressible>(*m); 
	
void initSod(Mesh *& m, Field<Compressible> *& W) {
	m = new Mesh(0.,1.,0.,0.05,20,1);
	W = new Field<Compressible>(*m); 
	
	for (int i=0; i<m.cell.size(); ++i) {
			if(T_x >= 0.5){
			p_pom = 1.0;
			W[i].rho = 1.0;
			W[i].rhoU = 0.0;
			W[i].e = p_pom / (kappa - 1.0) + 0.5 * W[i].rho * W[i].u() * W[i].u();
		}else{
			p_pom = 0.1;
			W[i].rho = 0.125;
			W[i].rhoU = 0.0;
			W[i].e = p_pom / (kappa - 1.0) + 0.5 * W[i].rho * W[i].u() * W[i].u();
		
	}
	
}	
	
}
