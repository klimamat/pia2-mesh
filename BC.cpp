#include "BC.h"
#include <iostream>

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
			W[cr].e = W[cr].eos_e_from_p(1);
		}
	}
}

void FreeBC::apply(Mesh const& m, Field<Compressible> & W) {
	for (auto const& e : m.edge) {
		if (isCorrectLocation(e.location)) {
			int cl = e.left();  // Internal cell index
			int cr = e.right(); // Ghost cell index
			
			W[cr].rho = W[cl].rho;
			W[cr].rhoU = W[cl].rhoU;
			W[cr].e = W[cl].e;
		}
	}
}

void puBC::apply(Mesh const& m, Field<Compressible> & W) {
	for (auto const& e : m.edge) {
		if (isCorrectLocation(e.location)) {
			int cl = e.left();  // Internal cell index
			int cr = e.right(); // Ghost cell index
			
			W[cr].rho = W[cl].rho;
			W[cr].rhoU = W[cl].rho*Vector2D(0.5,0.0);
			W[cr].e = W[cl].e;
		}
	}
}

void muBC::apply(Mesh const& m, Field<Compressible> & W) {
	for (auto const& e : m.edge) {
		if (isCorrectLocation(e.location)) {
			int cl = e.left();  // Internal cell index
			int cr = e.right(); // Ghost cell index
			
			W[cr].rho = W[cl].rho;
			W[cr].rhoU = W[cl].rho*Vector2D(-0.5,0.0);
			W[cr].e = W[cl].e;
		}
	}
}
