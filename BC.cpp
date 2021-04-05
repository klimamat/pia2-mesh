#include "BC.h"

void SlipWallBC::apply(Mesh const& m, Field<Compressible> & W) {
	for (auto const& e : m.edge) {
		if (isCorrectLocation(e.location)) {
			int cl = e.left();
			int cr = e.right();
			Vector2D en = e.unitNormal();
			Vector2D un = en*dot(W[cl].u(),en); // Normal component of the fluid velocity relative to the edge
			Vector2D ur = W[cl].u() - 2.0*un; // Velocity vector mirrored along the edge
			
			W[cr].rho = W[cl].rho;
			W[cr].rhoU = W[cl].rho * ur;
			W[cr].e = W[cl].e;
		}
	}
}
