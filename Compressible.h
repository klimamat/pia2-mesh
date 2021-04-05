#ifndef COMPRESSIBLE_H
#define COMPRESSIBLE_H

#include "Mesh.h"
#include "Field.h"
#include "Vector2D.h"

class Compressible {
public:
    Compressible() : rho(0.0), rhoU(0.0,0.0), e(0.0) {};
    Compressible(double _rho, double _rhoux, double _rhouy, double _e) : rho(_rho), rhoU(_rhoux,_rhouy), e(_e) {};
    double rho, e;
    Vector2D rhoU;
    double epsilon() const;
    double p() const;
};

Compressible fluxUpwind(Compressible Wl, Compressible Wr, Vector2D ne);
double timestep(Mesh const& m, Field<Compressible> const& W);
void FVMstep(Mesh const& m, Field<Compressible> & W, double dt);

#endif // COMPRESSIBLE_H
