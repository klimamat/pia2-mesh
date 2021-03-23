#ifndef COMPRESSIBLE_H
#define COMPRESSIBLE_H

#include "Vector2D.h"

class Compressible {
public:
    Compressible(double _rho, double _rhoux, double _rhouy, double _e) : rho(_rho), rhoU(_rhoux,_rhouy), e(_e) {};
    double rho, e;
    Vector2D rhoU;
    double epsilon() const;
    double p() const;
};

Compressible fluxUpwind(Compressible Wl, Compressible Wr, Vector2D ne);

#endif // COMPRESSIBLE_H
