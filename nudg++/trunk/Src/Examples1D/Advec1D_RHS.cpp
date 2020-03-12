// Advec1DRHS.m
// Advec1DRHS.m
//---------------------------------------------------------
#include "NDGLib_headers.h"
#include "Advec1D.h"


//---------------------------------------------------------
DMat& Advec1D::RHS(DMat& u, double time, double a)
//---------------------------------------------------------
{
    DMat* prhsu = new DMat(Np, K);
    DMat& rhsu = (*prhsu);
    double alpha = 1.;
    double uin;
    DMat du = zeros(Nfp*Nfaces, K);

    // Face fluxes
    du = (u(vmapM)-u(vmapP)).dm(a*nx-(1.-alpha)*abs(a*nx))/2.;

    // Boundaries
    // Inflow boundary
    uin = -sin(a*time);
    du(mapI) = (u(vmapI)-uin).dm(a*nx(mapI)-(1.-alpha)*abs(a*nx(mapI)))/2.;
    du(mapO) = 0.;

    rhsu = -a*rx.dm(Dr*u) + LIFT*(Fscale.dm(du));

    return rhsu;
}
