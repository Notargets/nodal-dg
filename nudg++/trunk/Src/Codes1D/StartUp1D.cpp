// StartUp1D.m
//---------------------------------------------------------
#include "NDGLib_headers.h"
#include "NDG1D.h"


void NDG1D::SimpleMesh1D(double xmin, double xmax, int K)
{
    int Nv = K+1;
    VX.resize(Nv);
    for (int i=1; i<=Nv; i++) {
        VX(i) = (xmax - xmin)*double(i-1)/double(K) +xmin;
    };
    // Calculate face vertex locations
    EToV.resize(K, 2);
    for (int i=1; i<=K; i++) {
        EToV(i, 1) = i;
        EToV(i, 2) = i+1;
    }
}

//---------------------------------------------------------
bool NDG1D::StartUp1D()
//---------------------------------------------------------
{
    // Purpose : Setup script, building operators, grid, metric,
    //           and connectivity tables.

    // Definition of constants
    Nfp = 1; Np = N+1; Nfaces=2; NODETOL = 1e-12;

    // Compute basic Legendre Gauss Lobatto Grid
    r = JacobiGL(0,0,N);

    // Build reference element matrices
    V = Vandermonde1D(N,r);
    GradVandermonde1D(N, r, Vr);
    Dr = Vr/V;

    LIFT = Lift1D(); // Compute surface lift terms
    Normals1D();

    IVec va = EToV(All, 1);
    IVec vb = EToV(All, 2);
    IVec vc;
    vc.range(1,K);
    DVec sT = VX(vb) - VX(va); // Allows for GC
    x = ones(Np)*(VX(va)) + 0.5*(r+1.)*(sT(vc));

    GeometricFactors1D();

    IVec fmask1 = find(abs(r+1.), '<', NODETOL);
    IVec fmask2 = find(abs(r-1.), '<', NODETOL);
    Fmask.reshape(1, 2);
    Fmask = concat(fmask1, fmask2);
    Fx = x(Fmask, All);

    Normals1D();

    DMat JJ = J(Fmask, All);
    Fscale = 1./JJ;

    Connect1D(EToV, EToE, EToF);

    BuildMaps1D();

    return true;
}
