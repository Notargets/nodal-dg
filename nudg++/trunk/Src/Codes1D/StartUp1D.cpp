// StartUp1D.m
// 
// 2007/06/06
//---------------------------------------------------------
#include "NDGLib_headers.h"
#include "NDG1D.h"


void NDG1D::SimpleMesh1D(double xmin, double xmax, int K)
{
    VX.resize(K+1);
    for (int i = 0; i<K+1; i++) {
        double val = (xmax - xmin)*double(i)/double(K) +xmin;
        VX(i+1) = val;
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
    invV = inv(V);
    Vr = GradVandermonde1D(N, r);
    Dr = Vr/V;

    LIFT = Lift1D(); // Compute surface lift terms
    Normals1D();

    IVec va = EToV(All, 1);
    IVec vb = EToV(All, 2);
    IVec vc = IVec(K);
    vc.range(1,K);
    x = ones(Np)*VX(va) + 0.5*(r+1.)*((VX(vb)-VX(va))(vc));

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
}
/*
% Compute basic Legendre Gauss Lobatto grid
r = JacobiGL(0,0,N);

% Build reference element matrices
V  = Vandermonde1D(N, r); invV = inv(V);
Dr = Dmatrix1D(N, r, V);

% Create surface integral terms
LIFT = Lift1D();

% build coordinates of all the nodes
va = EToV(:,1)'; vb = EToV(:,2)';
x = ones(N+1,1)*VX(va) + 0.5*(r+1)*(VX(vb)-VX(va));

% calculate geometric factors
[rx,J] = GeometricFactors1D(x,Dr);

% Compute masks for edge nodes
fmask1 = find( abs(r+1) < NODETOL)';
fmask2 = find( abs(r-1) < NODETOL)';
Fmask  = [fmask1;fmask2]';
Fx = x(Fmask(:), :);

% Build surface normals and inverse metric at surface
[nx] = Normals1D();
Fscale = 1./(J(Fmask,:));

% Build connectivity matrix
[EToE, EToF] = Connect1D(EToV);

% Build connectivity maps
[vmapM, vmapP, vmapB, mapB] = BuildMaps1D;
 */

