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
    VX.print();
    // Calculate face vertex locations
    EToV.resize(2, K);
    for (int i=0; i<K; i++) {
        EToV(1, i+1) = i;
        EToV(2, i+1) = i+1;
    }
    EToV.print();
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
}
