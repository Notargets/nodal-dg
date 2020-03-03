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
        VX[i] = val;
    };
    VX.print();
    // Calculate face vertex locations
    EToV.resize(K,2);
    for (int i=0; i<K; i++) {
        EToV[1][i+1] = i;
        EToV[2][i+1] = i+1;
    }
    EToV.print();
    exit(1);
}

//---------------------------------------------------------
bool NDG1D::StartUp1D()
//---------------------------------------------------------
{
  // Purpose : Setup script, building operators, grid, metric,
  //           and connectivity tables.

  // Definition of constants
  Nfp = 1; Np = N+1; Nfaces=2; NODETOL = 1e-12;
  K = 10;

  // Generate 1D equi-spaced mesh with K+1 vertices
  SimpleMesh1D(0, 2., K);

  // Compute basic Legendre Gauss Lobatto Grid
  r = JacobiGL(0,0,N);

  // Build reference element matrices
  V = Vandermonde1D(N,r); invV = inv(V);
  DVr = GradVandermonde1D(N, r);
  Dr = DVr/V;

  MassMatrix = trans(invV)*invV;

  /*
  // Create surface integral terms
  Lift1D();
  Normals1D();

  // build coordinates of all the nodes
  IVec va = EToV(All,1), vb = EToV(All,2), vc = EToV(All,3);

  // Note: outer products of (Vector,MappedRegion1D)
  x = 0.5 * (-(r+s)*VX(va) + (1.0+r)*VX(vb) + (1.0+s)*VX(vc));
  y = 0.5 * (-(r+s)*VY(va) + (1.0+r)*VY(vb) + (1.0+s)*VY(vc));

  // find all the nodes that lie on each edge
  IVec fmask1,fmask2,fmask3;
  fmask1 = find( abs(s+1.0), '<', NODETOL); 
  fmask2 = find( abs(r+s  ), '<', NODETOL);
  fmask3 = find( abs(r+1.0), '<', NODETOL);
  Fmask.resize(Nfp,3);                    // set shape (M,N) before concat()
  Fmask = concat(fmask1,fmask2,fmask3);   // load vector into shaped matrix

  Fx = x(Fmask, All); Fy = y(Fmask, All);


  // calculate geometric factors
  ::GeometricFactors1D(x,y,Dr,Ds,  rx,sx,ry,sy,J);

  // calculate geometric factors
  Fscale = sJ.dd(J(Fmask,All));


#if (0)
  OutputNodes(false); // volume nodes
  OutputNodes(true);  // face nodes
  umERROR("Exiting early", "Check {volume,face} nodes");
#endif

  // Build connectivity matrix
  tiConnect1D(EToV, EToE,EToF);

  // Build connectivity maps
  BuildMaps1D();

  // Compute weak operators (could be done in preprocessing to save time)
  DMat Vr,Vs;  GradVandermonde1D(N, r, s, Vr, Vs);
  VVT = V*trans(V);
  Drw = (V*trans(Vr))/VVT;  Dsw = (V*trans(Vs))/VVT;

  return true;
   */
}
