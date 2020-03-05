// Globals1D.h
// container for (1D) global data
//---------------------------------------------------------
#ifndef NDG__Globals_1D_H__INCLUDED
#define NDG__Globals_1D_H__INCLUDED

#include "Global_funcs.h"
#include "VecObj_Type.h"


//---------------------------------------------------------
class Globals1D
//---------------------------------------------------------
{
public:

  Globals1D();
  virtual ~Globals1D();

  void init();
  void reset();
  void clear();

  int N, Nfp, Np, K;
  double  NODETOL;
  int Nfaces;
  DVec    r, x, VX;
  DMat    Dr, LIFT;
  DMat    Fx, nx, Fscale;
  DMat    V, invV, Vr;
  IVec    vmapB, mapB, vmapM, vmapP;
  IVec    vmapI, vmapO, mapI, mapO;
  IMat    Fmask, EToE, EToF, EToV;
  DVec    rk4a, rk4b, rk4c;
 int     tstep, Nsteps;
 double  dt, time, FinalTime, RKtime, pi, eps;
/*
global N Nfp Np K
global Nfaces EToE EToF
global NODETOL
global r x  VX
global Dr LIFT
global nx Fx Fscale
global vmapM vmapP vmapB mapB Fmask
global vmapI vmapO mapI mapO
global rx J
global rk4a rk4b rk4c
global V invV
 */
};

#endif // NDG__Globals_1D_H__INCLUDED
