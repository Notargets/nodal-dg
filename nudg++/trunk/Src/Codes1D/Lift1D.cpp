// Lift1D.m
//---------------------------------------------------------
#include "NDGLib_headers.h"
#include "NDG1D.h"

//---------------------------------------------------------
DMat& NDG1D::Lift1D()
//---------------------------------------------------------
{
  // function [LIFT] = Lift1D()
  // Purpose  : Compute surface to volume lift term for DG formulation
  DMat Emat(Np, Nfaces*Nfp);
  Emat(1,1) = 1.0;
  Emat(Np, 2) = 1.0;
  return V*(trans(V)*Emat);
}
