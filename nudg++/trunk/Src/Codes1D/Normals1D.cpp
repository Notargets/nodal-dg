// Normals1D.m
//---------------------------------------------------------
#include "NDGLib_headers.h"
#include "NDG1D.h"


//---------------------------------------------------------
void NDG1D::Normals1D()
//---------------------------------------------------------
{
  // Purpose : Compute outward pointing normals at
  //	    elements faces as well as surface Jacobians
  nx.resize(Nfp*Nfaces, K);
  nx(1, All) = -1.0;
  nx(2, All) = 1.0;
}
