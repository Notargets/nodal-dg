// GradVandermonde1D.m
// function [Dvr] = GradVandermonde1D(N,r)
//---------------------------------------------------------
#include "NDGLib_headers.h"


//---------------------------------------------------------
void GradVandermonde1D
(
        int   N,      // [in]
  const DVec& r,      // [in]
  DMat& Vr      // [out]
)
//---------------------------------------------------------
{
// function [DVr] = GradVandermonde1D(N,r)
  // Purpose : Initialize the gradient of the modal basis (i)
  //		at (r) at order N

  Vr.resize(r.size(), N+1);
  // Initialize matrices
  for (int i=0; i<N; i++) {
      Vr(All,i+1)= GradJacobiP(r, 0, 0, i);
  }
}
