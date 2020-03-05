// GradVandermonde1D.m
// function [Dvr] = GradVandermonde1D(N,r)
//---------------------------------------------------------
#include "NDGLib_headers.h"


//---------------------------------------------------------
DMat& GradVandermonde1D
(
        int   N,      // [in]
  const DVec& r      // [in]
)
//---------------------------------------------------------
{
    DMat *pDVr = new DMat(r.size(), N+1);
    DMat& DVr = (*pDVr);
// function [DVr] = GradVandermonde1D(N,r)
  // Purpose : Initialize the gradient of the modal basis (i)
  //		at (r) at order N

  // Initialize matrices
  for (int i=0; i<N; i++) {
      DVr(All,i+1)= GradJacobiP(r, 0, 0, i);
  }
  return DVr;
}
