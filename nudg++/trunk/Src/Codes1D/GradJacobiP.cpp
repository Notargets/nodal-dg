// GradJacobiP.m
// function [dP] = GradJacobiP(x, alpha, beta, N);
// 2007/06/06
//---------------------------------------------------------
#include "NDGLib_headers.h"

//---------------------------------------------------------
DVec& GradJacobiP(const DVec& x, double alpha, double beta, int N)
//---------------------------------------------------------
{
  // function [dP] = GradJacobiP(z, alpha, beta, N);
  // Purpose: Evaluate the derivative of the orthonormal Jacobi
  //	   polynomial of type (alpha,beta)>-1, at points x
  //          for order N and returns dP[1:length(xp))]

  DVec* dP = new DVec(x.size(), 0.0, OBJ_temp, "dP");
  if (0 == N) {
    dP->fill(0.0);
  } else {
    (*dP) = sqrt(N*(N+alpha+beta+1))*JacobiP(x,alpha+1,beta+1, N-1);
  }
  return (*dP);
}
