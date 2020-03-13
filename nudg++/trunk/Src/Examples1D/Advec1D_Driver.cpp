// Advec1DDriver.m
// Driver script for solving the 2D vacuum Euler's equations 
// 2007/06/21
//---------------------------------------------------------
#include "NDGLib_headers.h"
#include "Advec1D.h"


//---------------------------------------------------------
void Advec1D::Driver()
//---------------------------------------------------------
{
  //--------------------------------------------------
  // select order of polynomial approximation (N) 
  //--------------------------------------------------
//N = 3;
//N = 4;
//N = 6;
//  N = 8;
N = 2;
//N = 10;

  //--------------------------------------------------
  // select mesh, initial conditions, and BC function
  //--------------------------------------------------

  try {
    Run();  // Solve Problem
  } catch (...) {
    umWARNING("Advec1D::Driver", "Caught exception from Run()");
    return;
  }
}
