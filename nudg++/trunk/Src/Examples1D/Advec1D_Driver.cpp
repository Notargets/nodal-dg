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
  N = 8;
//N = 10;

  //--------------------------------------------------
  // select simulation type
  //--------------------------------------------------
  sim_type = eIsentropicVortex;
//sim_type = eChannelFlow;
//sim_type = eCouetteFlow;

  //--------------------------------------------------
  // select flux type
  //--------------------------------------------------
//flux_type = FT_LaxF;
  flux_type = FT_Roe;
//flux_type = FT_HLL;

  //--------------------------------------------------
  // select mesh, initial conditions, and BC function
  //--------------------------------------------------

  switch (sim_type) {
  case eIsentropicVortex:
    FileName        = "Grid/Euler2D/vortexA04.neu";
    InitialSolution = &Advec1D::IsentropicVortexIC2D;
    ExactSolution   = &Advec1D::IsentropicVortexIC2D;
    BCSolution      = &Advec1D::IsentropicVortexBC2D;
    break;

  case eChannelFlow:
    FileName        = "Grid/Euler2D/Euler01.neu";
    InitialSolution = &Advec1D::ChannelIC2D;
    ExactSolution   = NULL;
    BCSolution      = &Advec1D::ChannelBC2D;
    break;

  case eCouetteFlow:
  //FileName        = "Grid/Euler2D/Couette_K082.neu";
    FileName        = "Grid/Euler2D/Couette_K242.neu";
  //FileName        = "Grid/Euler2D/Couette_K856.neu";
    InitialSolution = &Advec1D::CouetteIC2D;
    ExactSolution   = NULL;
    BCSolution      = &Advec1D::CouetteBC2D;
    break;

  default:
    umERROR("Advec1D::Driver()", "Simulation case unknown");
    return;
  }

  //--------------------------------------------------
  // read mesh: vertices, elements, materials, BC's
  //--------------------------------------------------
  if (!MeshReaderGambit2D(FileName)) {
    umWARNING("Advec1D::Driver", "Error loading mesh (file: %s)\nExiting.\n",FileName.c_str());
    return;
  }

  try {
  //FinalTime =  0.1;
    FinalTime = 10.0;
  //FinalTime =  2.0;
  //FinalTime =  1.0;
    Run();  // Solve Problem
  } catch (...) {
    umWARNING("Advec1D::Driver", "Caught exception from Run()");
    return;
  }
}
