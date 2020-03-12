// Advec1D.cpp
// member routines for class Advec1D
// 2007/06/21
//---------------------------------------------------------
#include "NDGLib_headers.h"
#include "Advec1D.h"


//---------------------------------------------------------
Advec1D::Advec1D()
//---------------------------------------------------------
{
  class_name = "Advec1D";
  // set simulation parameters
}


//---------------------------------------------------------
Advec1D::~Advec1D()
//---------------------------------------------------------
{
}


//---------------------------------------------------------
void Advec1D::Resize()
//---------------------------------------------------------
{
  // Allocate member arrays
  u.resize   (Np, K); // solution fields
  rhsu.resize(Np, K); // Runge-Kutta stage values
}

//---------------------------------------------------------
void Advec1D::SetIC()
//---------------------------------------------------------
{
  // compute initial condition (time=0)
  u = sin(x);
}


//---------------------------------------------------------
void Advec1D::SetStepSize()
//---------------------------------------------------------
{
    // Set the current time step
}

//---------------------------------------------------------
void Advec1D::InitRun()
//---------------------------------------------------------
{
  K = 10;
  FinalTime = 10;
  // Generate 1D equi-spaced mesh with K+1 vertices
  SimpleMesh1D(0, 2., K);
  StartUp1D();      // construct grid and metric
  Resize();         // allocate arrays
  SetIC();          // set initial conditions
  SetStepSize();    // compute initial timestep (using IC's)

  //---------------------------------------------
  // base class version sets counters and flags
  //---------------------------------------------
  NDG1D::InitRun();

//Nreport =   1;      // set frequency of reporting
//Nreport =   5;      // set frequency of reporting
//Nreport =  20;      // set frequency of reporting
  Nreport =  50;      // set frequency of reporting
//Nreport = 100;      // no reports when timing
//Nreport = 500;      // no reports when timing

  Nrender = Nreport;  // output frequency
//Nrender = 100;      // output frequency
//Nrender = 10000;    // output frequency

//NvtkInterp = 16;  // set output resolution
//NvtkInterp = 10;  // set output resolution
  NvtkInterp =  8;  // set output resolution
//NvtkInterp =  6;  // set output resolution


  Summary();          // Show simulation details
}


//---------------------------------------------------------
void Advec1D::Summary()
//---------------------------------------------------------
{
  // TODO: add details of operators and sparse solvers
  NDG1D::Summary();
}


//---------------------------------------------------------
void Advec1D::Report(bool bForce)
//---------------------------------------------------------
{
  if (1 == tstep) {
    // print header
    umLOG(1, "\n step   time    rho(min) rho(max)   En(min)  En(max)       dt   \n"
               "----------------------------------------------------------------\n");
  }

  if (!umMOD(tstep,Nreport) || bForce || (1==tstep)) 
  {
    double r_min=u.min_col_val(1), r_max=u.max_col_val(1);

    umLOG(1, "%5d  %6.3lf   %8.5lf %8.5lf   %8.5lf %8.5lf   %8.6lf\n", 
              tstep, time, r_min, r_max, dt);
  }

  if (!umMOD(tstep,Nrender) || bForce)
  {
    // Output solution
  }
}


//---------------------------------------------------------
void Advec1D::FinalReport()
//---------------------------------------------------------
{
  // force report on final step
  this->Report(true);

  umLOG(1, "\n  time for NDG work:  %0.2lf secs\n",  time_work);
  umLOG(1,   "           rhs work:  %0.2lf\n",       time_rhs);
  umLOG(1,   " time for main loop:  %0.2lf secs\n\n",time_total);
}
