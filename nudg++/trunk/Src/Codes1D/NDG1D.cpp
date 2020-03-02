// NDG1D.cpp
// 
// 2007/07/03
//---------------------------------------------------------
#include "NDGLib_headers.h"
#include "NDG1D.h"
#include "VecSort_Type.h"

NDG1D* g_D2 = NULL;        // global pointer

int     NDG1D::N1Dobjects=0;      // static counter
int     NDG1D::PlotNumber=0;      // accumulate plot count
double  NDG1D::TotalSimTime=0.0;  // accumulate sim time

//---------------------------------------------------------
NDG1D::NDG1D()
//---------------------------------------------------------
: 
  //-------------------------------------
  // initialize member data
  //-------------------------------------
  ti0(0.0), ti1(0.0)

{
  ++N1Dobjects;       // track allocation of NDG1D objects
  timer.start();      // initialize the timer
  g_D2 = this;        // set global simulation pointer

  // reporting parameters
  Nreport     = 10;   // frequency of reporting
  Nrender     = 10;   // frequency of rendering
  Nplotfield  = 1;    // selected field to plot
  NvtkInterp  = 6;    // Vtk interpolation order

  // get machine precision for relative tol. tests
  m_eps = 1e-12;
  double tol1 = 1.0 + m_eps;
  while (tol1 > 1.0) {
    m_eps = 0.5*m_eps;
    tol1 = 1.0 + m_eps;
  }

  // only valid if this->HasAnalyticSol();
  m_maxAbsError = -9.9e9;

  class_name = "NDG1DSim";
  umTRC(3, "Created NDG1D object (no. %d) \n", N1Dobjects);
}


//---------------------------------------------------------
NDG1D::~NDG1D()
//---------------------------------------------------------
{
  timer.stop();     // finished timing
  --N1Dobjects;     // decrement count of 1D simulators

  if (g_D2 == this) {
    // invalidate global pointer
    g_D2 = NULL;
  }

  umTRC(3, "Deleted NDG1D object (%d remain).\n", N1Dobjects);
}


//---------------------------------------------------------
void NDG1D::InitRun()
//---------------------------------------------------------
{
  // perform any initializations that should 
  // be done between multiple runs:

  time    =  0.0;     // reset time
  RKtime  =  0.0;     // reset RKtime
  tstep   =  1;       // reset current step
//Nsteps  =  0;       // reset number of steps
//Nreport =  2;       // set frequency of reporting (param)
//Nreport =  5;       // set frequency of reporting (param)
  Nreport = 20;       // set frequency of reporting (param)
//Nreport = 1000;     // no reports when timing
  Nrender = Nreport;  // output frequency           (param)

  // reset time counters
  ti0 = ti1   = 0.0;
  time_rhs    = 0.0;
  time_rhs_c  = 0.0;
  time_flux   = 0.0;
  time_upw    = 0.0;
  time_bc     = 0.0;
  time_source = 0.0;
  time_iter   = 0.0;
  time_limit  = 0.0;
  time_work   = 0.0;
  time_total  = 0.0;
}


//---------------------------------------------------------
void NDG1D::Summary()
//---------------------------------------------------------
{
  //-------------------------------------
  // default summary of current simulator
  //-------------------------------------
  umLOG(1, "\nNuDG++ 1D simulation:\n"
             "  Model type  = %s\n"
             "  Order (N)   = %d\n"
             "  Np          = %d\n"
             "  K           = %d\n"
             "  FileName    = %s\n"
             "  report freq = %d\n"
             "  render freq = %d\n"
             "  plot field  = %d\n"
             "  Vtk interp  = %d\n",
             this->GetClassName(), this->N, this->Np,
             Nreport, Nrender, Nplotfield, NvtkInterp);

    umLOG(1, "  Finaltime   = %0.2g\n"
             "  time-step   = %0.5g  (inital dt)\n"
             "  num. steps  = %d\n\n\n", 
             this->FinalTime, this->dt, this->Nsteps);
}


//---------------------------------------------------------
void NDG1D::Report(bool bForce)
//---------------------------------------------------------
{}

//---------------------------------------------------------
void NDG1D::FinalReport()
//---------------------------------------------------------
{}


//---------------------------------------------------------
double NDG1D::GetAnalyticError()
//---------------------------------------------------------
{
  // only valid if this->HasAnalyticSol();
  return m_maxAbsError;
}
