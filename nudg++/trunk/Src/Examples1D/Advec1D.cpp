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
  gamma = 1.4;
  gm1   = gamma - 1.0;
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
  Q.resize   (Np*K,  4);    // solution fields
  rhsQ.resize(Np*K,  4);    // Runge-Kutta stage values
  resQ.resize(Np*K,  4);    // Runge-Kutta residual 
}

//---------------------------------------------------------
void Advec1D::SetIC()
//---------------------------------------------------------
{
  // compute initial condition (time=0)
}


//---------------------------------------------------------
void Advec1D::SetStepSize()
//---------------------------------------------------------
{
  // function dt = Euler2Ddt(Q, gamma)
  // compute time step for the compressible Euler equations
  DVec rho,rhou,rhov,Ener, u,v,p,c,squv,Fscale_2,w_speeds, q1,q2,q3,q4;

  // since "self-mapping" of arrays is illegal, 
  // e.g. rho = rho(vmapM), use temp wrappers
  q1.borrow(Np*K,Q.pCol(1));  q2.borrow(Np*K,Q.pCol(2));
  q3.borrow(Np*K,Q.pCol(3));  q4.borrow(Np*K,Q.pCol(4));

  rho=q1(vmapM); rhou=q2(vmapM); rhov=q3(vmapM); Ener=q4(vmapM);
  u = rhou.dd(rho); v = rhov.dd(rho);  squv=sqr(u)+sqr(v);

  p = gm1 * (Ener - rho.dm(squv)/2.0);
  c = sqrt(abs(gamma*p.dd(rho)));

  Fscale_2 = 0.5*Fscale;
  w_speeds=SQ(N+1)*Fscale_2.dm(sqrt(squv)+c);
  dt = 1.0/w_speeds.max_val();

  Nsteps = (int)ceil(FinalTime/dt);
  dt = FinalTime/(double)Nsteps;
}


//---------------------------------------------------------
void Advec1D::InitRun()
//---------------------------------------------------------
{
  StartUp1D();      // construct grid and metric


  Resize();         // allocate arrays
  SetIC();          // set initial conditions
  SetStepSize();    // compute initial timestep (using IC's)

  // storage for residual at each time-step
  resid.resize(Nsteps+1);

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
    double r_min=Q.min_col_val(1), r_max=Q.max_col_val(1);
    double e_min=Q.min_col_val(4), e_max=Q.max_col_val(4);

    umLOG(1, "%5d  %6.3lf   %8.5lf %8.5lf   %8.5lf %8.5lf   %8.6lf\n", 
              tstep, time, r_min, r_max, e_min, e_max, dt);
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
