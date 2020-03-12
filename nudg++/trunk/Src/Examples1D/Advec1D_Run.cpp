// Advec1D.m
// function Q = Advec1D(Q, FinalTime, ExactSolution, ExactSolutionBC, fluxtype)
// 2007/06/21
//---------------------------------------------------------
#include "NDGLib_headers.h"
#include "Advec1D.h"


//---------------------------------------------------------
void Advec1D::Run()
//---------------------------------------------------------
{
  // function Q = Advec1D(Q, FinalTime, ExactSolution, ExactSolutionBC, fluxtype)
    // Purpose  : Integrate 2D Euler equations using a 4th order low storage RK
    // compute time step size
    double xmin, CFL;
    // advection speed
    double a = 2*pi;
    double timelocal;

    InitRun();
    time = 0;

    // Runge-Kutta residual storage
    resid = zeros(Np,K);

    xmin = (abs(x(1,All)-x(2,All))).min_val();
    CFL=0.75; dt   = CFL/(2*pi)*xmin; dt = .5*dt;
    Nsteps = ceil(FinalTime/dt); dt = FinalTime/Nsteps;

    // outer time step loop
    for (int tstep=1; tstep<=Nsteps; tstep++) {
        for (int INTRK=1; INTRK<=5; INTRK++) {
            timelocal = time + rk4c(INTRK) * dt;
            rhsu = this->RHS(u, timelocal, a);
            resid = rk4a(INTRK) * resid + dt * rhsu;
            u = u + rk4b(INTRK) * resid;
        }
        time = time+dt;
        umLOG(1, "max_resid[%d] = %g\n", tstep, resid.max_val());
    }
    // Increment time
/*
  ti0=timer.read();       // time simulation loop

  // outer time step loop 
  while (time<FinalTime) {

    if (time+dt > FinalTime) {dt=FinalTime-time;}
    tw1=timer.read();   // time NDG work

    // 3rd order SSP Runge-Kutta
    this->RHS(Q, time,BCSolution);  Q1 = Q + dt*rhsQ;
    this->RHS(Q1,time,BCSolution);  Q2 = (3.0*Q + Q1 + dt*rhsQ)/4.0;
    this->RHS(Q2,time,BCSolution);  Q  = (Q + 2.0*Q2 + 2.0*dt*rhsQ)/3.0;

    time += dt;       // increment time 
    SetStepSize();    // compute new timestep

    time_work += timer.read() - tw1;  // accumulate cost of NDG work
    Report();         // optional reporting
    ++tstep;          // increment timestep

  }

  time_total = timer.read()-ti0;  // stop timing
  FinalReport();                  // final report
  */
}
