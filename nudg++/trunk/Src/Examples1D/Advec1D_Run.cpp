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
    double a = 2.*pi; // advection speed
    double xmin, CFL;
    double timelocal;

    InitRun();

    // compute time step size
    xmin = (abs(x(1,All)-x(2,All))).min_val();
    CFL=0.75;
    dt   = .5 * (CFL/(2*pi)*xmin);
    umLOG(1, "dt = %g\n", dt);
    Nsteps = ceil(FinalTime/dt);
    dt = FinalTime/Nsteps;

    // outer time step loop
    resid = zeros(Np,K); // Runge-Kutta residual storage
    time = 0;
    Nsteps = 1;
    for (int tstep=1; tstep<=Nsteps; tstep++) {
        for (int INTRK=1; INTRK<=5; INTRK++) {
            timelocal = time + rk4c(INTRK) * dt;
            this->RHS(u, timelocal, a);
            resid = rk4a(INTRK) * resid + dt * rhsu;
            u += rk4b(INTRK) * resid;
        }
        time = time+dt;
        umLOG(1, "max_resid[%d] = %g, time = %g, dt = %g\n", tstep, resid.max_val(), time, dt);
        this->Report(true);
    }
    resid.print();
    u.print();
}
