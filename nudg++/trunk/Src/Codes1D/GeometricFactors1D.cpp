// GeometricFactors1D.m
//---------------------------------------------------------
#include "NDGLib_headers.h"
#include "NDG1D.h"

// member version
//---------------------------------------------------------
void NDG1D::GeometricFactors1D()
//---------------------------------------------------------
{
    // Calculate geometric factors
    DMat xr = Dr*x;
    J = xr;
    rx = 1./J;
}
