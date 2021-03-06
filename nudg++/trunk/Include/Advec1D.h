// Advec1D.h
//---------------------------------------------------------
#ifndef NDG__Advec1D_H__INCLUDED
#define NDG__Advec1D_H__INCLUDED

#include "NDG1D.h"

//---------------------------------------------------------
class Advec1D : public NDG1D
//---------------------------------------------------------
{
public:
  Advec1D();
  virtual ~Advec1D();
  virtual void Driver();

protected:

  virtual void Run();

  virtual void Resize();        // resize system arrays
  virtual void SetIC();
  virtual void SetStepSize();
  virtual void InitRun();

  virtual void Summary();
  virtual void Report(bool bForce=false);
  virtual void FinalReport();

  virtual void RHS(DMat& u, double time, double a);

protected:

  //-------------------------------------
  // member data
  //-------------------------------------

  DMat u, rhsu;
  DMat resid;
};

#endif  // NDG__Advec1D_H__INCLUDED
