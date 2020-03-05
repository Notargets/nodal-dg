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


  // http://www.codeproject.com/cpp/FastDelegate.asp

  // define pointers to functions
  typedef void (Advec1D::*fp_IC)(const DVec& xi, const DVec& yi, double ti, DMat& Qo);
  typedef void (Advec1D::*fp_ES)(const DVec& xi, const DVec& yi, double ti, DMat& Qo);
  typedef void (Advec1D::*fp_BC)(const DVec& xi, const DVec& yi, const DVec& nxi, const DVec& nyi, const IVec& tmapI, const IVec& tmapO, const IVec& tmapW, const IVec& tmapC, double ti, DMat& Qio);

  fp_IC InitialSolution;  // function pointer to initial solution
  fp_ES ExactSolution;    // function pointer to exact solution
  fp_BC BCSolution;       // function pointer to boundary conditions

  void MapGaussFaceData();
  void PreCalcBdryData();
  void RHS(DMat& Qin, double ti, fp_BC SolutionBC);

  void Fluxes(DMat& Qin, DMat& F, DMat& G);
  void Fluxes(DMat& Qin, double gamma, DMat& F, DMat& G, DVec& rho, DVec& u, DVec& v, DVec& p);

protected:

  //-------------------------------------
  // member data
  //-------------------------------------

  DMat Q, rhsQ, resQ;
  DVec resid;
};

#endif  // NDG__Advec1D_H__INCLUDED
