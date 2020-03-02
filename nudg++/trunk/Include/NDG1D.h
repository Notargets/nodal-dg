// NDG1D.h
// interface implemented by 1D solvers
//---------------------------------------------------------
#ifndef NDG__NDG_1D_H__INCLUDED
#define NDG__NDG_1D_H__INCLUDED

#include "Globals1D.h"
#include "Stopwatch.h"

// sparse matrix
#include "CS_Type.h"

//---------------------------------------------------------
class NDG1D : public Globals1D
//---------------------------------------------------------
{
public:
  NDG1D();
  virtual ~NDG1D();
  virtual void Driver()=0;  // All simulators implement Driver()

  const char* GetClassName() const      { return class_name.c_str(); }

protected:
  virtual void Run()=0;     // All simulators implement Run()
  virtual void InitRun();
  virtual void Summary();
  virtual void Report(bool bForce=false);
  virtual void FinalReport();
  virtual double GetAnalyticError();


  // Setup routines
  bool    StartUp1D();
  DMat&   Lift1D();
  void    Normals1D();
  void    BuildMaps1D();
  void    BuildBCMaps1D();
  void    BuildPeriodicMaps1D(double xperiod, double yperiod);

  void    Dmatrices1D();

  void    GeometricFactors1D();

  void    dtscale1D(DVec& dtscale);

  //-------------------------------------
  // Filters
  //-------------------------------------
  DMat&   Filter1D(int Norder, int Nc, double sp);
  DMat&   CutOffFilter1D(int Nc, double frac);
  void    filter_Q(const DMat& filter, DMat& Q);


  //-------------------------------------
  // Output functions
  //-------------------------------------
  // render selected solution fields
  void OutputVTK(const DMat& FData, int order, int zfield=0);

  // Utility routines
  void OutputSampleXYZ(int sample_N, DMat &newX, DMat &newY, DMat &newZ, 
                       const DMat &FData, DMat &newFData, int zfield=0);
  void OutputSampleELMT1D(int sample_N, IMat& ELMT);
  void OutputNodes(bool bFaceNodes=false);
  void OutputNodes_gauss();

protected:

  //-------------------------------------
  // member data
  //-------------------------------------
  stopwatch timer;          // timer class
  double  m_eps;            // machine precision
  double  ti0,ti1,tw1,trhs; // time data
  double  time_rhs,         // time for evaluating RHS
          time_rhs_c,       // time to adjust for curved elements
          time_flux,        // time to evaluate flux
          time_upw,         // time to evaluate eigen-flux
          time_bc,          // time to evaluate BCs
          time_source,      // time to evaluate source terms
          time_iter,        // time for 1 iteration
          time_limit,       // time used by "limiters"
          time_work,        // total NDG++ time
          time_total;       // total simulation time

  int     Nreport;          // frequency of reporting
  int     Nrender;          // frequency of rendering
  int     Nplotfield;       // selected field to plot
  int     NvtkInterp;       // Vtk interpolation order
  DMat    Q_plot;           // storage for fields to render
  DVec    dtscale;          // stepsize calculation
  DMat    m_Filter;         // filter

  // error analysis
  double  m_maxAbsError;
  DVec    m_ErrAnalytic;
  DVec    m_ErrEstimate;

  //-------------------------------------
  // static member data
  //-------------------------------------
  static int    N1Dobjects;     // count of active simulators 
  static int    PlotNumber;     // accumulated plot count
  static double TotalSimTime;   // accumulate time for multiple runs

  string  class_name;       // identify without RTTI
};

extern NDG1D* g_D1;        // global pointer

#endif  // NDG__NDG_1D_H__INCLUDED
