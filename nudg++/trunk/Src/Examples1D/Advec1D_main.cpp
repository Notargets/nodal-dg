// NDG.cpp: entry point for the NDG (console version)
// Note: reduced version (Curved Euler2D only)
// 2007/05/26
//---------------------------------------------------------
#include "NDGLib_headers.h"
#include "NDG_headers.h"
#include "Advec1D.h" 

//---------------------------------------------------------
int main(int argc, char* argv[])
//---------------------------------------------------------
{
  InitGlobalInfo();     // create global data and open logs

  umLOG(1, "\n");
  umLOG(1, "--------------------------------\n");
  umLOG(1, "              NuDG++            \n");
  umLOG(1, "  Nodal Discontinuous Galerkin  \n");
  umLOG(1, "     Method for non-linear      \n");
  umLOG(1, "          PDE systems           \n");
  umLOG(1, "--------------------------------\n\n");

  NDG1D *p = new Advec1D;

  if (p)
    {
        p->Driver();    // call driver
      delete p;       // delete simulator

      umLOG(1, "\nSimulation complete.\n\n");
    } else { 
      umWARNING("NDGDriver", "No simulator created"); 
    }

  FreeGlobalInfo();     // release global data and close logs
  return 0;
}
