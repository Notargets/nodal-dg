cmake_minimum_required(VERSION 3.15)
project(NDG)

set(CMAKE_CXX_STANDARD 14)

include_directories(../Include)

find_library(BLAS libBlasLapack.a
        PATHS /home/llonergan/src/github.com/notargets/nodal-dg/nudg++/trunk/Lib)

add_definitions(-DUNDERSCORE -fpermissive -Wno-write-strings)

add_executable(Advect1D
        Examples1D/Advec1D.cpp
        Examples1D/Advec1D_main.cpp
        Examples1D/Advec1D_Driver.cpp
        Examples1D/Advec1D_Run.cpp
        Examples1D/Advec1D_RHS.cpp
        Codes1D/BuildMaps1D.cpp
        Codes1D/Connect1D.cpp
        Codes1D/GradVandermonde1D.cpp
        Codes1D/Lift1D.cpp
        Codes1D/Normals1D.cpp
        Codes1D/GeometricFactors1D.cpp
        Codes1D/StartUp1D.cpp
        Codes1D/NDG1D.cpp
        Codes1D/GradJacobiP.cpp
        Codes1D/JacobiGL.cpp
        Codes1D/JacobiGQ.cpp
        Codes1D/JacobiP.cpp
        Codes1D/Vandermonde1D.cpp
        Codes1D/Globals1D.cpp
        ServiceRoutines/LOG.cpp
        ServiceRoutines/INIT.cpp
        ServiceRoutines/Global_funcs.cpp
        Arrays/Mat_COL.cpp
        Arrays/ArrayMacros.cpp)

target_link_libraries(Advect1D ${BLAS} gfortran quadmath)

add_executable(Euler2D
        Arrays/ArrayMacros.cpp
        Arrays/Mat_COL.cpp
        Arrays/Sort_Index.cpp
        Codes1D/GradJacobiP.cpp
        Codes1D/JacobiGL.cpp
        Codes1D/JacobiGQ.cpp
        Codes1D/JacobiP.cpp
        Codes1D/Vandermonde1D.cpp
        Codes2D/BuildBCMaps2D.cpp
        Codes2D/BuildCurvedOPS2D.cpp
        Codes2D/BuildMaps2D.cpp
        Codes2D/BuildPeriodicMaps2D.cpp
        Codes2D/ConformingHrefine2D.cpp
        Codes2D/Connect2D.cpp
        Codes2D/Cub2D.cpp
        Codes2D/Cubature2D.cpp
        Codes2D/CubatureData2D.cpp
        Codes2D/CubatureVolumeMesh2D.cpp
        Codes2D/Curl2D.cpp
        Codes2D/CurvedDGDiv2D.cpp
        Codes2D/CurvedDGGrad2D.cpp
        Codes2D/CurvedDGJump2D.cpp
        Codes2D/CurvedPoissonIPDG2D.cpp
        Codes2D/CurvedPoissonIPDGbc2D.cpp
        Codes2D/CutOffFilter2D.cpp
        Codes2D/Div2D.cpp
        Codes2D/Dmatrices2D.cpp
        Codes2D/dtscale2D.cpp
        Codes2D/Filter2D.cpp
        Codes2D/FindLocalCoords2D.cpp
        Codes2D/GaussFaceMesh2D.cpp
        Codes2D/GeometricFactors2D.cpp
        Codes2D/Globals2D.cpp
        Codes2D/Grad2D.cpp
        Codes2D/GradSimplex2DP.cpp
        Codes2D/GradVandermonde2D.cpp
        Codes2D/Hrefine2D.cpp
        Codes2D/InterpMatrix2D.cpp
        Codes2D/Lift2D.cpp
        Codes2D/MakeCylinder2D.cpp
        Codes2D/NDG2D.cpp
        Codes2D/NDG2D_Output.cpp
        Codes2D/Nodes2D.cpp
        Codes2D/Normals2D.cpp
        Codes2D/PhysDmatrices2D.cpp
        Codes2D/rstoab.cpp
        Codes2D/Sample2D.cpp
        Codes2D/Simplex2DP.cpp
        Codes2D/StartUp2D.cpp
        Codes2D/tiConnect2D.cpp
        Codes2D/Vandermonde2D.cpp
        Codes2D/Warpfactor.cpp
        Codes2D/xytors.cpp
        Examples2D/CurvedEuler2D/ChannelBC2D.cpp
        Examples2D/CurvedEuler2D/ChannelIC2D.cpp
        Examples2D/CurvedEuler2D/CouetteBC2D.cpp
        Examples2D/CurvedEuler2D/CouetteIC2D.cpp
        Examples2D/CurvedEuler2D/CurvedEuler2D.cpp
        Examples2D/CurvedEuler2D/CurvedEuler2D_Driver.cpp
        Examples2D/CurvedEuler2D/CurvedEuler2D_Fluxes.cpp
        Examples2D/CurvedEuler2D/CurvedEuler2D_main.cpp
        Examples2D/CurvedEuler2D/CurvedEuler2D_RHS.cpp
        Examples2D/CurvedEuler2D/CurvedEuler2D_Run.cpp
        Examples2D/CurvedEuler2D/EulerHLL2D.cpp
        Examples2D/CurvedEuler2D/EulerLF2D.cpp
        Examples2D/CurvedEuler2D/EulerRoe2D.cpp
        Examples2D/CurvedEuler2D/IsentropicVortexBC2D.cpp
        Examples2D/CurvedEuler2D/IsentropicVortexIC2D.cpp
        ServiceRoutines/Global_funcs.cpp
        ServiceRoutines/INIT.cpp
        ServiceRoutines/LOG.cpp
        ServiceRoutines/MeshReaderGambit2D.cpp
        ServiceRoutines/MeshReaderGambit3D.cpp
        ServiceRoutines/Tokenizer.cpp
        Sparse/CHOLMOD_solver.cpp
        Sparse/CS_Cholinc.cpp
        Sparse/CS_Solve.cpp
        Sparse/CS_Utils.cpp)

target_link_libraries(Euler2D ${BLAS} gfortran quadmath)


