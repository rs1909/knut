// ------------------------------------------------------------------------- //
//
// This is part of KNUT
// Copyright (c) 2009 by Robert Szalai
//
// For license, see the file COPYING in the package's root directory
//
// ------------------------------------------------------------------------- //

#ifndef NCOLLOC_H
#define NCOLLOC_H

#include "basecolloc.h"
#include "system.h"
#include "matrix.h"
#include "spmatrix.h"

class ODEColloc : public PerSolColloc
{
  public:

    ODEColloc(System& _sys, const int nint, const int ndeg);        // computes mesh, metric, metricD

    ~ODEColloc()
    {}

    void init(const Vector& sol, const Vector& par);   // computes time, kk, ee, dd, rr ...

    void interpolate(Array3D<double>& out, const Vector& sol);
    void interpolateComplex(Array3D<double>& outRe, Array3D<double>& outIm, const Vector& sol);
    void interpolateOnMesh(Array3D<double>& out, const Vector& sol);

    // continuation of solution

    void rightHandSide(Vector& rhs, const Vector& par, const Vector& sol);
    void rightHandSide_p(Vector& rhs, const Vector& par, const Vector& sol, int p);   // sol is currently not needed
    void rightHandSide_x(SpMatrix& A, const Vector& par, const Vector& sol);          // sol is currently not needed
    void jacobianOfStability(SpMatrix& A, const Vector& par);

  private:
    
    template <bool periodic> void RHS_jacobian(SpMatrix& A, const Vector& par);

    // rINT and rDEG is included in idx. Only rDIM is necessary
    inline int& WRIDX(SpMatrix& A, int idx, int rDIM, int cDEG, int cDIM)
    {
      return A.WrLi(ndim + ndim*idx + rDIM, cDIM + ndim*cDEG);
    }

    inline double& WRDAT(SpMatrix& A, int idx, int rDIM, int cDEG, int cDIM)
    {
      return A.WrLx(ndim + ndim*idx + rDIM, cDIM + ndim*cDEG);
    }

    // for CharJac_x
    inline int& WRIDXCPLX(SpMatrix& A, int idx, int rDIM, int rIM, int cDEG, int cDIM, int cIM)
    {
      return A.WrLi(rIM + 2*(ndim + ndim*idx + rDIM), cIM + 2*(cDIM + ndim*cDEG));
    }

    inline double& WRDATCPLX(SpMatrix& A, int idx, int rDIM, int rIM, int cDEG, int cDIM, int cIM)
    {
      return A.WrLx(rIM + 2*(ndim + ndim*idx + rDIM), cIM + 2*(cDIM + ndim*cDEG));
    }


    // it stores all the collocation data
    Array3D<double> tt;       // interpolation at the collocation points
    Array3D<double> ttMSH;    // interpolation at the representation points

    // for rhs, and derivatives

    Array3D<double>  solData;
    Array2D<double>  p_fx;
    Array3D<double>  p_dfx;
    Array3D<double>  p_dfp;
    Array2D<double>& p_fxRe;
    Array2D<double>  p_fxIm;
    Array3D<double>& p_dfxRe;
    Array3D<double>  p_dfxIm;
    Array2D<double>  p_tauMSH;
    Array2D<double>  p_dtauMSH;
    Array2D<double>  p_fxMSH;
    Array3D<double>  p_dfxMSH;
    Array3D<double>  p_dfpMSH;
    Array3D<double>  p_dummy;
};

#endif
