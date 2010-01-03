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

class KNOdeBvpCollocation : public KNAbstractBvpCollocation
{
  public:

    KNOdeBvpCollocation(KNSystem& _sys, const int nint, const int ndeg);        // computes mesh, metric, metricD

    ~KNOdeBvpCollocation()
    {}

    void init(const KNVector& sol, const KNVector& par);   // computes time, kk, ee, dd, rr ...

    void interpolate(KNArray3D<double>& out, const KNVector& sol);
    void interpolateComplex(KNArray3D<double>& outRe, KNArray3D<double>& outIm, const KNVector& sol);
    void interpolateOnMesh(KNArray3D<double>& out, const KNVector& sol);

    // continuation of solution

    void rightHandSide(KNVector& rhs, const KNVector& par, const KNVector& sol);
    void rightHandSide_p(KNVector& rhs, const KNVector& par, const KNVector& sol, int p);   // sol is currently not needed
    void rightHandSide_x(KNSparseMatrix& A, const KNVector& par, const KNVector& sol);          // sol is currently not needed
    void jacobianOfStability(KNSparseMatrix& A, const KNVector& par);

  private:
    
    template <bool periodic> void RHS_jacobian(KNSparseMatrix& A, const KNVector& par);

    // rINT and rDEG is included in idx. Only rDIM is necessary
    inline int& WRIDX(KNSparseMatrix& A, int idx, int rDIM, int cDEG, int cDIM)
    {
      return A.writeIndex(ndim + ndim*idx + rDIM, cDIM + ndim*cDEG);
    }

    inline double& WRDAT(KNSparseMatrix& A, int idx, int rDIM, int cDEG, int cDIM)
    {
      return A.writeData(ndim + ndim*idx + rDIM, cDIM + ndim*cDEG);
    }

    // for CharJac_x
    inline int& WRIDXCPLX(KNSparseMatrix& A, int idx, int rDIM, int rIM, int cDEG, int cDIM, int cIM)
    {
      return A.writeIndex(rIM + 2*(ndim + ndim*idx + rDIM), cIM + 2*(cDIM + ndim*cDEG));
    }

    inline double& WRDATCPLX(KNSparseMatrix& A, int idx, int rDIM, int rIM, int cDEG, int cDIM, int cIM)
    {
      return A.writeData(rIM + 2*(ndim + ndim*idx + rDIM), cIM + 2*(cDIM + ndim*cDEG));
    }


    // it stores all the collocation data
    KNArray3D<double> tt;       // interpolation at the collocation points
    KNArray3D<double> ttMSH;    // interpolation at the representation points

    // for rhs, and derivatives

    KNArray3D<double>  solData;
    KNArray2D<double>  p_fx;
    KNArray3D<double>  p_dfx;
    KNArray3D<double>  p_dfp;
    KNArray2D<double>& p_fxRe;
    KNArray2D<double>  p_fxIm;
    KNArray3D<double>& p_dfxRe;
    KNArray3D<double>  p_dfxIm;
    KNArray2D<double>  p_tauMSH;
    KNArray2D<double>  p_dtauMSH;
    KNArray2D<double>  p_fxMSH;
    KNArray3D<double>  p_dfxMSH;
    KNArray3D<double>  p_dfpMSH;
    KNArray3D<double>  p_dummy;
};

#endif
