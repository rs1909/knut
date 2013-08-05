// ------------------------------------------------------------------------- //
//
// This is part of KNUT
// Copyright (c) 2013 by Robert Szalai
//
// For license, see the file COPYING in the package's root directory
//
// ------------------------------------------------------------------------- //

#ifndef STCOLLOC_H
#define STCOLLOC_H

#include "matrix.h"
#include "basecolloc.h"

class KNExprSystem;

class KNSteadyStateJacobian : public KNAbstractCollocation
{
  public:
    KNSteadyStateJacobian(KNExprSystem& sys_);
    virtual ~KNSteadyStateJacobian() {}
    
    // init(...) is called before any attempt to solve to calculate the jacobians
    void   init(const KNVector& sol, const KNVector& par);
    // this provides adjoint vectors
    void   star(KNVector& out, const KNVector& sol) { out = sol; }
    // this defines the inner product on the space
    double integrate(const KNVector& v1, const KNVector& v2) { return v1*v2; }
    double integrateWithCp(const KNVector& v1, const KNVector& v2, const KNVector& v3)
    {
      double res = 0.0;
      for (size_t k=0; k < v1.size(); k++) res += v1(k)*(v2(k) - v3(k));
      return res;
    }
    // this adapts the mesh of the numerical scheme
    void   meshAdapt(KNVector& newprofile, const KNVector& profile, KNVector& newtangent, const KNVector& tangent) { }
    
    void rightHandSide(KNVector& rhs, const KNVector& par, const KNVector& sol);
    void rightHandSide_p(KNVector& rhs, const KNVector& par, const KNVector& sol, size_t p);
    void rightHandSide_x(KNSparseMatrix& A, const KNVector& par, const KNVector& sol);
    
    // make these abstract in the base class
    // supplementary
    inline size_t nDim() const
    {
      return ndim;
    }
    inline size_t nPar() const
    {
      return npar;
    }

  private:

    // the equations
    KNExprSystem* sys;

    const size_t ndim;
    const size_t npar;
    const size_t ntau;
    
    KNVector time;

    // for rhs, and derivatives
    KNArray3D<double>  solData;
    KNArray2D<double>  p_fx;
    KNArray3D<double>  p_dfx;
    KNArray3D<double>  p_dfp;
    KNArray2D<double>& p_fxRe;
    KNArray2D<double>  p_fxIm;
    KNArray3D<double>& p_dfxRe;
    KNArray3D<double>  p_dfxIm;
    KNArray3D<double>  p_dummy;

};

#endif
