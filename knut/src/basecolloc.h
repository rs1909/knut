// ------------------------------------------------------------------------- //
//
// This is part of KNUT
// Copyright (c) 2009 by Robert Szalai
//
// For license, see the file COPYING in the package's root directory
//
// ------------------------------------------------------------------------- //

#ifndef BASECOLLOC_H
#define BASECOLLOC_H

#include "matrix.h"

class KNExprSystem;

class KNAbstractCollocation
{
  public:
    virtual ~KNAbstractCollocation() {}

    // init(...) is called before any attempt to solve to calculate the jacobians
    virtual void   init(const KNVector& sol, const KNVector& par) = 0;
    // this provides adjoint vectors
    virtual void   star(KNVector& out, const KNVector& sol) = 0;
    // this defines the inner product on the space
    virtual double integrate(const KNVector& v1, const KNVector& v2) = 0;
    // this adapts the mesh of the numerical scheme
    virtual void   meshAdapt(KNVector& newprofile, const KNVector& profile, KNVector& newtangent, const KNVector& tangent) = 0;
};

// It can contain mesh adaptation, importing exporting reading writing to file
// the data structures, mesh, meshINT, col are the same.
class KNAbstractBvpCollocation : public KNAbstractCollocation
{
  public:

    KNAbstractBvpCollocation(KNExprSystem& sys, const size_t nint, const size_t ndeg);

    virtual ~KNAbstractBvpCollocation() {}

    virtual void init(const KNVector& sol, const KNVector& par) = 0;
    void meshAdapt(KNVector& newprofile, const KNVector& profile, KNVector& newtangent, const KNVector& tangent);

    virtual void interpolate(KNArray3D<double>& out, const KNVector& sol) = 0;
    virtual void interpolateComplex(KNArray3D<double>& outRe, KNArray3D<double>& outIm, const KNVector& sol) = 0;
    virtual void interpolateOnMesh(KNArray3D<double>& out, const KNVector& sol) = 0;

    static void   getMetric(KNMatrix& mt, const KNVector& t);
    static void   getDiffMetric(KNMatrix& mt, const KNVector& t);
    static void   star(KNVector& out, const KNVector& in, const KNMatrix& mt, const KNVector& msh, size_t dim);
    static double integrate(const KNVector& v1, const KNVector& v2, const KNMatrix& mt, const KNVector& msh, size_t dim);
    static size_t meshlookup_left(const KNVector& mesh, double t);
    static size_t meshlookup_right(const KNVector& mesh, double t);

    void   star(KNVector& out, const KNVector& sol);
    double integrate(const KNVector& v1, const KNVector& v2);
    double integrateWithDerivative(const KNVector& v1, const KNVector& v2);
    double integrateWithCp(const KNVector& v1, const KNVector& v2, const KNVector& v3);

    void   phaseStar(KNVector& V1, const KNVector& V2);
    void   phaseRotationStar(KNVector& V1, const KNVector& V2, const KNArray1D<size_t>& Re, const KNArray1D<size_t>& Im);

    void   importProfile(KNVector& out, const KNVector& in, const KNVector& mesh, size_t deg_, bool adapt);
    void   exportProfile(KNVector& out, const KNVector& mshint, const KNVector& mshdeg, const KNVector& in);
    void   pdMeshConvert(KNVector& newprofile, KNVector& newtangent, const KNVector& oldprofile, const KNVector& oldtangent);
    
    void   resetProfile(KNExprSystem& sys, KNVector& profile);

    // computing the Jacobians, right-hand sides, characteristic matrices

    // continuation of solution

    virtual void rightHandSide(KNVector& rhs, const KNVector& par, const KNVector& sol) = 0;
    virtual void rightHandSide_p(KNVector& rhs, const KNVector& par, const KNVector& sol, size_t p) = 0;   // sol is currently not needed
    virtual void rightHandSide_x(KNSparseMatrix& A, const KNVector& par, const KNVector& sol) = 0;          // sol is currently not needed

    // supplementary
    inline size_t nDim() const
    {
      return ndim;
    }
    inline size_t nPar() const
    {
      return npar;
    }
    inline size_t nInt() const
    {
      return nint;
    }
    inline size_t nDeg() const
    {
      return ndeg;
    }

    inline const KNVector& getElem()
    {
      return meshINT;
    }
    void setMesh(const KNVector& msh)
    {
      P_ASSERT_X(msh.size() == nint + 1, "Error in KNAbstractBvpCollocation::setMesh : Bad dimensions.");
      mesh = msh;
    }
    inline const KNVector& getMesh()
    {
      return mesh;
    }
    inline double Profile(size_t i, size_t d)
    {
      return mesh(i) + meshINT(d)*(mesh(i + 1) - mesh(i));
    }

  protected:

    // helper functions
    void meshAdapt_internal( KNVector& newmesh, const KNVector& profile );

    // the equations
    KNExprSystem* sys;
    KNArray1D<double>  mass;

    const size_t ndim;
    const size_t npar;

    const size_t nint;
    const size_t ndeg;

    KNVector    time;
    KNArray1D<double> timeMSH;  // the representation points

    // matrix for integration
    KNMatrix metric;
    // integration with derivatives (for phase conditions)
    KNMatrix metricPhase;

    KNVector mesh;
    KNVector meshINT;
    KNVector col;
    KNArray1D< KNArray1D<double> > lgr;
    // internal use for the initialization
    KNVector out;
};

#endif /*BASECOLLOC_H*/
