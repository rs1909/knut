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

class System;

class BaseColloc
{
  public:
    virtual void   Init(const Vector& sol, const Vector& par) = 0;
    virtual void   Star(Vector& out, const Vector& sol) = 0;
    virtual double Integrate(const Vector& v1, const Vector& v2) = 0;
    virtual void   meshAdapt(Vector& newprofile, const Vector& profile, Vector& newtangent, const Vector& tangent) = 0;
};

// It can contain mesh adaptation, importing exporting reading writing to file
// the data structures, mesh, meshINT, col are the same.
class PerSolColloc : public BaseColloc
{
  public:

    PerSolColloc(System& sys, const int nint, const int ndeg);

    virtual ~PerSolColloc() {}

    virtual void Init(const Vector& sol, const Vector& par) = 0;
    void meshAdapt(Vector& newprofile, const Vector& profile, Vector& newtangent, const Vector& tangent);

    virtual void Interpolate(Array3D<double>& out, const Vector& sol) = 0;
    virtual void InterpolateCPLX(Array3D<double>& outRe, Array3D<double>& outIm, const Vector& sol) = 0;
    virtual void InterpolateMSH(Array3D<double>& out, const Vector& sol) = 0;

    static void   getMetric(Matrix& mt, const Vector& t);
    static void   getDiffMetric(Matrix& mt, const Vector& t);
    static void   star(Vector& out, const Vector& in, const Matrix& mt, const Vector& msh, int dim);
    static double integrate(const Vector& v1, const Vector& v2, const Matrix& mt, const Vector& msh, int dim);
    static int meshlookup(const Vector& mesh, double t);

    void   Star(Vector& out, const Vector& sol);
    double Integrate(const Vector& v1, const Vector& v2);
    double IntegrateDerivative(const Vector& v1, const Vector& v2);
    double IntegrateCont(const Vector& v1, const Vector& v2, const Vector& v3);

    void   PhaseStar(Vector& V1, const Vector& V2);
    void   PhaseRotStar(Vector& V1, const Vector& V2, const Array1D<int>& Re, const Array1D<int>& Im);

    void   Import(Vector& out, const Vector& in, const Vector& mesh, int deg_, bool adapt);
    void   Export(Vector& out, const Vector& mshint, const Vector& mshdeg, const Vector& in);
    void   pdMeshConvert(Vector& newprofile, Vector& newtangent, const Vector& oldprofile, const Vector& oldtangent);

    // computing the Jacobians, right-hand sides, characteristic matrices

    // continuation of solution

    virtual void RHS(Vector& rhs, const Vector& par, const Vector& sol) = 0;
    virtual void RHS_p(Vector& rhs, const Vector& par, const Vector& sol, int p) = 0;   // sol is currently not needed
    virtual void RHS_x(SpMatrix& A, const Vector& par, const Vector& sol) = 0;          // sol is currently not needed

    // supplementary
    inline int Ndim() const
    {
      return ndim;
    }
    inline int Npar() const
    {
      return npar;
    }
    inline int Nint() const
    {
      return nint;
    }
    inline int Ndeg() const
    {
      return ndeg;
    }

    inline const Vector& getElem()
    {
      return meshINT;
    }
    void setMesh(const Vector& msh)
    {
      P_ASSERT_X(msh.size() == nint + 1, "Error in PerSolColloc::setMesh : Bad dimensions.");
      mesh = msh;
    }
    inline const Vector& getMesh()
    {
      return mesh;
    }
    inline double Profile(int i, int d)
    {
      return mesh(i) + meshINT(d)*(mesh(i + 1) - mesh(i));
    }

  protected:

    // helper functions
    void meshAdapt_internal( Vector& newmesh, const Vector& profile );

    // the equations
    System* sys;

    const int ndim;
    const int npar;

    const int nint;
    const int ndeg;

    Vector    time;
    Array1D<double> timeMSH;  // the representation points

    // matrix for integration
    Matrix metric;
    // integration with derivatives (for phase conditions)
    Matrix metricPhase;

    Vector mesh;
    Vector meshINT;
    Vector col;
    Array1D< Array1D<double> > lgr;
    // internal use for the initialization
    Vector out;
};

#endif /*BASECOLLOC_H*/
