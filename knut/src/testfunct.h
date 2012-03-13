// ------------------------------------------------------------------------- //
//
// This is part of KNUT
// Copyright (c) 2005 by Robert Szalai
//
// For license, see the file COPYING in the package's root directory
//
// ------------------------------------------------------------------------- //

#ifndef TESTFUNCT_H
#define TESTFUNCT_H

#include <hypermatrix.h>

// forward declarations
class KNDdeBvpCollocation;

#define TF_NKERNITER 20
#define TF_KERNEPS 1.0e-12

class KNAbstractTestFunctional
{
  public:
    KNAbstractTestFunctional() { kernEps = TF_KERNEPS; kernIter = TF_NKERNITER; }
    virtual        ~KNAbstractTestFunctional()
    {}
    virtual void   init(KNDdeBvpCollocation& col, const KNVector& par, const KNVector& sol) = 0;
    virtual double initStep() = 0;
    virtual double funct(KNDdeBvpCollocation& col, const KNVector& par, const KNVector& sol) = 0;
    virtual double funct_p(KNDdeBvpCollocation& col, const KNVector& par, const KNVector& sol, size_t alpha) = 0;
    virtual void   funct_x(KNVector& func, KNDdeBvpCollocation& col, const KNVector& sol, const KNVector& par) = 0;

    virtual void init(KNDdeBvpCollocation& col, const KNVector& par, const KNVector& sol,
                      double Re, double Im) = 0;
    virtual void funct(double& f1, double& f2,
                       KNDdeBvpCollocation& col, const KNVector& par, const KNVector& sol, double Re, double Im) = 0;
    virtual void funct_p(double& f1, double& f2,
                         KNDdeBvpCollocation& col, const KNVector& par, const KNVector& sol,
                         size_t alpha) = 0;
    virtual void funct_z(double& f1, double& f2,
                         KNDdeBvpCollocation& col, const KNVector& par, const KNVector& sol) = 0;
    virtual void funct_x(KNVector& func1, KNVector& func2,
                         KNDdeBvpCollocation& col, const KNVector& par, const KNVector& sol) = 0;
    virtual void kernel(KNVector& phi) = 0;
    virtual void setKernelTolerance(double eps, size_t iter) { kernEps = eps; kernIter = iter; }
  protected:
    double kernEps;
    size_t kernIter;
};

class KNTestFunctional : public KNAbstractTestFunctional
{
  public:
    KNTestFunctional(KNDdeBvpCollocation& col, double Z);
    ~KNTestFunctional();
    void   init(KNDdeBvpCollocation& col, const KNVector& par, const KNVector& sol);
    double initStep();
    double funct(KNDdeBvpCollocation& col, const KNVector& par, const KNVector& sol);
    double funct_p(KNDdeBvpCollocation& col, const KNVector& par, const KNVector& sol, size_t alpha);
    void   funct_x(KNVector& func, KNDdeBvpCollocation& col, const KNVector& par, const KNVector& sol);

    void init(KNDdeBvpCollocation&, const KNVector&, const KNVector&, double, double)
    {}
    void funct(double&, double&, KNDdeBvpCollocation&, const KNVector&, const KNVector&, double, double)
    {}
    void funct_p(double&, double&, KNDdeBvpCollocation&, const KNVector&, const KNVector&, size_t)
    {}
    void funct_z(double&, double&, KNDdeBvpCollocation&, const KNVector&, const KNVector&)
    {}
    void funct_x(KNVector&, KNVector&, KNDdeBvpCollocation&, const KNVector&, const KNVector&)
    {}
    void kernel(KNVector& phi);

  private:
    bool        first;
    double      ZZ;
    KNSparseBlockMatrix AHAT;
    KNVector      A_p;
    KNSparseMatrix    A_x;
    KNVector      rhs;
    KNVector      uu;
    KNVector      vv;
    KNVector      uudiff;
    KNVector      vvdiff;
    KNArray3D<double> vvData;
};

class KNComplexTestFunctional : public KNAbstractTestFunctional
{
  public:
    KNComplexTestFunctional(KNDdeBvpCollocation& col);
    ~KNComplexTestFunctional();
    void   init(KNDdeBvpCollocation&, const KNVector&, const KNVector&)
    {}
    double funct(KNDdeBvpCollocation&, const KNVector&, const KNVector&)
    {
      return 0.0;
    }
    double funct_p(KNDdeBvpCollocation&, const KNVector&, const KNVector&, size_t)
    {
      return 0.0;
    }
    void   funct_x(KNVector&, KNDdeBvpCollocation&, const KNVector&, const KNVector&)
    {}

    void init(KNDdeBvpCollocation& col, const KNVector& par, const KNVector& sol,
              double Re, double Im);
    double initStep();
    void funct(double& f1, double& f2,
               KNDdeBvpCollocation& col, const KNVector& par, const KNVector& sol, double Re, double Im);
    void funct_p(double& f1, double& f2,
                 KNDdeBvpCollocation& col, const KNVector& par, const KNVector& sol,
                 size_t alpha);
    void funct_z(double& f1, double& f2,
                 KNDdeBvpCollocation& col, const KNVector& par, const KNVector& sol);
    void funct_x(KNVector& func1, KNVector& func2,
                 KNDdeBvpCollocation& col, const KNVector& par, const KNVector& sol);
    inline void kernel(KNVector&)
    {}
    void   kernel(KNVector& Re, KNVector& Im, double& alpha);
    double kernelComplex(double& newperiod, KNVector& Re, KNVector& Im, KNDdeBvpCollocation& col, const KNVector& par);

  private:
    bool        first;
    double      ZRe, ZIm;
    KNSparseBlockMatrix AHAT;
    KNVector      A_p;
    KNSparseMatrix    A_x;
    KNVector      rhs;
    KNVector      one;
    KNVector      uu;
    KNVector      vv;
    KNVector      uudiff;
    KNVector      vvdiff;
    KNVector      gg;
    KNVector      hh;
    KNVector      ggdiff;
    KNVector      hhdiff;
    KNArray3D<double> vvDataRe;
    KNArray3D<double> vvDataIm;
};

class KNLpAutTestFunctional : public KNAbstractTestFunctional
{
  public:
    KNLpAutTestFunctional(KNDdeBvpCollocation& col, double Z);
    ~KNLpAutTestFunctional();
    void   init(KNDdeBvpCollocation& col, const KNVector& par, const KNVector& sol);
    double initStep();
    double funct(KNDdeBvpCollocation& col, const KNVector& par, const KNVector& sol);
    double funct_p(KNDdeBvpCollocation& col, const KNVector& par, const KNVector& sol, size_t alpha);
    void   funct_x(KNVector& func, KNDdeBvpCollocation& col, const KNVector& sol, const KNVector& par);

    void init(KNDdeBvpCollocation&, const KNVector&, const KNVector&, double, double)
    {}
    void funct(double&, double&, KNDdeBvpCollocation&, const KNVector&, const KNVector&, double, double)
    {}
    void funct_p(double&, double&, KNDdeBvpCollocation&, const KNVector&, const KNVector&, size_t)
    {}
    void funct_z(double&, double&, KNDdeBvpCollocation&, const KNVector&, const KNVector&)
    {}
    void funct_x(KNVector&, KNVector&, KNDdeBvpCollocation&, const KNVector&, const KNVector&)
    {}
    void kernel(KNVector& phi);

  private:
    bool        first;
    double      ZZ;
    KNSparseBlockMatrix AHAT;
    KNVector      A_p;
    KNSparseMatrix    A_x;
    KNSparseMatrix    mB;
    KNVector      mB_p;
    KNSparseMatrix    mB_x;
    KNVector      rhs;   // this is zero all the time
    KNVector      phi;
    KNVector      temp;
    KNVector      DpPhi;
    KNVector      DxPhi;
    KNVector      uu2;   // for computing the generalized eigenvector
    KNVector      vv2;
    KNVector      gg2;
    KNVector      hh2;
    KNVector      one2;  // this is ( 0.0, 1.0 )
    KNArray3D<double> phiData;
    KNArray3D<double> vv2Data;
    KNArray3D<double> solMSHData;
};

class KNLpAutRotTestFunctional : public KNAbstractTestFunctional
{
  public:
    KNLpAutRotTestFunctional(KNDdeBvpCollocation& col, KNArray1D<size_t> CRe, KNArray1D<size_t> CIm, double Z);
    ~KNLpAutRotTestFunctional();
    void   init(KNDdeBvpCollocation& col, const KNVector& par, const KNVector& sol);
    double initStep();
    double funct(KNDdeBvpCollocation& col, const KNVector& par, const KNVector& sol);
    double funct_p(KNDdeBvpCollocation& col, const KNVector& par, const KNVector& sol, size_t alpha);
    void   funct_x(KNVector& func, KNDdeBvpCollocation& col, const KNVector& par, const KNVector& sol);

    void init(KNDdeBvpCollocation&, const KNVector&, const KNVector&, double, double)
    {}
    void funct(double&, double&, KNDdeBvpCollocation&, const KNVector&, const KNVector&, double, double)
    {}
    void funct_p(double&, double&, KNDdeBvpCollocation&, const KNVector&, const KNVector&, size_t)
    {}
    void funct_z(double&, double&, KNDdeBvpCollocation&, const KNVector&, const KNVector&)
    {}
    void funct_x(KNVector&, KNVector&, KNDdeBvpCollocation&, const KNVector&, const KNVector&)
    {}
    void kernel(KNVector& phi);

  private:
    bool         first;
    double       ZZ;
    KNArray1D<size_t> Re, Im;
    KNSparseBlockMatrix  AHAT;
    KNVector       A_p;
    KNSparseMatrix     A_x;
    KNSparseMatrix     mB;
    KNVector       mB_p;
    KNSparseMatrix     mB_x;
    KNVector       phi;
    KNVector       DpPhi;
    KNVector       DxPhi;
    KNVector       LAM;
    KNVector       DxLAM;
    KNVector       uu3;   // for computing the generalized eigenvector
    KNVector       vv3;
    KNVector       gg3;
    KNVector       hh3;
    KNVector       rhs3;
    KNVector       one3;  // this is ( 0.0, 1.0 )
    KNVector       temp;
    KNArray3D<double>  phiData;
    KNArray3D<double>  LAMData;
    KNArray3D<double>  vv3Data;
    KNArray3D<double>  solMSHData;
};


class KNLpAutRotTestFunctional2 : public KNAbstractTestFunctional
{
  public:
    KNLpAutRotTestFunctional2(KNDdeBvpCollocation& col, KNArray1D<size_t> CRe, KNArray1D<size_t> CIm, double Z);
    ~KNLpAutRotTestFunctional2();
    void   init(KNDdeBvpCollocation& col, const KNVector& par, const KNVector& sol);
    double funct(KNDdeBvpCollocation& col, const KNVector& par, const KNVector& sol);
    double funct_p(KNDdeBvpCollocation& col, const KNVector& par, const KNVector& sol, size_t alpha);
    void   funct_x(KNVector& func, KNDdeBvpCollocation& col, const KNVector& par, const KNVector& sol);

    void init(KNDdeBvpCollocation&, const KNVector&, const KNVector&, double, double)
    {}
    double initStep()
    { return 0.0; }
    void funct(double&, double&, KNDdeBvpCollocation&, const KNVector&, const KNVector&, double, double)
    {}
    void funct_p(double&, double&, KNDdeBvpCollocation&, const KNVector&, const KNVector&, size_t)
    {}
    void funct_z(double&, double&, KNDdeBvpCollocation&, const KNVector&, const KNVector&)
    {}
    void funct_x(KNVector&, KNVector&, KNDdeBvpCollocation&, const KNVector&, const KNVector&)
    {}
    void kernel(KNVector& phi);

  private:
    bool         first;
    double       ZZ;
    KNArray1D<size_t> Re, Im;
    KNSparseBlockMatrix  AHAT;
    KNVector       A_p;
    KNSparseMatrix     A_x;
    KNSparseMatrix     mB;
    KNVector       mB_p;
    KNSparseMatrix     mB_x;
    KNVector       phi;
    KNVector       DpPhi;
    KNVector       DxPhi;
    KNVector       LAM;
    KNVector       DxLAM;
    KNVector       uu3;   // for computing the generalized eigenvector
    KNVector       vv3;
    KNVector       gg3;
    KNVector       hh3;
    KNVector       one3;  // this is ( 0.0, 0.0, 1.0 )
    KNVector       uu1;   // for computing the generalized eigenvector
    KNVector       vv1;
    KNVector       gg1;
    KNVector       hh1;
    KNVector       one1;  // this is ( 1.0, 0.0 )
    KNVector       uu2;   // for computing the generalized eigenvector
    KNVector       vv2;
    KNVector       gg2;
    KNVector       hh2;
    KNVector       one2;  // this is ( 0.0, 1.0 )
    KNVector       rhs;   // this is always zero, so one instance is enogh
    KNVector       temp;
    KNArray3D<double>  vv1Data;
    KNArray3D<double>  vv2Data;
    KNArray3D<double>  vv3Data;
    KNArray3D<double>  solMSHData;
};

class TestFunctIntersect : public KNAbstractTestFunctional
{
  public:
    TestFunctIntersect() { }
    virtual        ~TestFunctIntersect()
    {}
    // only first derivatives
    virtual void   init(KNDdeBvpCollocation& col, const KNVector& par, const KNVector& sol);
    virtual double initStep();
    // return value
    virtual double funct(KNDdeBvpCollocation& col, const KNVector& par, const KNVector& sol);
    // derivative w.r.t. alpha
    virtual double funct_p(KNDdeBvpCollocation& col, const KNVector& par, const KNVector& sol, size_t alpha);
    // derivative w.r.t. state variables
    virtual void   funct_x(KNVector& func, KNDdeBvpCollocation& col, const KNVector& sol, const KNVector& par);

    // these are not implemented
    virtual void init(KNDdeBvpCollocation& , const KNVector& , const KNVector& ,
                      double , double ) { }
    virtual void funct(double& , double& ,
                       KNDdeBvpCollocation& , const KNVector& , const KNVector& , double , double ) { }
    virtual void funct_p(double& , double& ,
                         KNDdeBvpCollocation& , const KNVector& , const KNVector& ,
                         size_t ) { }
    virtual void funct_z(double& , double& ,
                         KNDdeBvpCollocation& , const KNVector& , const KNVector& ) { }
    virtual void funct_x(KNVector& , KNVector& ,
                         KNDdeBvpCollocation& , const KNVector& , const KNVector& ) { }
    virtual void kernel(KNVector& ) { }
};

class TestFunctGrazing : public KNAbstractTestFunctional
{
  public:
    TestFunctGrazing() { }
    virtual        ~TestFunctGrazing()
    {}
    // set up second derivatives
    virtual void   init(KNDdeBvpCollocation& col, const KNVector& par, const KNVector& sol);
    virtual double initStep();
    // return value
    virtual double funct(KNDdeBvpCollocation& col, const KNVector& par, const KNVector& sol);
    // derivative w.r.t. alpha
    virtual double funct_p(KNDdeBvpCollocation& col, const KNVector& par, const KNVector& sol, size_t alpha);
    // derivative w.r.t. state variables
    virtual void   funct_x(KNVector& func, KNDdeBvpCollocation& col, const KNVector& sol, const KNVector& par);

    // these are not implemented
    virtual void init(KNDdeBvpCollocation& , const KNVector& , const KNVector& ,
                      double , double ) { }
    virtual void funct(double& , double& ,
                       KNDdeBvpCollocation& , const KNVector& , const KNVector& , double , double ) { }
    virtual void funct_p(double& , double& ,
                         KNDdeBvpCollocation& , const KNVector& , const KNVector& ,
                         size_t ) { }
    virtual void funct_z(double& , double& ,
                         KNDdeBvpCollocation& , const KNVector& , const KNVector& ) { }
    virtual void funct_x(KNVector& , KNVector& ,
                         KNDdeBvpCollocation& , const KNVector& , const KNVector& ) { }
    virtual void kernel(KNVector& ) { }
};

#endif
