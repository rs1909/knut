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
    virtual double initStep(KNDdeBvpCollocation&) = 0;
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
    ~KNTestFunctional() override;
    void   init(KNDdeBvpCollocation& col, const KNVector& par, const KNVector& sol) override;
    double initStep(KNDdeBvpCollocation&) override;
    double funct(KNDdeBvpCollocation& col, const KNVector& par, const KNVector& sol) override;
    double funct_p(KNDdeBvpCollocation& col, const KNVector& par, const KNVector& sol, size_t alpha) override;
    void   funct_x(KNVector& func, KNDdeBvpCollocation& col, const KNVector& par, const KNVector& sol) override;

    void init(KNDdeBvpCollocation&, const KNVector&, const KNVector&, double, double) override
    {}
    void funct(double&, double&, KNDdeBvpCollocation&, const KNVector&, const KNVector&, double, double) override
    {}
    void funct_p(double&, double&, KNDdeBvpCollocation&, const KNVector&, const KNVector&, size_t) override
    {}
    void funct_z(double&, double&, KNDdeBvpCollocation&, const KNVector&, const KNVector&) override
    {}
    void funct_x(KNVector&, KNVector&, KNDdeBvpCollocation&, const KNVector&, const KNVector&) override
    {}
    void kernel(KNVector& phi) override;

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
    ~KNComplexTestFunctional() override;
    void   init(KNDdeBvpCollocation&, const KNVector&, const KNVector&) override
    {}
    double funct(KNDdeBvpCollocation&, const KNVector&, const KNVector&) override
    {
      return 0.0;
    }
    double funct_p(KNDdeBvpCollocation&, const KNVector&, const KNVector&, size_t) override
    {
      return 0.0;
    }
    void   funct_x(KNVector&, KNDdeBvpCollocation&, const KNVector&, const KNVector&) override
    {}

    void init(KNDdeBvpCollocation& col, const KNVector& par, const KNVector& sol,
              double Re, double Im) override;
    double initStep(KNDdeBvpCollocation&) override;
    void funct(double& f1, double& f2,
               KNDdeBvpCollocation& col, const KNVector& par, const KNVector& sol, double Re, double Im) override;
    void funct_p(double& f1, double& f2,
                 KNDdeBvpCollocation& col, const KNVector& par, const KNVector& sol,
                 size_t alpha) override;
    void funct_z(double& f1, double& f2,
                 KNDdeBvpCollocation& col, const KNVector& par, const KNVector& sol) override;
    void funct_x(KNVector& func1, KNVector& func2,
                 KNDdeBvpCollocation& col, const KNVector& par, const KNVector& sol) override;
    inline void kernel(KNVector&) override
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
    ~KNLpAutTestFunctional() override;
    void   init(KNDdeBvpCollocation& col, const KNVector& par, const KNVector& sol) override;
    double initStep(KNDdeBvpCollocation&) override;
    double funct(KNDdeBvpCollocation& col, const KNVector& par, const KNVector& sol) override;
    double funct_p(KNDdeBvpCollocation& col, const KNVector& par, const KNVector& sol, size_t alpha) override;
    void   funct_x(KNVector& func, KNDdeBvpCollocation& col, const KNVector& sol, const KNVector& par) override;

    void init(KNDdeBvpCollocation&, const KNVector&, const KNVector&, double, double) override
    {}
    void funct(double&, double&, KNDdeBvpCollocation&, const KNVector&, const KNVector&, double, double) override
    {}
    void funct_p(double&, double&, KNDdeBvpCollocation&, const KNVector&, const KNVector&, size_t) override
    {}
    void funct_z(double&, double&, KNDdeBvpCollocation&, const KNVector&, const KNVector&) override
    {}
    void funct_x(KNVector&, KNVector&, KNDdeBvpCollocation&, const KNVector&, const KNVector&) override
    {}
    void kernel(KNVector& phi) override;

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
    ~KNLpAutRotTestFunctional() override;
    void   init(KNDdeBvpCollocation& col, const KNVector& par, const KNVector& sol) override;
    double initStep(KNDdeBvpCollocation&) override;
    double funct(KNDdeBvpCollocation& col, const KNVector& par, const KNVector& sol) override;
    double funct_p(KNDdeBvpCollocation& col, const KNVector& par, const KNVector& sol, size_t alpha) override;
    void   funct_x(KNVector& func, KNDdeBvpCollocation& col, const KNVector& par, const KNVector& sol) override;

    void init(KNDdeBvpCollocation&, const KNVector&, const KNVector&, double, double) override
    {}
    void funct(double&, double&, KNDdeBvpCollocation&, const KNVector&, const KNVector&, double, double) override
    {}
    void funct_p(double&, double&, KNDdeBvpCollocation&, const KNVector&, const KNVector&, size_t) override
    {}
    void funct_z(double&, double&, KNDdeBvpCollocation&, const KNVector&, const KNVector&) override
    {}
    void funct_x(KNVector&, KNVector&, KNDdeBvpCollocation&, const KNVector&, const KNVector&) override
    {}
    void kernel(KNVector& phi) override;

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
    ~KNLpAutRotTestFunctional2() override;
    void   init(KNDdeBvpCollocation& col, const KNVector& par, const KNVector& sol) override;
    double funct(KNDdeBvpCollocation& col, const KNVector& par, const KNVector& sol) override;
    double funct_p(KNDdeBvpCollocation& col, const KNVector& par, const KNVector& sol, size_t alpha) override;
    void   funct_x(KNVector& func, KNDdeBvpCollocation& col, const KNVector& par, const KNVector& sol) override;

    void init(KNDdeBvpCollocation&, const KNVector&, const KNVector&, double, double) override
    {}
    double initStep(KNDdeBvpCollocation&) override
    { return 0.0; }
    void funct(double&, double&, KNDdeBvpCollocation&, const KNVector&, const KNVector&, double, double) override
    {}
    void funct_p(double&, double&, KNDdeBvpCollocation&, const KNVector&, const KNVector&, size_t) override
    {}
    void funct_z(double&, double&, KNDdeBvpCollocation&, const KNVector&, const KNVector&) override
    {}
    void funct_x(KNVector&, KNVector&, KNDdeBvpCollocation&, const KNVector&, const KNVector&) override
    {}
    void kernel(KNVector& phi) override;

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
           ~TestFunctIntersect() override
    {}
    // only first derivatives
    void   init(KNDdeBvpCollocation& col, const KNVector& par, const KNVector& sol) override;
    double initStep(KNDdeBvpCollocation&) override;
    // return value
    double funct(KNDdeBvpCollocation& col, const KNVector& par, const KNVector& sol) override;
    // derivative w.r.t. alpha
    double funct_p(KNDdeBvpCollocation& col, const KNVector& par, const KNVector& sol, size_t alpha) override;
    // derivative w.r.t. state variables
    void   funct_x(KNVector& func, KNDdeBvpCollocation& col, const KNVector& sol, const KNVector& par) override;

    // these are not implemented
    void init(KNDdeBvpCollocation& , const KNVector& , const KNVector& ,
                      double , double ) override { }
    void funct(double& , double& ,
                       KNDdeBvpCollocation& , const KNVector& , const KNVector& , double , double ) override { }
    void funct_p(double& , double& ,
                         KNDdeBvpCollocation& , const KNVector& , const KNVector& ,
                         size_t ) override { }
    void funct_z(double& , double& ,
                         KNDdeBvpCollocation& , const KNVector& , const KNVector& ) override { }
    void funct_x(KNVector& , KNVector& ,
                         KNDdeBvpCollocation& , const KNVector& , const KNVector& ) override { }
    void kernel(KNVector& ) override { }
};

class TestFunctGrazing : public KNAbstractTestFunctional
{
  public:
    TestFunctGrazing() { }
           ~TestFunctGrazing() override
    {}
    // set up second derivatives
    void   init(KNDdeBvpCollocation& col, const KNVector& par, const KNVector& sol) override;
    double initStep(KNDdeBvpCollocation&) override;
    // return value
    double funct(KNDdeBvpCollocation& col, const KNVector& par, const KNVector& sol) override;
    // derivative w.r.t. alpha
    double funct_p(KNDdeBvpCollocation& col, const KNVector& par, const KNVector& sol, size_t alpha) override;
    // derivative w.r.t. state variables
    void   funct_x(KNVector& func, KNDdeBvpCollocation& col, const KNVector& sol, const KNVector& par) override;

    // these are not implemented
    void init(KNDdeBvpCollocation& , const KNVector& , const KNVector& ,
                      double , double ) override { }
    void funct(double& , double& ,
                       KNDdeBvpCollocation& , const KNVector& , const KNVector& , double , double ) override { }
    void funct_p(double& , double& ,
                         KNDdeBvpCollocation& , const KNVector& , const KNVector& ,
                         size_t ) override { }
    void funct_z(double& , double& ,
                         KNDdeBvpCollocation& , const KNVector& , const KNVector& ) override { }
    void funct_x(KNVector& , KNVector& ,
                         KNDdeBvpCollocation& , const KNVector& , const KNVector& ) override { }
    void kernel(KNVector& ) override { }
};

#endif
