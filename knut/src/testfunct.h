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
class NColloc;

#define TF_NKERNITER 20
#define TF_KERNEPS 1.0e-12

class baseTestFunct
{
  public:
    baseTestFunct() { kernEps = TF_KERNEPS; kernIter = TF_NKERNITER; }
    virtual        ~baseTestFunct()
    {}
    virtual void   init(NColloc& col, const Vector& par, const Vector& sol) = 0;
    virtual double Funct(NColloc& col, const Vector& par, const Vector& sol) = 0;
    virtual double Funct_p(NColloc& col, const Vector& par, const Vector& sol, int alpha) = 0;
    virtual void   Funct_x(Vector& func, NColloc& col, const Vector& sol, const Vector& par) = 0;

    virtual void init(NColloc& col, const Vector& par, const Vector& sol,
                      double Re, double Im) = 0;
    virtual void Funct(double& f1, double& f2,
                       NColloc& col, const Vector& par, const Vector& sol, double Re, double Im) = 0;
    virtual void Funct_p(double& f1, double& f2,
                         NColloc& col, const Vector& par, const Vector& sol,
                         int alpha) = 0;
    virtual void Funct_z(double& f1, double& f2,
                         NColloc& col, const Vector& par, const Vector& sol) = 0;
    virtual void Funct_x(Vector& func1, Vector& func2,
                         NColloc& col, const Vector& par, const Vector& sol) = 0;
    virtual void Switch(Vector& phi) = 0;
    virtual void setKernelTolerance(double eps, int iter) { kernEps = eps; kernIter = iter; }
  protected:
    double kernEps;
    int    kernIter;
};

class TestFunct : public baseTestFunct
{
  public:
    TestFunct(NColloc& col, double Z);
    ~TestFunct();
    void   init(NColloc& col, const Vector& par, const Vector& sol);
    double Funct(NColloc& col, const Vector& par, const Vector& sol);
    double Funct_p(NColloc& col, const Vector& par, const Vector& sol, int alpha);
    void   Funct_x(Vector& func, NColloc& col, const Vector& par, const Vector& sol);

    void init(NColloc&, const Vector&, const Vector&, double, double)
    {}
    void Funct(double&, double&, NColloc&, const Vector&, const Vector&, double, double)
    {}
    void Funct_p(double&, double&, NColloc&, const Vector&, const Vector&, int)
    {}
    void Funct_z(double&, double&, NColloc&, const Vector&, const Vector&)
    {}
    void Funct_x(Vector&, Vector&, NColloc&, const Vector&, const Vector&)
    {}
    void Switch(Vector& phi);

  private:
    bool        first;
    double      ZZ;
    HyperMatrix AHAT;
    Vector      A_p;
    SpMatrix    A_x;
    Vector      rhs;
    Vector      uu;
    Vector      vv;
    Vector      uudiff;
    Vector      vvdiff;
    Array3D<double> vvData;
};

class TestFunctCPLX : public baseTestFunct
{
  public:
    TestFunctCPLX(NColloc& col);
    ~TestFunctCPLX();
    void   init(NColloc&, const Vector&, const Vector&)
    {}
    double Funct(NColloc&, const Vector&, const Vector&)
    {
      return 0.0;
    }
    double Funct_p(NColloc&, const Vector&, const Vector&, int)
    {
      return 0.0;
    }
    void   Funct_x(Vector&, NColloc&, const Vector&, const Vector&)
    {}

    void init(NColloc& col, const Vector& par, const Vector& sol,
              double Re, double Im);
    void Funct(double& f1, double& f2,
               NColloc& col, const Vector& par, const Vector& sol, double Re, double Im);
    void Funct_p(double& f1, double& f2,
                 NColloc& col, const Vector& par, const Vector& sol,
                 int alpha);
    void Funct_z(double& f1, double& f2,
                 NColloc& col, const Vector& par, const Vector& sol);
    void Funct_x(Vector& func1, Vector& func2,
                 NColloc& col, const Vector& par, const Vector& sol);
    inline void Switch(Vector&)
    {}
    void   Switch(Vector& Re, Vector& Im, double& alpha);
    double SwitchHB(Vector& Re, Vector& Im, NColloc& col, const Vector& par);

  private:
    bool        first;
    double      ZRe, ZIm;
    HyperMatrix AHAT;
    Vector      A_p;
    SpMatrix    A_x;
    Vector      rhs;
    Vector      one;
    Vector      uu;
    Vector      vv;
    Vector      gg;
    Vector      hh;
    Array3D<double> vvDataRe;
    Array3D<double> vvDataIm;
};

class TestFunctLPAUT : public baseTestFunct
{
  public:
    TestFunctLPAUT(NColloc& col, double Z);
    ~TestFunctLPAUT();
    void   init(NColloc& col, const Vector& par, const Vector& sol);
    double Funct(NColloc& col, const Vector& par, const Vector& sol);
    double Funct_p(NColloc& col, const Vector& par, const Vector& sol, int alpha);
    void   Funct_x(Vector& func, NColloc& col, const Vector& sol, const Vector& par);

    void init(NColloc&, const Vector&, const Vector&, double, double)
    {}
    void Funct(double&, double&, NColloc&, const Vector&, const Vector&, double, double)
    {}
    void Funct_p(double&, double&, NColloc&, const Vector&, const Vector&, int)
    {}
    void Funct_z(double&, double&, NColloc&, const Vector&, const Vector&)
    {}
    void Funct_x(Vector&, Vector&, NColloc&, const Vector&, const Vector&)
    {}
    void Switch(Vector& phi);

  private:
    bool        first;
    double      ZZ;
    HyperMatrix AHAT;
    Vector      A_p;
    SpMatrix    A_x;
    SpMatrix    mB;
    Vector      mB_p;
    SpMatrix    mB_x;
    Vector      rhs;   // this is zero all the time
    Vector      phi;
    Vector      temp;
    Vector      DpPhi;
    Vector      DxPhi;
    Vector      uu2;   // for computing the generalized eigenvector
    Vector      vv2;
    Vector      gg2;
    Vector      hh2;
    Vector      one2;  // this is ( 0.0, 1.0 )
    Array3D<double> phiData;
    Array3D<double> vv2Data;
    Array3D<double> solMSHData;
};

class TestFunctLPAUTROT : public baseTestFunct
{
  public:
    TestFunctLPAUTROT(NColloc& col, Array1D<int> CRe, Array1D<int> CIm, double Z);
    ~TestFunctLPAUTROT();
    void   init(NColloc& col, const Vector& par, const Vector& sol);
    double Funct(NColloc& col, const Vector& par, const Vector& sol);
    double Funct_p(NColloc& col, const Vector& par, const Vector& sol, int alpha);
    void   Funct_x(Vector& func, NColloc& col, const Vector& par, const Vector& sol);

    void init(NColloc&, const Vector&, const Vector&, double, double)
    {}
    void Funct(double&, double&, NColloc&, const Vector&, const Vector&, double, double)
    {}
    void Funct_p(double&, double&, NColloc&, const Vector&, const Vector&, int)
    {}
    void Funct_z(double&, double&, NColloc&, const Vector&, const Vector&)
    {}
    void Funct_x(Vector&, Vector&, NColloc&, const Vector&, const Vector&)
    {}
    void Switch(Vector& phi);

  private:
    bool         first;
    double       ZZ;
    Array1D<int> Re, Im;
    HyperMatrix  AHAT;
    Vector       A_p;
    SpMatrix     A_x;
    SpMatrix     mB;
    Vector       mB_p;
    SpMatrix     mB_x;
    Vector       phi;
    Vector       DpPhi;
    Vector       DxPhi;
    Vector       LAM;
    Vector       DxLAM;
    Vector       uu3;   // for computing the generalized eigenvector
    Vector       vv3;
    Vector       gg3;
    Vector       hh3;
    Vector       rhs3;
    Vector       one3;  // this is ( 0.0, 1.0 )
    Vector       temp;
    Array3D<double>  phiData;
    Array3D<double>  LAMData;
    Array3D<double>  vv3Data;
    Array3D<double>  solMSHData;
};


class TestFunctLPAUTROT_X : public baseTestFunct
{
  public:
    TestFunctLPAUTROT_X(NColloc& col, Array1D<int> CRe, Array1D<int> CIm, double Z);
    ~TestFunctLPAUTROT_X();
    void   init(NColloc& col, const Vector& par, const Vector& sol);
    double Funct(NColloc& col, const Vector& par, const Vector& sol);
    double Funct_p(NColloc& col, const Vector& par, const Vector& sol, int alpha);
    void   Funct_x(Vector& func, NColloc& col, const Vector& par, const Vector& sol);

    void init(NColloc&, const Vector&, const Vector&, double, double)
    {}
    void Funct(double&, double&, NColloc&, const Vector&, const Vector&, double, double)
    {}
    void Funct_p(double&, double&, NColloc&, const Vector&, const Vector&, int)
    {}
    void Funct_z(double&, double&, NColloc&, const Vector&, const Vector&)
    {}
    void Funct_x(Vector&, Vector&, NColloc&, const Vector&, const Vector&)
    {}
    void Switch(Vector& phi);

  private:
    bool         first;
    double       ZZ;
    Array1D<int> Re, Im;
    HyperMatrix  AHAT;
    Vector       A_p;
    SpMatrix     A_x;
    SpMatrix     mB;
    Vector       mB_p;
    SpMatrix     mB_x;
    Vector       phi;
    Vector       DpPhi;
    Vector       DxPhi;
    Vector       LAM;
    Vector       DxLAM;
    Vector       uu3;   // for computing the generalized eigenvector
    Vector       vv3;
    Vector       gg3;
    Vector       hh3;
    Vector       one3;  // this is ( 0.0, 0.0, 1.0 )
    Vector       uu1;   // for computing the generalized eigenvector
    Vector       vv1;
    Vector       gg1;
    Vector       hh1;
    Vector       one1;  // this is ( 1.0, 0.0 )
    Vector       uu2;   // for computing the generalized eigenvector
    Vector       vv2;
    Vector       gg2;
    Vector       hh2;
    Vector       one2;  // this is ( 0.0, 1.0 )
    Vector       rhs;   // this is always zero, so one instance is enogh
    Vector       temp;
    Array3D<double>  vv1Data;
    Array3D<double>  vv2Data;
    Array3D<double>  vv3Data;
    Array3D<double>  solMSHData;
};

class TestFunctIntersect : public baseTestFunct
{
  public:
    TestFunctIntersect() { }
    virtual        ~TestFunctIntersect()
    {}
    // only first derivatives
    virtual void   init(NColloc& col, const Vector& par, const Vector& sol);
    // return value
    virtual double Funct(NColloc& col, const Vector& par, const Vector& sol);
    // derivative w.r.t. alpha
    virtual double Funct_p(NColloc& col, const Vector& par, const Vector& sol, int alpha);
    // derivative w.r.t. state variables
    virtual void   Funct_x(Vector& func, NColloc& col, const Vector& sol, const Vector& par);

    // these are not implemented
    virtual void init(NColloc& col, const Vector& par, const Vector& sol,
                      double Re, double Im) { }
    virtual void Funct(double& f1, double& f2,
                       NColloc& col, const Vector& par, const Vector& sol, double Re, double Im) { }
    virtual void Funct_p(double& f1, double& f2,
                         NColloc& col, const Vector& par, const Vector& sol,
                         int alpha) { }
    virtual void Funct_z(double& f1, double& f2,
                         NColloc& col, const Vector& par, const Vector& sol) { }
    virtual void Funct_x(Vector& func1, Vector& func2,
                         NColloc& col, const Vector& par, const Vector& sol) { }
    virtual void Switch(Vector& phi) { }
};

class TestFunctGrazing : public baseTestFunct
{
  public:
    TestFunctGrazing() { }
    virtual        ~TestFunctGrazing()
    {}
    // set up second derivatives
    virtual void   init(NColloc& col, const Vector& par, const Vector& sol);
    // return value
    virtual double Funct(NColloc& col, const Vector& par, const Vector& sol);
    // derivative w.r.t. alpha
    virtual double Funct_p(NColloc& col, const Vector& par, const Vector& sol, int alpha);
    // derivative w.r.t. state variables
    virtual void   Funct_x(Vector& func, NColloc& col, const Vector& sol, const Vector& par);

    // these are not implemented
    virtual void init(NColloc& col, const Vector& par, const Vector& sol,
                      double Re, double Im) { }
    virtual void Funct(double& f1, double& f2,
                       NColloc& col, const Vector& par, const Vector& sol, double Re, double Im) { }
    virtual void Funct_p(double& f1, double& f2,
                         NColloc& col, const Vector& par, const Vector& sol,
                         int alpha) { }
    virtual void Funct_z(double& f1, double& f2,
                         NColloc& col, const Vector& par, const Vector& sol) { }
    virtual void Funct_x(Vector& func1, Vector& func2,
                         NColloc& col, const Vector& par, const Vector& sol) { }
    virtual void Switch(Vector& phi) { }
};

#endif
