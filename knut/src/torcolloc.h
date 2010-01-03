// ------------------------------------------------------------------------- //
//
// This is part of KNUT
// Copyright (c) 2002, 2003, 2004, 2005 by Robert Szalai
//
// For license, see the file COPYING in the package's root directory
//
// ------------------------------------------------------------------------- //

#ifndef TORCOLLOC_H
#define TORCOLLOC_H

#include "matrix.h"
#include "basecolloc.h"

class KNSystem;
class KNSparseMatrix;

class KNDdeTorusCollocation : public KNAbstractCollocation
{
  public:
    KNDdeTorusCollocation(KNSystem& sys_, int ndeg1_, int ndeg2_, int nint1_, int nint2_);
    // this provides the jacobian, the right hand side and the derivatives w.r.t. var
    // the difficulty is with the derivative w.r.t. the frequencies
    void init(const KNVector& sol, const KNVector& par);
    void jacobian(KNSparseMatrix& A, KNArray1D< KNVector* > Avar, KNVector& rhs, KNVector& par, KNVector& sol, KNArray1D<int>& var);
    // not yet implemented
    // these are easy
    void PhaseONE(KNVector& ph, KNVector& presol);                                       // implemented
    void PhaseBOTH(KNVector& ph0, KNVector& ph1, KNVector& presol);                        // implemented
    double integrate(const KNVector& ph1, const KNVector& ph2);
    double IntegrateDIFF(KNVector& ph1, KNVector& ph2, KNVector& ph3);
    void star(KNVector& ph1, const KNVector& ph2);

    void importSolution(KNVector& out, KNVector& in);
    void importTangent(KNVector& out, KNVector& Re, KNVector& Im, double alpha);
    void Save(const char* dat, const char* idx, const KNVector& in);
    void meshAdapt(KNVector& newsol, const KNVector& sol, KNVector& newtan, const KNVector& tan) { newsol = sol; newtan = tan; }

    inline const KNVector& getMesh1() const
    {
      return mesh1;
    }
    inline const KNVector& getMesh2() const
    {
      return mesh2;
    }

    inline int Ndeg1()
    {
      return ndeg1;
    }
    inline int Ndeg2()
    {
      return ndeg2;
    }
    inline int Nint1()
    {
      return nint1;
    }
    inline int Nint2()
    {
      return nint2;
    }

    //utils
    inline int idxmap(int j1, int j2, int i1, int i2);
    inline int idxkk(int j1, int j2, int k);
    inline KNSystem& system()
    {
      return *sys;
    }

  private:
    KNSystem* sys;
    const int ndim, ntau, npar;
    const int ndeg1, ndeg2;
    const int nint1, nint2;
    KNVector col1,  col2;
    KNVector mesh1, mesh2;
    KNArray1D< KNArray1D<double> > lgr1, lgr2; // 1. meshpoint 2. polynom
    KNArray1D< KNArray1D<double> > dlg1, dlg2; // 1. meshpoint 2. polynom

    KNMatrix I1, ID1, I2, ID2;
    KNVector mlg1, mlg2, mlgd1, mlgd2, ilg1, ilg2, ilgd1, ilgd2;
    // for the vectorization
    KNVector          time1, time2;   // ndeg1*ndeg2*nint1*nint2
    KNArray2D<int>    kk, ee, rr;     // (ntau+1)*(ndeg1+1)*(ndeg2+1) X ndeg1*ndeg2*nint1*nint2
    KNArray2D<double> p_tau, p_dtau;  // ntau X ndeg1*ndeg2*nint1*nint2
    KNArray3D<double> p_xx;           // ndim X ntau+1 X ndeg1*ndeg2*nint1*nint2
    KNArray2D<double> p_fx;           // ndim X ndeg1*ndeg2*nint1*nint2
    KNArray3D<double> p_dfp;          // ndim X 1 X ndeg1*ndeg2*nint1*nint2
    KNArray3D<double> p_dfx;          // ndim X ndim X ndeg1*ndeg2*nint1*nint2
    KNArray3D<double> p_dummy;
    // functions
};

inline int KNDdeTorusCollocation::idxmap(int j1, int j2, int i1, int i2)
{
  if (j1 < ndeg1)
  {
    if (j2 < ndeg2)
    {
      // j1 < NDEG1
      // j2 < NDEG2
      return j1 + j2*ndeg1 + i1*ndeg1*ndeg2 + i2*ndeg1*ndeg2*nint1;
    }
    else
    {
      // j1 < NDEG1
      // j2 = NDEG2
      return j1 + 0*ndeg1 + i1*ndeg1*ndeg2 + ((i2 + 1) % nint2)*ndeg1*ndeg2*nint1;
    }
  }
  else
  {
    if (j2 < ndeg2)
    {
      // j1 = NDEG1
      // j2 < NDEG2
      return 0 + j2*ndeg1 + ((i1 + 1) % nint1)*ndeg1*ndeg2 + i2*ndeg1*ndeg2*nint1;
    }
    else
    {
      // j1 = NDEG1
      // j2 = NDEG2
      if ((i1 == nint1 - 1) && (i2 == nint2 - 1)) return 0;
      else return 0 + 0*ndeg1 + ((i1 + 1) % nint1)*ndeg1*ndeg2 + ((i2 + 1) % nint2)*ndeg1*ndeg2*nint1;
    }
  }
}

inline int KNDdeTorusCollocation::idxkk(int j1, int j2, int k)
{
  return j1 + j2*(ndeg1 + 1) + k*(ndeg1 + 1)*(ndeg2 + 1);
}

#endif
