// ------------------------------------------------------------------------- //
//
// This is part of KNUT
// Copyright (c) 2002, 2003, 2004, 2005 by Robert Szalai
//
// For license, see the file COPYING in the package's root directory
//
// ------------------------------------------------------------------------- //

#ifndef NCOLLOC_H
#define NCOLLOC_H

#include "basecolloc.h"
#include "matrix.h"
#include "spmatrix.h"
#include "exprsystem.h"

class KNDdeBvpCollocation : public KNAbstractBvpCollocation
{
  public:

    KNDdeBvpCollocation(KNExprSystem& _sys, const size_t nint, const size_t ndeg);        // computes mesh, metric, metricD

    ~KNDdeBvpCollocation() override {}

    void init(const KNVector& sol, const KNVector& par) override;   // computes time, kk, ee, dd, rr ...

    void interpolate(KNArray3D<double>& out, const KNVector& sol) override;
    void interpolateComplex(KNArray3D<double>& outRe, KNArray3D<double>& outIm, const KNVector& sol) override;
    void interpolateOnMesh(KNArray3D<double>& out, const KNVector& sol) override;

    // computing the Jacobians, right-hand sides, characteristic matrices

    // continuation of solution

    void rightHandSide(KNVector& rhs, const KNVector& par, const KNVector& sol) override;
    void rightHandSide_p(KNVector& rhs, const KNVector& par, const KNVector& sol, size_t p) override;   // sol is currently not needed
    void rightHandSide_x(KNSparseMatrix& A, const KNVector& par, const KNVector& sol) override;          // sol is currently not needed

    // for stability computation
    void jacobianOfStability(KNSparseMatrixPolynomial& AB, const KNVector& par);
    // this computes its own matrix structures, because it is nowhere else needed: kkSI, eeSI, rrSI, ddSI, etc.
    // However, the variables will be contained within the KNDdeBvpCollocation class

    // continuation of bifurcations -> characteristic matrices

    inline size_t nTau() const { return ntau; }
    inline size_t nMat() const { return nmat; }

    // Jacobian of the test funtional: jotf
    void jotf_x(KNSparseMatrix& A, const KNVector& par, double Z);
    void jotf_x(KNSparseMatrix& A, const KNVector& par, double ZRe, double ZIm);

    void jotf_x_p(KNVector& V, const KNVector& par, const KNArray3D<double>& phiData, double Z, size_t p);
    void jotf_x_p(KNVector& V, const KNVector& par,
      const KNArray3D<double>& phiDataRe, const KNArray3D<double>& phiDataIm, double ZRe, double ZIm, size_t p);

    void jotf_x_x(KNSparseMatrix& A, const KNVector& par, const KNArray3D<double>& phiData, double Z);
    void jotf_x_x(KNSparseMatrix& A, const KNVector& par,
      const KNArray3D<double>& phiDataRe, const KNArray3D<double>& phiDataIm, double Re, double Im);

    void jotf_x_z(KNVector& V, const KNVector& par, const KNVector& phi,
      const KNArray3D<double>& phiDataRe, const KNArray3D<double>& phiDataIm, double Re, double Im);
    // we do need Re, Im

    // for the autonomous FOLD bifurcation
    void jotf_mB(KNSparseMatrix& A, const KNVector& par, double Z);
    void jotf_mB_p(KNVector& V,   const KNVector& par, const KNArray3D<double>& phiData, double Z, size_t alpha);
    void jotf_mB_x(KNSparseMatrix& A, const KNVector& par, const KNArray3D<double>& phiData, double Z);

    // autonomous trivial eigenvector
    void jotf_trivialKernelOnMesh(KNVector& V, const KNVector& par, const KNArray3D<double>& solMSHData);
    template<bool trans> void jotf_trivialKernelOnMesh_x(KNVector& V, const KNVector& par, const KNArray3D<double>& solMSHData, const KNVector& phi);
    void jotf_trivialKernelOnMesh_p(KNVector& V, const KNVector& par, const KNArray3D<double>& solMSHData, size_t alpha);

  private:

    // rINT and rDEG is included in idx. Only rDIM is necessary
    inline int& WRIDX(KNSparseMatrix& A, size_t idx, size_t rDIM, size_t cTAU, size_t cDEG, size_t cDIM)
    {
      return A.writeIndex(ndim + ndim*(idx) + rDIM, cDIM + ndim*(cDEG - dd(cTAU, idx) + rr(cTAU, idx)*(ndeg + 1)));
    }

    inline double& WRDAT(KNSparseMatrix& A, size_t idx, size_t rDIM, size_t cTAU, size_t cDEG, size_t cDIM)
    {
      return A.writeData(ndim + ndim*(idx) + rDIM, cDIM + ndim*(cDEG - dd(cTAU, idx) + rr(cTAU, idx)*(ndeg + 1)));
    }

    // for CharJac_x
    inline int& WRIDXCPLX(KNSparseMatrix& A, size_t idx, size_t rDIM, size_t rIM, size_t cTAU, size_t cDEG, size_t cDIM, size_t cIM)
    {
      return A.writeIndex(rIM + 2*(ndim + ndim*(idx) + rDIM), cIM + 2*(cDIM + ndim*(cDEG - dd(cTAU, idx) + rr(cTAU, idx)*(ndeg + 1))));
    }

    inline double& WRDATCPLX(KNSparseMatrix& A, size_t idx, size_t rDIM, size_t rIM, size_t cTAU, size_t cDEG, size_t cDIM, size_t cIM)
    {
      return A.writeData(rIM + 2*(ndim + ndim*(idx) + rDIM), cIM + 2*(cDIM + ndim*(cDEG - dd(cTAU, idx) + rr(cTAU, idx)*(ndeg + 1))));
    }

    // for CharJac_x_x
    inline int& WRIDXCPLXM(KNSparseMatrix& A, size_t idx, size_t rDIM, size_t rIM, size_t cTAU, size_t cDEG, size_t cDIM)
    {
      return A.writeIndex(rIM + 2*(ndim + ndim*(idx) + rDIM), cDIM + ndim*(cDEG - dd(cTAU, idx) + rr(cTAU, idx)*(ndeg + 1)));
    }

    inline double& WRDATCPLXM(KNSparseMatrix& A, size_t idx, size_t rDIM, size_t rIM, size_t cTAU, size_t cDEG, size_t cDIM)
    {
      return A.writeData(rIM + 2*(ndim + ndim*(idx) + rDIM), cDIM + ndim*(cDEG - dd(cTAU, idx) + rr(cTAU, idx)*(ndeg + 1)));
    }

    inline int& WRIDXS(KNSparseMatrix& A, size_t idx, size_t rDIM, size_t cTAU, size_t cDEG, size_t cDIM)
    {
      return A.writeIndex(ndim + ndim*(idx) + rDIM, cDIM + ndim*(cDEG - ddS(cTAU, idx) + rrS(cTAU, idx)*(ndeg + 1)));
    }

    double& WRDATS(KNSparseMatrix& A, size_t idx, size_t rDIM, size_t cTAU, size_t cDEG, size_t cDIM)
    {
      return A.writeData(ndim + ndim*(idx) + rDIM, cDIM + ndim*(cDEG - ddS(cTAU, idx) + rrS(cTAU, idx)*(ndeg + 1)));
    }

    inline int& WRIDXI(KNSparseMatrix& A, size_t idx, size_t rDIM, size_t cTAU, size_t cDEG, size_t cDIM)
    {
      return A.writeIndex(ndim + ndim*(idx) + rDIM, cDIM + ndim*(cDEG - ddI(cTAU, idx) + rrI(cTAU, idx)*(ndeg + 1)));
    }

    inline double& WRDATI(KNSparseMatrix& A, size_t idx, size_t rDIM, size_t cTAU, size_t cDEG, size_t cDIM)
    {
      return A.writeData(ndim + ndim*(idx) + rDIM, cDIM + ndim*(cDEG - ddI(cTAU, idx) + rrI(cTAU, idx)*(ndeg + 1)));
    }

    const size_t ntau;
    size_t nmat; // is the ceil( max(delays)/period )

    // these store the structure of the sparse matrix NO SEPARATION
    KNArray2D<size_t> kk;   // dim(NTAU+1,NDEG*NINT) : which delay to which interval
    KNArray2D<size_t> ee;   // dim(NTAU+1,NDEG*NINT) : the ordering of kk
    KNArray2D<size_t> dd;   // dim(NTAU+1,NDEG*NINT) : how many neighbouring intervals up to its index
    KNArray2D<size_t> rr;   // dim(NTAU+1,NDEG*NINT) : how many non overlapping intervals up to its index

    // these store the structure of the sparse matrix WITH SEPARATION
    KNArray2D<size_t> kkS;  // same as kk, but not folded back
    KNArray2D<size_t> eeS;  //  ...
    KNArray2D<size_t> rrS;  //  ...
    KNArray2D<size_t> ddS;  //  ...
    KNArray2D<size_t> mmS;  // which matrix are we in
    KNArray2D<size_t> szS;  // szI(mmI( . ,idx),idx) size of the line within the matrix mmS

    // For the stability matrices
    KNArray2D<size_t> kkI;
    KNArray2D<size_t> eeI;
    KNArray2D<size_t> rrI;
    KNArray2D<size_t> ddI;
    KNArray2D<size_t> mmI;
    KNArray2D<size_t> szI;

    KNArray2D<size_t> kkMSH; // dim(NTAU+1,NDEG*NINT+1) : which delay to which interval

    KNArray3D<double> tt;       // interpolation at the collocation points
    KNArray3D<double> ttMSH;    // interpolation at the representation points. I guess its the identity

    // for rhs, and derivatives
    // KNMatrix fx, dfx, t_dfx, dummy
    KNArray3D<double>  solData;
    KNArray2D<double>  p_tau;
    KNArray2D<double>  p_dtau;
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

#ifdef DEBUG
    // DIAGNOSTICS
    size_t count_RHS;
    size_t count_RHS_p;
    size_t count_RHS_x;
    size_t count_StabJac;
    size_t count_CharJac_x;
    size_t count_CharJac_x_p;
    size_t count_CharJac_x_x;
    size_t count_CharJac_x_z;
    size_t count_CharJac_mB;
    size_t count_CharJac_mB_p;
    size_t count_CharJac_mB_x;
    size_t count_CharJac_MSHphi;
    size_t count_CharJac_MSHphi_p;
    size_t count_CharJac_MSHphi_x;
    // DIAG FUNCTIONS
    void count_reset()
    {
      count_RHS = 0;
      count_RHS_p = 0;
      count_RHS_x = 0;
      count_StabJac = 0;
      count_CharJac_x = 0;
      count_CharJac_x_p = 0;
      count_CharJac_x_x = 0;
      count_CharJac_x_z = 0;
      count_CharJac_mB = 0;
      count_CharJac_mB_p = 0;
      count_CharJac_mB_x = 0;
      count_CharJac_MSHphi = 0;
      count_CharJac_MSHphi_p = 0;
      count_CharJac_MSHphi_x = 0;
    }
    void count_print()
    {
      if (count_RHS) std::cout << "count_RHS " << count_RHS << "\n";
      if (count_RHS_p) std::cout << "count_RHS_p " << count_RHS_p << "\n";
      if (count_RHS_x) std::cout << "count_RHS_x " << count_RHS_x << "\n";
      if (count_StabJac) std::cout << "count_StabJac " << count_StabJac << "\n";
      if (count_CharJac_x) std::cout << "count_CharJac_x " << count_CharJac_x << "\n";
      if (count_CharJac_x_p) std::cout << "count_CharJac_x_p " << count_CharJac_x_p << "\n";
      if (count_CharJac_x_x) std::cout << "count_CharJac_x_x " << count_CharJac_x_x << "\n";
      if (count_CharJac_x_z) std::cout << "count_CharJac_x_z " << count_CharJac_x_z << "\n";
      if (count_CharJac_mB) std::cout << "count_CharJac_mB " << count_CharJac_mB << "\n";
      if (count_CharJac_mB_p) std::cout << "count_CharJac_mB_p " << count_CharJac_mB_p << "\n";
      if (count_CharJac_mB_x) std::cout << "count_CharJac_mB_x " << count_CharJac_mB_x << "\n";
      if (count_CharJac_MSHphi) std::cout << "count_CharJac_MSHphi " << count_CharJac_MSHphi << "\n";
      if (count_CharJac_MSHphi_p) std::cout << "count_CharJac_MSHphi_p " << count_CharJac_MSHphi_p << "\n";
      if (count_CharJac_MSHphi_x) std::cout << "count_CharJac_MSHphi_x " << count_CharJac_MSHphi_x << "\n";
    }
#endif // DEBUG
};

#define NDIM ndim
#define NPAR npar
#define NTAU ntau
#define NINT nint
#define NDEG ndeg
#define NMAT nmat

template<bool trans> void KNDdeBvpCollocation::jotf_trivialKernelOnMesh_x(KNVector& V, const KNVector& par, const KNArray3D<double>& solMSHData, const KNVector& phi)
{
#ifdef DEBUG
  count_CharJac_MSHphi_x++;
#endif
  V.clear(); // it is not cleared otherwise!!!!

  for (size_t k = 0; k < NTAU; k++)
  {
    size_t nx = 1, vx = k, np = 0, vp = 0;
    sys->p_deri(p_dfxMSH, timeMSH, solMSHData, par, 0, nx, &vx, np, &vp, p_dummy);
    for (size_t idx = 0; idx < NINT*NDEG + 1; idx++) // i: interval; j: which collocation point
    {
      for (size_t l = 0; l < NDEG + 1; l++) // degree
      {
        for (size_t p = 0; p < NDIM; p++)   // row
        {
          for (size_t q = 0; q < NDIM; q++)   //column
          {
            if (trans) V(q + NDIM*(l + NDEG*kkMSH(k, idx))) -=
                p_dfxMSH(p, q, idx) * ttMSH(k, l, idx) * phi(p + NDIM * idx);
            else      V(p + NDIM*idx) -=
                p_dfxMSH(p, q, idx) * ttMSH(k, l, idx) * phi(q + NDIM * (l + NDEG * kkMSH(k, idx)));
          }
        }
      }
    }
  }
}

#undef NDIM
#undef NPAR
#undef NTAU
#undef NINT
#undef NDEG
#undef NMAT

#endif
