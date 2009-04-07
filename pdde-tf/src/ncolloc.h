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
#include "system.h"
#include "matrix.h"
#include "spmatrix.h"

class NColloc : public BaseColloc
{
  public:

    NColloc(System& _sys, const int nint, const int ndeg, int nmat);        // computes mesh, metric, metricD

    ~NColloc()
    {}

    void Init(const Vector& sol, const Vector& par);   // computes time, kk, ee, dd, rr ...
    void meshAdapt(Vector& newprofile, const Vector& profile, Vector& newtangent, const Vector& tangent);

    void Interpolate(Array3D<double>& out, const Vector& sol);
    void InterpolateCPLX(Array3D<double>& outRe, Array3D<double>& outIm, const Vector& sol);
    void InterpolateMSH(Array3D<double>& out, const Vector& sol);

    static void   getMetric(Matrix& mt, const Vector& t);
    static void   getDiffMetric(Matrix& mt, const Vector& t);
    static void   star(Vector& out, const Vector& in, const Matrix& mt, const Vector& msh, int dim);
    static double integrate(const Vector& v1, const Vector& v2, const Matrix& mt, const Vector& msh, int dim);

    void   Star(Vector& out, const Vector& sol);
    double Integrate(const Vector& v1, const Vector& v2);
    double IntegrateDerivative(const Vector& v1, const Vector& v2);
    double IntegrateCont(const Vector& v1, const Vector& v2, const Vector& v3);

    void   PhaseStar(Vector& V1, const Vector& V2);
    void   PhaseRotStar(Vector& V1, const Vector& V2, const Array1D<int>& Re, const Array1D<int>& Im);

    void   Import(Vector& out, const Vector& in, const Vector& mesh, int deg_);
    void   Export(Vector& out, const Vector& mshint, const Vector& mshdeg, const Vector& in);
    void   pdMeshConvert(Vector& newprofile, Vector& newtangent, const Vector& oldprofile, const Vector& oldtangent);

    // computing the Jacobians, right-hand sides, characteristic matrices

    // continuation of solution

    void RHS(Vector& rhs, const Vector& par, const Vector& sol);
    void RHS_p(Vector& rhs, const Vector& par, const Vector& sol, int p);   // sol is currently not needed
    void RHS_x(SpMatrix& A, const Vector& par, const Vector& sol);          // sol is currently not needed

    // for stability computation
    void StabJac(StabMatrix& AB, const Vector& par);
    // this computes its own matrix structures, because it is nowhere else needed: kkSI, eeSI, rrSI, ddSI, etc.
    // However, the variables will be contained within the NColloc class

    // continuation of bifurcations -> characteristic matrices

    void CharJac_x(SpMatrix& A, const Vector& par, double Z);
    void CharJac_x(SpMatrix& A, const Vector& par, double ZRe, double ZIm);

    void CharJac_x_p(Vector& V, const Vector& par, const Array3D<double>& phiData, double Z, int p);
    void CharJac_x_p(Vector& V, const Vector& par,
      const Array3D<double>& phiDataRe, const Array3D<double>& phiDataIm, double ZRe, double ZIm, int p);

    void CharJac_x_x(SpMatrix& A, const Vector& par, const Array3D<double>& phiData, double Z);
    void CharJac_x_x(SpMatrix& A, const Vector& par,
      const Array3D<double>& phiDataRe, const Array3D<double>& phiDataIm, double Re, double Im);

    void CharJac_x_z(Vector& V, const Vector& par, const Vector& phi,
      const Array3D<double>& phiDataRe, const Array3D<double>& phiDataIm, double Re, double Im);
    // we do need Re, Im

    // for the autonomous FOLD bifurcation
    void CharJac_mB(SpMatrix& A, const Vector& par, double Z);
    void CharJac_mB_p(Vector& V,   const Vector& par, const Array3D<double>& phiData, double Z, int alpha);
    void CharJac_mB_x(SpMatrix& A, const Vector& par, const Array3D<double>& phiData, double Z);

    // autonomous trivial eigenvector
    void CharJac_MSHphi(Vector& V, const Vector& par, const Array3D<double>& solMSHData);
    template<bool trans> void CharJac_MSHphi_x(Vector& V, const Vector& par, const Array3D<double>& solMSHData, const Vector& phi);
    void CharJac_MSHphi_p(Vector& V, const Vector& par, const Array3D<double>& solMSHData, int alpha);

    // supplementary
    inline int Ndim() const
    {
      return ndim;
    }
    inline int Npar() const
    {
      return npar;
    }
    inline int Ntau() const
    {
      return ntau;
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
      P_ASSERT_X(msh.Size() == nint + 1, "Error in NColloc::setMesh : Bad dimensions\n");
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

  private:

    // helper functions
    void meshAdapt_internal( Vector& newmesh, const Vector& profile );

    // rINT and rDEG is included in idx. Only rDIM is necessary
    inline int& WRIDX(SpMatrix& A, int idx, int rDIM, int cTAU, int cDEG, int cDIM)
    {
      return A.WrLi(ndim + ndim*(idx) + rDIM, cDIM + ndim*(cDEG - dd(cTAU, idx) + rr(cTAU, idx)*(ndeg + 1)));
    }

    inline double& WRDAT(SpMatrix& A, int idx, int rDIM, int cTAU, int cDEG, int cDIM)
    {
      return A.WrLx(ndim + ndim*(idx) + rDIM, cDIM + ndim*(cDEG - dd(cTAU, idx) + rr(cTAU, idx)*(ndeg + 1)));
    }

    // for CharJac_x
    inline int& WRIDXCPLX(SpMatrix& A, int idx, int rDIM, int rIM, int cTAU, int cDEG, int cDIM, int cIM)
    {
      return A.WrLi(rIM + 2*(ndim + ndim*(idx) + rDIM), cIM + 2*(cDIM + ndim*(cDEG - dd(cTAU, idx) + rr(cTAU, idx)*(ndeg + 1))));
    }

    inline double& WRDATCPLX(SpMatrix& A, int idx, int rDIM, int rIM, int cTAU, int cDEG, int cDIM, int cIM)
    {
      return A.WrLx(rIM + 2*(ndim + ndim*(idx) + rDIM), cIM + 2*(cDIM + ndim*(cDEG - dd(cTAU, idx) + rr(cTAU, idx)*(ndeg + 1))));
    }

    // for CharJac_x_x
    inline int& WRIDXCPLXM(SpMatrix& A, int idx, int rDIM, int rIM, int cTAU, int cDEG, int cDIM)
    {
      return A.WrLi(rIM + 2*(ndim + ndim*(idx) + rDIM), cDIM + ndim*(cDEG - dd(cTAU, idx) + rr(cTAU, idx)*(ndeg + 1)));
    }

    inline double& WRDATCPLXM(SpMatrix& A, int idx, int rDIM, int rIM, int cTAU, int cDEG, int cDIM)
    {
      return A.WrLx(rIM + 2*(ndim + ndim*(idx) + rDIM), cDIM + ndim*(cDEG - dd(cTAU, idx) + rr(cTAU, idx)*(ndeg + 1)));
    }

    inline int& WRIDXS(SpMatrix& A, int idx, int rDIM, int cTAU, int cDEG, int cDIM)
    {
      return A.WrLi(ndim + ndim*(idx) + rDIM, cDIM + ndim*(cDEG - ddS(cTAU, idx) + rrS(cTAU, idx)*(ndeg + 1)));
    }

    double& WRDATS(SpMatrix& A, int idx, int rDIM, int cTAU, int cDEG, int cDIM)
    {
      return A.WrLx(ndim + ndim*(idx) + rDIM, cDIM + ndim*(cDEG - ddS(cTAU, idx) + rrS(cTAU, idx)*(ndeg + 1)));
    }

    inline int& WRIDXI(SpMatrix& A, int idx, int rDIM, int cTAU, int cDEG, int cDIM)
    {
      return A.WrLi(ndim + ndim*(idx) + rDIM, cDIM + ndim*(cDEG - ddI(cTAU, idx) + rrI(cTAU, idx)*(ndeg + 1)));
    }

    inline double& WRDATI(SpMatrix& A, int idx, int rDIM, int cTAU, int cDEG, int cDIM)
    {
      return A.WrLx(ndim + ndim*(idx) + rDIM, cDIM + ndim*(cDEG - ddI(cTAU, idx) + rrI(cTAU, idx)*(ndeg + 1)));
    }

    // the equations
    System* sys;

    const int ndim;
    const int npar;
    const int ntau;

    const int nint;
    const int ndeg;

    const int nmat;

    Vector mesh;
    Vector time;

    // these store the structure of the sparse matrix NO SEPARATION
    Array2D<int> kk;   // dim(NTAU+1,NDEG*NINT) : which delay to which interval
    Array2D<int> ee;   // dim(NTAU+1,NDEG*NINT) : the ordering of kk
    Array2D<int> dd;   // dim(NTAU+1,NDEG*NINT) : how many neighbouring intervals up to its index
    Array2D<int> rr;   // dim(NTAU+1,NDEG*NINT) : how many non overlapping intervals up to its index

    // these store the structure of the sparse matrix WITH SEPARATION
    Array2D<int> kkS;  // same as kk, but not folded back
    Array2D<int> eeS;  //  ...
    Array2D<int> rrS;  //  ...
    Array2D<int> ddS;  //  ...
    Array2D<int> mmS;  // which matrix are we in
    Array2D<int> szS;  // szI(mmI( . ,idx),idx) size of the line within the matrix mmS

    // For the stability matrices
    Array2D<int> kkI;
    Array2D<int> eeI;
    Array2D<int> rrI;
    Array2D<int> ddI;
    Array2D<int> mmI;
    Array2D<int> szI;

    // it stores all the collocation data
    Array3D<double> tt;       // interpolation at the collocation points
    Array1D<double> timeMSH;  // the representation points
    Array3D<double> ttMSH;    // interpolation at the representation points
    Array2D<int>    kkMSH;

    // matrix for integration
    Matrix metric;
    // integration with derivatives (for phase conditions)
    Matrix metricPhase;

    // internal use for the initialization
    Vector col;
    Vector out;
    Vector meshINT;
    Array1D< Array1D<double> > lgr;

    // for rhs, and derivatives
    // Matrix fx, dfx, t_dfx, dummy
    Array3D<double>  solData;
    Array2D<double>  p_tau;
    Array2D<double>  p_dtau;
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

#ifdef DEBUG
    // DIAGNOSTICS
    int count_RHS;
    int count_RHS_p;
    int count_RHS_x;
    int count_StabJac;
    int count_CharJac_x;
    int count_CharJac_x_p;
    int count_CharJac_x_x;
    int count_CharJac_x_z;
    int count_CharJac_mB;
    int count_CharJac_mB_p;
    int count_CharJac_mB_x;
    int count_CharJac_MSHphi;
    int count_CharJac_MSHphi_p;
    int count_CharJac_MSHphi_x;
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
      std::cout << "count_RHS " << count_RHS << "\n";
      std::cout << "count_RHS_p " << count_RHS_p << "\n";
      std::cout << "count_RHS_x " << count_RHS_x << "\n";
      std::cout << "count_StabJac " << count_StabJac << "\n";
      std::cout << "count_CharJac_x " << count_CharJac_x << "\n";
      std::cout << "count_CharJac_x_p " << count_CharJac_x_p << "\n";
      std::cout << "count_CharJac_x_x " << count_CharJac_x_x << "\n";
      std::cout << "count_CharJac_x_z " << count_CharJac_x_z << "\n";
      std::cout << "count_CharJac_mB " << count_CharJac_mB << "\n";
      std::cout << "count_CharJac_mB_p " << count_CharJac_mB_p << "\n";
      std::cout << "count_CharJac_mB_x " << count_CharJac_mB_x << "\n";
      std::cout << "count_CharJac_MSHphi " << count_CharJac_MSHphi << "\n";
      std::cout << "count_CharJac_MSHphi_p " << count_CharJac_MSHphi_p << "\n";
      std::cout << "count_CharJac_MSHphi_x " << count_CharJac_MSHphi_x << "\n";
    }
#endif // DEBUG
};

#define NDIM ndim
#define NPAR npar
#define NTAU ntau
#define NINT nint
#define NDEG ndeg
#define NMAT nmat

template<bool trans> void NColloc::CharJac_MSHphi_x(Vector& V, const Vector& par, const Array3D<double>& solMSHData, const Vector& phi)
{
#ifdef DEBUG
  count_CharJac_MSHphi_x++;
#endif
  V.Clear(); // it is not cleared otherwise!!!!

  for (int k = 0; k < NTAU; k++)
  {
    int nx = 1, vx = k, np = 0, vp;
    sys->p_deri(p_dfxMSH, timeMSH, solMSHData, par, nx, &vx, np, &vp, p_dummy);
    for (int idx = 0; idx < NINT*NDEG + 1; idx++) // i: interval; j: which collocation point
    {
      for (int l = 0; l < NDEG + 1; l++) // degree
      {
        for (int p = 0; p < NDIM; p++)   // row
        {
          for (int q = 0; q < NDIM; q++)   //column
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
