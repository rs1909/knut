// ------------------------------------------------------------------------- //
//
// This is part of KNUT
// Copyright (c) 2002, 2003, 2004, 2005 by Robert Szalai
//
// For license, see the file COPYING in the package's root directory
//
// ------------------------------------------------------------------------- //

#ifndef TORPOINT_H
#define TORPOINT_H

class KNExprSystem;
class KNDataFile;

#include "matrix.h"
#include "hypermatrix.h"
#include "torcolloc.h"
#include "pointtype.h"
#include "basepoint.h"

class KNDdeTorusSolution : public KNAbstractPoint
{
  public:
    KNDdeTorusSolution(KNAbstractContinuation* cnt, KNExprSystem& sys_, 
      KNArray1D<Eqn>& eqn_, KNArray1D<Var>& var_, size_t ndeg1_, size_t ndeg2_, size_t nint1_, size_t nint2_);
    ~KNDdeTorusSolution() override;

    inline void    setRho(double rho)
    {
      par(VarToIndex(VarRot,colloc->system().npar())) = rho;
      std::cout << "RHO: " << rho << "\n";
      std::cout << "par(0): " << par(0) << "\n";
    }

    // making starting solution from a periodic solution
    inline void    importSolution(KNVector& Sol) // "Sol" is already a torus with zero diameter
    {
      colloc->importSolution(xx->getV1(), Sol);
    }
    inline void    importTangent(KNVector& Re, KNVector& Im, double alpha)
    {
      colloc->importTangent(xxDot->getV1(), Re, Im, alpha);
      xxDot->getV3().clear();
//       p1Dot = 0.0; // should be xxDot->getV3(varMap.size()) ???
    }
    inline void    startingPoint(double ds)
    {
      for (size_t i = 0; i < xx->getV1().size(); i++) xx->getV1()(i) += ds * (xxDot->getV1()(i));
    }
//     inline int     npar()
//     {
//       return colloc->system().npar();
//     }
    void           loadPoint(const KNAbstractData& data, size_t n);
    void           savePoint(KNAbstractData& data, size_t n);
    void           BinaryRead(const KNAbstractData& data, size_t n) override { loadPoint(data, n); }
    void           BinaryWrite(KNAbstractData& data, BifType bif, size_t n) override { savePoint(data, n); }


  private:

//     void ToJacVar(KNArray1D<KNVector*>& A13, KNArray1D<int>& JacVar, bool cont);
//     void JacobianFixed(KNVector& sol, KNVector& presol, KNVector& par, KNVector& prepar, bool cont);
//     void JacobianVarying(KNArray1D<KNVector*> A13, KNArray1D<int>& JacVar,
//                          KNVector& sol, KNVector& presol, KNVector& par, KNVector& prepar, bool cont, double ds);
    void jacobian(
      KNSparseBlockMatrix& AA, KNBlockVector& RHS, // output
      KNVector& parPrev, KNVector& par,      // parameters
      KNVector& solPrev, KNVector& sol,      // solution
      KNArray1D<size_t>& varMap,           // contains the variables. If cont => contains the P0 too.
      double ds, bool cont               // ds stepsize, cont: true if continuation
    ) override;
    void postProcess() override {}

    inline const   KNVector& getMesh1() const
    {
      return colloc->getMesh1();
    }
    inline const   KNVector& getMesh2() const
    {
      return colloc->getMesh2();
    }

    KNDdeTorusCollocation* colloc;
};

#endif
