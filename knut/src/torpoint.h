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

class System;
class mat4Data;

#include "matrix.h"
#include "hypermatrix.h"
#include "torcolloc.h"
#include "pointtype.h"
#include "basepoint.h"

class PointTR : public BasePoint
{
  public:
    PointTR(System& sys_, Array1D<Eqn>& eqn_, Array1D<Var>& var_, int ndeg1_, int ndeg2_, int nint1_, int nint2_);
    ~PointTR();

    inline void    setRho(double rho)
    {
      par(colloc->Sys().npar() + ParRot) = rho;
      std::cout << "RHO: " << rho << "\n";
      std::cout << "par(0): " << par(0) << "\n";
    }

    // making starting solution from a periodic solution
    inline void    ImportSol(Vector& Sol) // "Sol" is already a torus with zero diameter
    {
      colloc->ImportSol(xx->getV1(), Sol);
    }
    inline void    ImportTan(Vector& Re, Vector& Im, double alpha)
    {
      colloc->ImportTan(xxDot->getV1(), Re, Im, alpha);
      xxDot->getV3().Clear();
//       p1Dot = 0.0; // should be xxDot->getV3(varMap.size()) ???
    }
    inline void    Start(double ds)
    {
      for (int i = 0; i < xx->getV1().size(); i++) xx->getV1()(i) += ds * (xxDot->getV1()(i));
    }
//     inline int     npar()
//     {
//       return colloc->Sys().npar();
//     }
    void           ReadBinary(mat4Data& data, int n);
    void           WriteBinary(mat4Data& data, int n);

  private:

//     void ToJacVar(Array1D<Vector*>& A13, Array1D<int>& JacVar, bool cont);
//     void JacobianFixed(Vector& sol, Vector& presol, Vector& par, Vector& prepar, bool cont);
//     void JacobianVarying(Array1D<Vector*> A13, Array1D<int>& JacVar,
//                          Vector& sol, Vector& presol, Vector& par, Vector& prepar, bool cont, double ds);
    void Jacobian
    (
      HyperMatrix& AA, HyperVector& RHS, // output
      Vector& parPrev, Vector& par,      // parameters
      Vector& solPrev, Vector& sol,      // solution
      Array1D<int>&    varMap,           // contains the variables. If cont => contains the P0 too.
      double ds, bool cont               // ds stepsize, cont: true if continuation
    );

    inline const   Vector& getMesh1() const
    {
      return colloc->getMesh1();
    }
    inline const   Vector& getMesh2() const
    {
      return colloc->getMesh2();
    }

    CollocTR* colloc;
};

#endif
