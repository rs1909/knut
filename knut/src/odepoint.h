// ------------------------------------------------------------------------- //
//
// This is part of KNUT
// Copyright (c) 2009 by Robert Szalai
//
// For license, see the file COPYING in the package's root directory
//
// ------------------------------------------------------------------------- //

#ifndef ODEPOINT_H
#define ODEPOINT_H

#include "pointtype.h"
#include "basepoint.h"

class System;
class Vector;
class HyperVector;
class ODEColloc;
template< class T > class Array1D;

class ODEPoint : public PerSolPoint
{
  public:

    ODEPoint(System& sys, Array1D<Eqn>& eqn_, Array1D<Var>& var_, int nint, int ndeg);
    ~ODEPoint();
    
    virtual void Stability();
    virtual void SwitchTFLP(BranchSW type, double ds) { P_MESSAGE1("Not implemented"); }     // switches branch with testFunct
    virtual void SwitchTFPD(double ds) { P_MESSAGE1("Not implemented"); }                    // switches branch with testFunct
    virtual void SwitchTFHB(double ds, std::ostream& out) { P_MESSAGE1("Not implemented"); } // switches branch with testFunct
  private:
    virtual void jacobian(
      HyperMatrix& AA, HyperVector& RHS, // output
      Vector& parPrev, Vector& par,      // parameters
      Vector& solPrev, Vector& sol,      // solution
      Array1D<int>&    varMap,           // contains the variables. If cont => contains the P0 too.
      double ds, bool cont               // ds stepsize, cont: true if continuation
    );
    
    ODEColloc* colloc;
    SpFact     jacStab;
    Matrix     matrixInitialCondition;
    Matrix     matrixSolution;
    Matrix     monodromyMatrix;
};

#endif // ODEPOINT_H
