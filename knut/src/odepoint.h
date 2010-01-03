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

class KNSystem;
class KNVector;
class KNBlockVector;
class KNOdeBvpCollocation;
template< class T > class KNArray1D;

class KNOdePeriodicSolution : public KNAbstractPeriodicSolution
{
  public:

    KNOdePeriodicSolution(KNSystem& sys, KNArray1D<Eqn>& eqn_, KNArray1D<Var>& var_, int nint, int ndeg);
    ~KNOdePeriodicSolution();
    
    virtual void Stability();
    virtual void SwitchTFLP(BranchSW type, double ds) { P_MESSAGE1("Not implemented"); }     // switches branch with testFunct
    virtual void SwitchTFPD(double ds) { P_MESSAGE1("Not implemented"); }                    // switches branch with testFunct
    virtual void SwitchTFHB(double ds, std::ostream& out) { P_MESSAGE1("Not implemented"); } // switches branch with testFunct
  private:
    virtual void jacobian(
      KNSparseBlockMatrix& AA, KNBlockVector& RHS, // output
      KNVector& parPrev, KNVector& par,      // parameters
      KNVector& solPrev, KNVector& sol,      // solution
      KNArray1D<int>&    varMap,           // contains the variables. If cont => contains the P0 too.
      double ds, bool cont               // ds stepsize, cont: true if continuation
    );
    
    KNOdeBvpCollocation* colloc;
    KNLuSparseMatrix     jacStab;
    KNMatrix     matrixInitialCondition;
    KNMatrix     matrixSolution;
    KNMatrix     monodromyMatrix;
};

#endif // ODEPOINT_H
