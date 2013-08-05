// ------------------------------------------------------------------------- //
//
// This is part of KNUT
// Copyright (c) 2013 by Robert Szalai
//
// For license, see the file COPYING in the package's root directory
//
// ------------------------------------------------------------------------- //

#ifndef STPOINT_H
#define STPOINT_H


#include "basepoint.h"

class KNExprSystem;
class KNVector;
class KNBlockVector;
class KNSteadyStateJacobian;
template< class T > class KNArray1D;

class KNSteadyStateSolution : public KNAbstractPoint
{
  public:

    KNSteadyStateSolution(KNAbstractContinuation* cnt, KNExprSystem& sys, 
      KNArray1D<Eqn>& eqn_, KNArray1D<Var>& var_);
    ~KNSteadyStateSolution();
    
    virtual void BinaryRead(KNDataFile& data, size_t n);
    virtual void BinaryWrite(KNDataFile& data, BifType bif, size_t n);

    virtual void SwitchTFLP(BranchSW , double ) { P_MESSAGE1("Not implemented"); }     // switches branch with testFunct
    virtual void SwitchTFPD(double ) { P_MESSAGE1("Not implemented"); }                    // switches branch with testFunct
    virtual void SwitchTFHB(double ) { P_MESSAGE1("Not implemented"); } // switches branch with testFunct
  private:
    virtual void jacobian(
      KNSparseBlockMatrix& AA, KNBlockVector& RHS, // output
      KNVector& parPrev, KNVector& par,      // parameters
      KNVector& solPrev, KNVector& sol,      // solution
      KNArray1D<size_t>& varMap,           // contains the variables. If cont => contains the P0 too.
      double ds, bool cont               // ds stepsize, cont: true if continuation
    );
    virtual void postProcess() {}
    
    KNSteadyStateJacobian* colloc;
    KNLuSparseMatrix     jacStab;
};

#endif
