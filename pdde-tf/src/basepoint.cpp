#include "basepoint.h"
#include "system.h"

// specified in the constants file
#define REFEPS    (1E-5)
#define NREFITER  (20)
#define CONTEPS   (1E-5)
#define NCONTITER (20)
#define KERNEPS   (1E-10)
#define NKERNITER (20)

BasePoint::BasePoint(System& sys_, const Array1D<Eqn>& eqn_, const Array1D<Var>& var_, 
          const int solsize, const int nz_jac_) :
    var(var_), eqn(eqn_), varMap(var_.Size()), varMapCont(var_.Size() + 1),
    sol(solsize), par(sys_.npar() + ParEnd),
    solNu(solsize), parNu(sys_.npar() + ParEnd)
{
  dim1   = solsize;
  nz_jac = nz_jac_;
  RefEps   = REFEPS;
  RefIter  = NREFITER;
  ContEps  = CONTEPS;
  ContIter = NCONTITER;
  KernEps  = KERNEPS;
  KernIter = NKERNITER;

  par(sys_.npar() + ParAngle) = 0.0;
  par(sys_.npar() + ParRot) = 0.0;
  // DO NOT CALL Construct! It will be called from its children!
}

BasePoint::~BasePoint()
{
  // DO NOT CALL Destruct! It will be called from its children!
}

// public
// remember that this ereases the state variables except sol, qq, par xxDot
void BasePoint::Reset(const Array1D<Eqn>& eqn_, const Array1D<Var>& var_)

{
  HyperVector* xxDot_temp = 0;
  if (xxDot) xxDot_temp = new HyperVector(*xxDot);
  Destruct();
  eqn.Init(eqn_.Size());
  eqn = eqn_;
  var.Init(var_.Size());
  var = var_;
  varMap.Init(var_.Size());
  varMapCont.Init(var_.Size() + 1);
  Construct();
  xxDot->getV1() = xxDot_temp->getV1();
  for (int i = 0; i < std::min<int>(xxDot_temp->getV3().Size(), xxDot->getV3().Size()); ++i)
    xxDot->getV3()(i) = xxDot_temp->getV3()(i);
  delete xxDot_temp;
}

void BasePoint::Construct()
{
  P_ERROR_X1((eqn.Size() != 0) && (var.Size() != 0) && (eqn.Size() == var.Size()), "Number of equations and variables do not agree.");
  dim3 = eqn.Size() - 1;

  for (int i = 1; i < var.Size(); i++)
  {
    P_ERROR_X4((var(i) - VarPAR0 >= 0) && (var(i) - VarPAR0 < par.Size()), "Non-existing parameter P", var(i) - VarPAR0, " at position ", i);
    varMap(i) = var(i) - VarPAR0;
  }
  for (int i = 0; i < var.Size(); i++) varMapCont(i) = varMap(i);

  xxDot   = new HyperVector(dim1, 0, dim3 + 1);
  xxDotNu = new HyperVector(dim1, 0, dim3 + 1);

  xx      = new HyperVector(dim1, 0, dim3 + 1);

  rhs     = new HyperVector(dim1, 0, dim3 + 1);

  jac     = new HyperMatrix(dim1, 0, dim3 + 1, nz_jac);
}

// private
void BasePoint::Destruct()
{
  delete jac;

  delete rhs;
  delete xx;

  delete xxDot;
  delete xxDotNu;
}

// this will adapt the mesh whenever it is specified
int BasePoint::Refine(std::ostream& out, bool adapt)
{
  // here solNu is the previous solution
  if (adapt)
  {
    basecolloc->meshAdapt(solNu, sol, xxDotNu->getV1(), xxDot->getV1());
    sol = solNu;
    xxDot->getV1() = xxDotNu->getV1();
  }
  solNu = sol;
  parNu = par;

  xx->getV3().Clear();

  if(!adapt) out << "IT\tERR\t\tSOLnorm\t\tDIFFnorm\n";

  int it = 0;
  double Xnorm, Dnorm;
  do
  {
    basecolloc->Init(sol, par);

    if (!adapt)
    {
      Jacobian(*jac, *rhs, parNu, par, solNu, sol, varMap, 0.0, false);
      jac->Solve(*xx, *rhs, dim3);
      Update(*xx);
    } else
    {
      Jacobian(*jac, *rhs, parNu, par, solNu, sol, varMapCont, 0.0, true);
      jac->Solve(*xx, *rhs, dim3+1);
      AdaptUpdate(*xx);
    }
    // computing norms to determine convergence
    Xnorm = sqrt(basecolloc->Integrate(sol, sol));
    Dnorm = sqrt(basecolloc->Integrate(xx->getV1(), xx->getV1()) + (xx->getV3()) * (xx->getV3()));
    if(!adapt)
    {
      out << " " << it << "\t" << Dnorm / (1.0 + Xnorm) << "\t" << Xnorm << "\t" << Dnorm << '\n';
      out.flush();
    }
  }
  while ((Dnorm / (1.0 + Xnorm) >= RefEps) && (it++ < RefIter));
  if (it >= RefIter) std::cout << "Warning: refinement did not converge. "
    << "CritNorm: " << Dnorm / (1.0 + Xnorm) << " SolNorm: " << Xnorm << " DiffNorm: " << Dnorm << '\n';

  return it;
}

int BasePoint::Tangent(bool adapt)
{
  double norm;

  basecolloc->Init(sol, par);

  if (!adapt)
  {
    // setting up a random tangent
    xxDot->getV1().Rand();
    xxDot->getV3().Rand();
    norm = sqrt(basecolloc->Integrate(xxDot->getV1(), xxDot->getV1()) + (xxDot->getV3()) * (xxDot->getV3()));
    xxDot->getV1() /= norm;
    xxDot->getV3() /= norm;
  }
  // az RHS-t feleslegesen szamolja ki && the first qq should be qq0
  Jacobian(*jac, *rhs, par, par, sol, sol, varMapCont, 0.0, true);

  double diffnorm = 1.0;
  int it = 0;
  do
  {
    jac->Multiply<false>(*rhs, *xxDot, dim3 + 1);
    rhs->getV3()(dim3) -= 1.0;
    jac->Solve(*xx, *rhs);
    xxDot->getV1() -= xx->getV1();
    xxDot->getV3() -= xx->getV3();
    diffnorm = sqrt(basecolloc->Integrate(xx->getV1(), xx->getV1()) + (xx->getV3()) * (xx->getV3()));
    norm = sqrt(basecolloc->Integrate(xxDot->getV1(), xxDot->getV1()) + (xxDot->getV3()) * (xxDot->getV3()));
    xxDot->getV1() /= norm;
    xxDot->getV3() /= norm;
    // putting back the tangent...
    if (dim1 != 0) basecolloc->Star(jac->getA31(dim3), xxDot->getV1());
    for (int i = 0; i < dim3 + 1; i++) jac->getA33()(dim3, i) = xxDot->getV3()(i);
  }
  while ((++it < KernIter) && (diffnorm > KernEps));
  if (diffnorm > KernEps) std::cout << "Point::Tangent: warning: No convergence in finding the singular vector. Residual = " << diffnorm << ", steps " << it << "\n";
  if (!adapt && (xxDot->getV3()(dim3) < 0.0))
  {
    xxDot->getV1() *= -1.0;
    xxDot->getV3() *= -1.0;
  }

  return it;
}

int BasePoint::Continue(double ds, bool jacstep)
{
  double Xnorm, Dnorm, Rnorm, Tnorm;

  parNu = par;
  for (int i = 0; i < solNu.Size(); i++)  solNu(i)           = sol(i)           + ds * xxDot->getV1()(i);
  for (int i = 1; i < varMapCont.Size(); i++) parNu(varMapCont(i)) = par(varMapCont(i)) + ds * xxDot->getV3()(i - 1);
  xxDotNu->getV1() = xxDot->getV1();
  xxDotNu->getV3() = xxDot->getV3();

  int  it = 0;
  bool conv;
  do
  {
    basecolloc->Init(solNu, parNu);

    Jacobian(*jac, *rhs, par, parNu, sol, solNu, varMapCont, 0.0, true);

    jac->Solve(*xx, *rhs);

    ContUpdate(*xx);

    Rnorm = sqrt(basecolloc->Integrate(rhs->getV1(), rhs->getV1()) + (rhs->getV3()) * (rhs->getV3()));
    Xnorm = sqrt(basecolloc->Integrate(solNu, solNu));
    Dnorm = sqrt(basecolloc->Integrate(xx->getV1(), xx->getV1()) + (xx->getV3()) * (xx->getV3()));
    conv = (Dnorm / (1.0 + Xnorm) >= ContEps) || (Rnorm / (1.0 + Xnorm) >= ContEps);

#ifdef DEBUG
    std::cout << "Dnorm: " << Dnorm << " Rnorm: " << Rnorm << " Xnorm: " << Xnorm << "\n";
    std::cout.flush();
#endif /*DEBUG*/
    // updating the tangent
    if (!jacstep)
    {
      jac->Multiply<false>(*rhs, *xxDotNu, dim3 + 1);
      rhs->getV3()(dim3) -= 1.0;
      jac->Solve(*xx, *rhs);
      xxDotNu->getV1() -= xx->getV1();
      xxDotNu->getV3() -= xx->getV3();
      Tnorm = sqrt(basecolloc->Integrate(xxDotNu->getV1(), xxDotNu->getV1()) + (xxDotNu->getV3()) * (xxDotNu->getV3()));
      xxDotNu->getV1() /= Tnorm;
      xxDotNu->getV3() /= Tnorm;
    }
    // end updating tangent
  }
  while (conv /*&& (Dnorm/(1.0+Xnorm) < 1.0)*/ && (++it < ContIter));
  if (!conv)
  {
#ifdef DEBUG
    /// checking the tangent and the secant
    double Pnorm = sqrt(xxDotNu->getV3()(dim3) * xxDotNu->getV3()(dim3));
    double Xnorm = sqrt(basecolloc->Integrate(xxDotNu->getV1(), xxDotNu->getV1())), Onorm = sqrt((xxDotNu->getV3()) * (xxDotNu->getV3()));
    std::cout << "Cnorm: " << Tnorm << "\nDot Pnorm: " << Pnorm << " Xnorm: " << Xnorm << " Onorm: " << Onorm;
    for (int i = 1; i < varMap.Size(); i++) std::cout << " O" << varMap(i) << ": " << xxDotNu->getV3()(i - 1);
    std::cout << '\n';

    xx->getV1() = solNu;
    xx->getV1() -= sol;
    for (int i = 1; i < varMapCont.Size(); i++) xx->getV3()(i - 1) = parNu(varMapCont(i)) - par(varMapCont(i));

    Pnorm = sqrt(xx->getV3()(dim3) * xx->getV3()(dim3)) / ds;
    Xnorm = sqrt(basecolloc->Integrate(xx->getV1(), xx->getV1())) / ds;
    Onorm = 0;
    for (int i = 0; i < dim3 + 1; i++) Onorm += (xx->getV3()(i)) * (xx->getV3()(i));
    Onorm = sqrt(Onorm) / ds;
    std::cout << "Dif Pnorm: " << Pnorm << " Xnorm: " << Xnorm << " Onorm: " << Onorm;
    for (int i = 1; i < varMap.Size(); i++) std::cout << " O" << varMap(i) << ": " << xx->getV3()(i - 1) / ds;
    std::cout << '\n';
    /// END OF CHECKING
#endif

    // copying back the solution
    sol = solNu;
    par = parNu;
    xxDot->getV1() = xxDotNu->getV1();
    xxDot->getV3() = xxDotNu->getV3();
  }
  else
  {
    std::cout << "\n\n\n ------------------- NO CONVERGENCE -------------------\n\n\n\n";
    // P_MESSAGE("");
  }

  return it;
}
