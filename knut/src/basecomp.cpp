#include "basecomp.h"
#include "knerror.h"
#include "point.h"
#include "odepoint.h"
#include "stpoint.h"
#include "torpoint.h"
#include "exprsystem.h"

#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <string>
#include <cmath>

#define ITSTEPS 5
#define ITLIM1 3
#define ITLIM2 5

// void KNAbstractContinuation::run(const char* branchFile)
// {
//   KNSystem *sys = 0;
//   try
//   {
//     sys = new KNSystem(params->getSysName());
//   }
//   catch (KNException ex)
//   {
//     delete sys; sys = 0;
//     raiseException(ex);
//     return;
//   }
//   if (sys->ndim() == 0) P_MESSAGE1("Number of dimensions are set to zero.");
//   params->initDimensions(sys);
//   run(sys, branchFile);
// }

void KNAbstractContinuation::run(KNExprSystem* sys, KNConstants* params,
                                 const KNAbstractData* inputData)
{
  P_ERROR_X1 (sys != nullptr, "Invalid system");
  if (sys->ndim() == 0) P_MESSAGE1("Number of dimensions are set to zero.");
  params->initDimensions(sys);

  setStopFlag(false);
  KNVector par(VarToIndex(VarEnd,sys->npar()));

  //-----------------------------------------------------------------------------------------------------------
  //
  // Initialize the system definition for the three cases of
  // a) refinement b) branch switching c) continuation
  //
  //-----------------------------------------------------------------------------------------------------------

  KNArray1D<Eqn> eqn; // used for continuation
  KNArray1D<Var> var; // used for continuation
  KNArray1D<Eqn> eqn_refine;
  KNArray1D<Var> var_refine;
  KNArray1D<Eqn> eqn_start;
  KNArray1D<Var> var_start;
  bool         findangle;

  //-----------------------------------------------------------------------------------------------------------
  //
  // END of initialization
  //
  //-----------------------------------------------------------------------------------------------------------

  KNAbstractPoint* pt_ptr = nullptr;
  KNAbstractPeriodicSolution* pt_per_ptr = nullptr;
  try
  {
    const bool needFN = params->toEqnVar(*sys, eqn, var, eqn_refine, var_refine, eqn_start, var_start, findangle);
    const size_t npar = sys->npar();
    std::vector<std::string> parNames = params->getParNames();

    screenout << std::scientific;
    screenout.precision(6);

//    screenout << "- REFINE -\n";
//    for( int i=0; i<eqn_refine.size(); i++ ) screenout<<EqnToStr( eqn_refine(i) )<<", ";
//    screenout << '\n';
//    for( int i=0; i<var_refine.size(); i++ ) screenout<<VarToStr( var_refine(i) )<<", ";
//    screenout << "\n- START -\n";
//    for( int i=0; i<eqn_start.size(); i++ ) screenout<<EqnToStr( eqn_start(i) )<<", ";
//    screenout << '\n';
//    for( int i=0; i<var_start.size(); i++ ) screenout<<VarToStr( var_start(i) )<<", ";
//    screenout << "\n- CONTINUE -\n";
//    for( int i=0; i<eqn.size(); i++ ) screenout<<EqnToStr( eqn(i) )<<", ";
//    screenout << '\n';
//    for( int i=0; i<var.size(); i++ ) screenout<<VarToStr( var(i) )<<", ";
//    screenout << "\n--\n";
//    printStream();

    // If it is a periodic solution OR switching from per. solotion to torus
    // Note that eqn_refine is _NEVER_ set to torus so we cannot rely on that.
    if (eqn_start(0) != EqnTORSol)
    {
      if (eqn_start(0) == EqnSol) pt_per_ptr = new KNDdePeriodicSolution(this, *sys, eqn_refine, var_refine, params->getNInt(), params->getNDeg(), params->getNMul());
      else if (eqn_start(0) == EqnODESol) pt_per_ptr = new KNOdePeriodicSolution(this, *sys, eqn_refine, var_refine, params->getNInt(), params->getNDeg());
      else if (eqn_start(0) == EqnSteady) pt_ptr = new KNSteadyStateSolution(this, *sys, eqn_refine, var_refine);
      else P_MESSAGE3("There is no such solution type: ", eqn_start(0), ".");
      if (pt_per_ptr) pt_ptr = pt_per_ptr;
      KNAbstractPoint& pt = *pt_ptr;

      pt.setContIter(params->getNItC());
      pt.setRefIter(params->getNItR());
      pt.setKernIter(params->getNItK());
      pt.setRefEps(params->getEpsR());
      pt.setContEps(params->getEpsC());
      pt.setKernEps(params->getEpsK());
      pt.setContCurvature(params->getCAngle());
      pt.setCont(params->getCp());

      // setting the symmetric components
      KNArray1D<size_t> sre(params->getSymReSize());
      KNArray1D<size_t> sim(params->getSymReSize());
      for (size_t i = 0; i < sre.size(); ++i)
      {
        sre(i) = params->getSymRe(i);
        sim(i) = params->getSymIm(i);
      }
      if (params->getSymReSize() != 0 && pt_per_ptr) pt_per_ptr->setSym(sre, sim); // !!

      // load the initial guess
      if (params->getLabel() != 0)
      {
        if (inputData != nullptr)
        {
          if (params->getFromType() != BifType::BifNone)
          {
            size_t pos = inputData->findType(params->getFromType(), params->getLabel() );
//             std::cout << " getfrom Type " << params->getFromType() << " lab " << params->getLabel() << " pos " << pos << "\n";
            if (pos < inputData->getNCols()) pt.BinaryRead(*inputData, pos); // !!
            else P_MESSAGE1("No such point in the input");
          } else
          {
            pt.BinaryRead(*inputData, params->getLabel()-1); // !!
          }
        } else P_MESSAGE1("Missing input data.");
      }

      if (needFN)
      {
        screenout << "\n--- Finding the bifurcation point (TF) ---\n";
        pt.reset(eqn_start, var_start);
        pt.setCont(params->getCp());
        if (findangle && pt_per_ptr) pt_per_ptr->findAngle(params->getSSM());   // it only computes the angle from the test functional // !!
        pt.refine();
      } else
      {
        screenout   << "\n---      Refine supplied solution      ---\n";
        printStream();
        pt.refine();
      }
    }
    // start the continuation!
    if (eqn(0) != EqnTORSol)
    {
      if (pt_ptr == nullptr) P_MESSAGE1("Internal error, cannot start the continuation.");
      KNAbstractPoint& pt = *pt_ptr;
      createData (params->getOutputFile(), pt_per_ptr != nullptr ? DataType::LC : DataType::ST, sys->ndim(), sys->npar(), params);

      screenout   << "\n---     Starting the continuation      ---\n";

      for (size_t j = 0; j < par.size(); j++) par(j) = pt.getPar()(j);
      //
      parNamePrint(screenout, npar, params->getCp(), var, parNames);
      screenout << "\n";
      parValuePrint(screenout, par, params->getCp(), var, 0, BifType::BifNone, pt.norm(), 0, 0);
      screenout << "\n";
      printStream();

      // making tangents
      switch (params->getBranchSW())
      {
        case BranchSW::TFPDSwitch:
          screenout << "\nSwitching to the period two branch (TF).\n";
          printStream();
          if (pt_per_ptr) pt_per_ptr->SwitchTFPD(params->getDsStart()); // !!
          break;
        case BranchSW::TFSSMSwitch:
            pt_per_ptr->setSSM();
             screenout << "\nIt is a subcentre manifold (TF).\n";
        case BranchSW::TFHBSwitch:
          screenout << "\nSwitching to the periodic solution branch at the HOPF point (TF).\n";
          printStream();
          if (pt_per_ptr) pt_per_ptr->SwitchTFHB(params->getDsStart()); // !!
          break;
        case BranchSW::TFBRSwitch:
        case BranchSW::TFBRAUTSwitch:
        case BranchSW::TFBRAUTROTSwitch:
          screenout << "\nSwitching to the other branch (TF). " << params->getBranchSW() << "\n";
          printStream();
          if (pt_per_ptr) pt_per_ptr->SwitchTFLP(static_cast<BranchSW>(params->getBranchSW()), params->getDsStart()); // !!
          break;
        default:
          screenout << "\nFinding the tangent.\n";
          printStream();
          pt.setCont(params->getCp());
          pt.tangent();
          break;
      }
      pt.reset(eqn, var);
      pt.setCont(params->getCp());

      // if no stability computation then clear the previously computed multipliers
      if (!(params->getStab()) && pt_per_ptr) pt_per_ptr->clearStability(); // !!

      screenout << '\n';
      printStream();
      storeCursor();
      size_t ustab = 0, ustabprev = 0;
      double norm = 0.0;
//       const size_t ithist = ITSTEPS;
      size_t printedln = 0;
      double ds = params->getDs();
      // setting trivial multipliers
      if( params->getBranchSW() == BranchSW::TFSSMSwitch )  pt_per_ptr->setSSM();
      for (size_t i = 0; i < params->getSteps(); i++)  // 35
      {
        //
        double angle;
        const size_t itc = pt.nextStep(ds, angle, ((i == 0) && (params->getBranchSW() != BranchSW::NOSwitch)) ? IterateTangent::no : IterateTangent::yes);
        // step size adaptation
        double dsmul1 = 1.0, dsmul2 = 1.0;
        if      (itc < 2) dsmul1 = 2.0;         // 2
        else if (itc < 3) dsmul1 = sqrt(2);     // 1
        else if (itc < 4) dsmul1 = 1.0;         // 0
        else if (itc < 5) dsmul1 = 1.0/sqrt(2); // -1
        else              dsmul1 = 1.0/2.0;     // -2
        const double ATOL = params->getCAngle();
        if      (angle < ATOL/4) dsmul2 = 2.0;
        else if (angle < ATOL/3) dsmul2 = sqrt(2);
        else if (angle < ATOL/2) dsmul2 = 1.0;
        else if (angle < ATOL)   dsmul2 = 1.0/sqrt(2);
        else                     dsmul2 = 1.0/2.0;
        const double dsmul = fmin(dsmul1,dsmul2);
        // whether we need to stop
        if (stopFlag)
        {
          deleteData();
          delete pt_ptr; pt_ptr = nullptr;
//          delete sys; sys = 0;
          return;
        }
        //
        if (pt_per_ptr)
        {
          pt_per_ptr->storeMultiplier(); // !!
          if (params->getStab()) pt_per_ptr->Stability(false); // we just stepped: no init // !!
          ustabprev = ustab;
          ustab = pt_per_ptr->UStab(); // !!
        }
        for (size_t j = 0; j < par.size(); j++) par(j) = pt.getPar()(j);
        norm = pt.norm();

        // console output
        const bool endpoint = i == params->getSteps()-1;
        const bool stabchange = (i != 0) && (ustab != ustabprev);
        const bool toprint = (((i+1) % params->getNPr()) == 0)||(i == 0);
        // stability output
        BifType bif = BifType::BifNone;
        if (endpoint) bif = BifType::BifEndPoint;
        if (stabchange && pt_per_ptr != nullptr) bif =  pt_per_ptr->testBif(); // !!

        if (toprint || stabchange || endpoint)
        {
          // ereases everyting from the last saved cursor
          clearLastLine();
          if (printedln % 24 == 0)
          {
            parNamePrint(screenout, npar, params->getCp(), var, parNames);
            screenout << "\n";
          }
          parValuePrint(screenout, par, params->getCp(), var, i, bif, norm, ustab, itc);
          if ( (fabs(ds)*dsmul < params->getDsMin()) )
          {
            screenout << "-0";
          } else if ( (fabs(ds)*dsmul > params->getDsMax()) )
          {
            screenout << "+0";
          } else
          {
            screenout << std::showpos << static_cast<int>(round(2*log(dsmul)/log(2))) << std::noshowpos;
          }
          printStream();
          ++printedln;
        }
        // adapt mesh if necessary
        if ((params->getIad() != 0) && (((i+1) % params->getIad()) == 0))
        {
          pt.adapt();
          const size_t itad = pt.nextStep(0.0, angle, IterateTangent::yes);
          if (itad >= params->getNItC()) break;
        }
        if (toprint || stabchange)
        {
          screenout << "\n";
          printStream();
          storeCursor();
        }
        // file output
        pt.BinaryWrite(data(), bif, i);
        if (toprint || endpoint) dataUpdated();

        if ( (fabs(ds)*dsmul >= params->getDsMin()) && (fabs(ds)*dsmul <= params->getDsMax()) ) ds *= dsmul;
        else if ((itc >= params->getNItC()) && (fabs(ds)*dsmul < params->getDsMin()))
        {
          parValuePrint(screenout, par, params->getCp(), var, i, BifType::BifNoConvergence, norm, ustab, itc);
          screenout << '\n';
          printStream();
          break;
        }
        // stop continuation if CP has reached the bounds
        if ((par(VarToIndex(params->getCp(),sys->npar())) < params->getCpMin())||
            (par(VarToIndex(params->getCp(),sys->npar())) > params->getCpMax()))
        {
          parValuePrint(screenout, par, params->getCp(), var, i, BifType::BifMax, norm, ustab, itc);
          screenout << '\n';
          printStream();
          break;
        }
      }
      deleteData();
    }
    else
    {
      createData (params->getOutputFile(), DataType::TR, sys->ndim(), sys->npar(), params);
//      setData(new KNDataFile(params->getOutputFile(), params->getParNames(),
//                           params->getSteps(), sys->ndim(), VarToIndex(VarEnd,sys->npar()),
//                           params->getNInt1(), params->getNInt2(), params->getNDeg1(), params->getNDeg2()));

      screenout << "ENTERING THE TORUS CODE!\n";
      printStream();

      // construct the initial torus
      double alpha;
      KNVector Sol;
      KNVector TRe, TIm;
      if (params->getBranchSW() == BranchSW::TFTRSwitch)
      {
        if (dynamic_cast<KNDdePeriodicSolution*>(pt_ptr) == nullptr) P_MESSAGE1("Cannot switch to torus because it is not a delay equation.");
        KNDdePeriodicSolution& pt = *dynamic_cast<KNDdePeriodicSolution*>(pt_ptr);
        Sol.init(pt.getSol().size());
        TRe.init(pt.getSol().size());
        TIm.init(pt.getSol().size());

        // making the mesh for the conversion
        KNVector meshint(params->getNInt1() + 1), meshdeg(params->getNDeg1() + 1);
        for (size_t i = 0; i < meshint.size(); i++) meshint(i) = (double)i / (params->getNInt1());
        for (size_t i = 0; i < meshdeg.size(); i++) meshdeg(i) = (double)i / (params->getNDeg1());

        // getting the sol and tangents
        screenout << "\nSwitching to the torus.\n";
        printStream();
        pt.SwitchTRSol(Sol, meshint, meshdeg);
        pt.SwitchTFTRTan(TRe, TIm, alpha, meshint, meshdeg);

        // getting the parameters
        for (size_t j = 0; j < par.size(); j++) par(j) = pt.getPar()(j);

        // destroy point, construct KNDdeTorusSolution
        delete pt_ptr;
        pt_ptr = nullptr;
      }
      KNDdeTorusSolution pttr(this, *sys, eqn, var, params->getNDeg1(), params->getNDeg2(), params->getNInt1(), params->getNInt2());
      // construct the solution tangent from the eigenvectors
      // these next three functions could be only one
      if (params->getBranchSW() == BranchSW::TFTRSwitch)
      {
        pttr.importSolution(Sol);
        pttr.importTangent(TRe, TIm, alpha);
        pttr.startingPoint(params->getDsStart());
        pttr.setPar(par);
        pttr.setRho(alpha / (2.0*M_PI));
        pttr.setCont(params->getCp());
      }
      else if (params->getBranchSW() == BranchSW::NOSwitch)
      {
        if (params->getLabel() != 0)
        {
          KNDataFile istr(params->getInputFile());
          pttr.loadPoint(istr, params->getLabel()-1);
          screenout << "\nFinding the tangent.\n";
          printStream();
          pttr.setCont(params->getCp());
          pttr.tangent();
        } else
        {
          P_MESSAGE1("A torus cannot be started from scratch.");
        }
      } else
      {
        P_MESSAGE1("Invalid branch switch.");
      }

      double ds = params->getDs();
      for (size_t i = 0; i < params->getSteps(); i++)
      {
        if (stopFlag)
        {
          deleteData();
//          delete sys; sys = 0;
          return;
        }
        // same as for periodic orbits
        double angle;
        size_t it = pttr.nextStep(ds, angle, ((i == 0) && (params->getBranchSW() != BranchSW::NOSwitch)) ? IterateTangent::no : IterateTangent::yes);
        for (size_t j = 0; j < par.size(); j++) par(j) = pttr.getPar()(j);
        double norm = pttr.norm();
        if (i % 24 == 0)
        {
          parNamePrint(screenout, npar, params->getCp(), var, parNames);
          screenout << "\n";
        }
        printStream();
        // console output
        parValuePrint(screenout, par, params->getCp(), var, i, BifType::BifNone, norm, 0, it);
        screenout << "\n";
        printStream();
        pttr.savePoint(data(), i);
        dataUpdated();
      }
      deleteData();
    }
    // **********************************************************************************************************
    delete pt_ptr; pt_ptr = nullptr;
//     delete sys; sys = 0;
  }
  catch (KNException ex)
  {
    delete pt_ptr; pt_ptr = nullptr;
//     delete sys; sys = 0;
    deleteData();
    raiseException(ex);
  }
}

// create a data file for different type of solutions
KNDataFile* KNAbstractContinuation::createDataStatic (const std::string& fileName, DataType t, size_t ndim, size_t npar, KNConstants* prms)
{
  if (t == DataType::ST) {
    return new KNDataFile(prms->getOutputFile(), prms->getParNames(),
                          prms->getSteps(), ndim, VarToIndex(VarEnd,npar),
                          0, 0, prms->getNMul());
  } else if (t == DataType::LC) {
    return new KNDataFile(prms->getOutputFile(), prms->getParNames(),
                          prms->getSteps(), ndim, VarToIndex(VarEnd,npar),
                          prms->getNInt(), prms->getNDeg(), prms->getNMul());
  } else if (t == DataType::TR) {
    return new KNDataFile(prms->getOutputFile(), prms->getParNames(),
                          prms->getSteps(), ndim, VarToIndex(VarEnd,npar),
                          prms->getNInt1(), prms->getNInt2(), prms->getNDeg1(), prms->getNDeg2());
  }
  return nullptr;
}
