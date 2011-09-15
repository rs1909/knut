#include "basecomp.h"
#include "knerror.h"
#include "system.h"
#include "point.h"
#include "odepoint.h"
#include "torpoint.h"

#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <string>
#include <cmath>

#define ITSTEPS 5
#define ITLIM1 3
#define ITLIM2 5

void KNAbstractContinuation::run(const char* branchFile)
{
  KNSystem *sys = 0;
  try
  {
    sys = new KNSystem(params->getSysName(), params->getNDeri());
  }
  catch (KNException ex)
  {
    raiseException(ex);
    delete sys;
    return;
  }

  
  if (sys->ndim() == 0) P_MESSAGE1("Number of dimensions are set to zero.");
  params->initDimensions(sys);

  KNVector par(VarToIndex(VarEnd,sys->npar()));
  std::ofstream ff;
  if (branchFile)
  {
    ff.open(branchFile);
    ff << std::scientific;
    ff.precision(12);
  }

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

  KNAbstractPeriodicSolution* pt_ptr = 0;
  try
  {
    const bool needFN = params->toEqnVar(*sys, eqn, var, eqn_refine, var_refine, eqn_start, var_start, findangle);
    const int npar = sys->npar();
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
      if (eqn_start(0) == EqnSol) pt_ptr = new KNDdePeriodicSolution(this, *sys, eqn_refine, var_refine, params->getNInt(), params->getNDeg(), params->getNMul());
      else if (eqn_start(0) == EqnODESol) pt_ptr = new KNOdePeriodicSolution(this, *sys, eqn_refine, var_refine, params->getNInt(), params->getNDeg());
      else P_MESSAGE3("There is no such solution type: ", eqn_start(0), ".");
      KNAbstractPeriodicSolution& pt = *pt_ptr;
  
      pt.setContIter(params->getNItC());
      pt.setRefIter(params->getNItR());
      pt.setKernIter(params->getNItK());
      pt.setRefEps(params->getEpsR());
      pt.setContEps(params->getEpsC());
      pt.setKernEps(params->getEpsK());
      pt.setContCurvature(params->getCAngle());
      pt.setCont(params->getCp());
  
      // setting the symmetric components
      KNArray1D<int> sre(params->getSymReSize());
      KNArray1D<int> sim(params->getSymReSize());
      for (int i = 0; i < sre.size(); ++i)
      {
        sre(i) = params->getSymRe(i);
        sim(i) = params->getSymIm(i);
      }
      if (params->getSymReSize() != 0) pt.setSym(sre, sim);

      // load the initial guess
      if (params->getLabel() != 0)
      {
        KNDataFile istr(params->getInputFile());
        if (params->getFromType() != BifNone)
        {
          int pos = istr.findType(params->getFromType(), params->getLabel() );
          std::cout << " getfrom Type " << params->getFromType() << " lab " << params->getLabel() << " pos " << pos << "\n";
          if (pos > 0) pt.BinaryRead(istr, pos);
          else P_MESSAGE1("No such point in the input");
        } else
        {
          pt.BinaryRead(istr, params->getLabel()-1);
        }
      }

      screenout   << "\n---      Refine supplied solution      ---\n";
      printStream();
      pt.refine();
      if (needFN)
      {
        screenout << "\n--- Finding the bifurcation point (TF) ---\n";
        pt.reset(eqn_start, var_start);
        pt.setCont(params->getCp());
        if (findangle) pt.findAngle();   // it only computes the angle from the test functional
        pt.refine();
      }
    }
    // start the continuation!
    if (eqn(0) != EqnTORSol)
    {
      KNAbstractPeriodicSolution& pt = *pt_ptr;
      setData(new KNDataFile(params->getOutputFile(), params->getParNames(),
                           params->getSteps(), sys->ndim(), VarToIndex(VarEnd,sys->npar()),
                           params->getNInt(), params->getNDeg(), params->getNMul()));

      screenout   << "\n---     Starting the continuation      ---\n";

      for (int j = 0; j < par.size(); j++) par(j) = pt.getPar()(j);
      //
      parNamePrint(screenout, npar, params->getCp(), var, parNames);
      screenout << "\n";
      parValuePrint(screenout, par, params->getCp(), var, -1, BifNone, pt.norm(), 0, 0);
      screenout << "\n";
      printStream();

      // making tangents
      switch (params->getBranchSW())
      {
        case TFPDSwitch:
          screenout << "\nSwitching to the period two branch (TF).\n";
          printStream();
          pt.SwitchTFPD(params->getDsStart());
          break;
        case TFHBSwitch:
          screenout << "\nSwitching to the periodic solution branch at the HOPF point (TF).\n";
          printStream();
          pt.SwitchTFHB(params->getDsStart());
          break;
        case TFBRSwitch:
        case TFBRAUTSwitch:
        case TFBRAUTROTSwitch:
          screenout << "\nSwitching to the other branch (TF).\n";
          printStream();
          pt.SwitchTFLP(static_cast<BranchSW>(params->getBranchSW()), params->getDsStart());
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
      if (!(params->getStab())) pt.clearStability();

      screenout << '\n';
      printStream();
      storeCursor();
      int ustab = 0, ustabprev = 0;
      double norm = 0.0;
      const int ithist = ITSTEPS;
      int printedln = 0;
      double ds = params->getDs();
      for (int i = 0; i < params->getSteps(); i++)  // 35
      {
        //
        double angle;
        const int itc = pt.nextStep(ds, angle, (i == 0) && (params->getBranchSW() != NOSwitch));
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
          delete pt_ptr;
          delete sys;
          return;
        }
        //
        pt.storeMultiplier();
        if (params->getStab()) pt.Stability(false); // we just stepped: no init
        ustabprev = ustab;
        ustab = pt.UStab();
        for (int j = 0; j < par.size(); j++) par(j) = pt.getPar()(j);
        norm = pt.norm();
        // console output
        const bool endpoint = i == params->getSteps()-1;
        const bool stabchange = (i != 0) && (ustab != ustabprev);
        const bool toprint = (((i+1) % params->getNPr()) == 0)||(i == 0);
        // stability output
        BifType bif = BifNone;
        if (endpoint) bif = BifEndPoint;
        if (stabchange) bif = pt.testBif();
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
            screenout << "--";
          } else if ( (fabs(ds)*dsmul > params->getDsMax()) )
          {
            screenout << "++";
          } else
          {
            screenout << " " << std::showpos << static_cast<int>(round(2*log(dsmul)/log(2))) << std::noshowpos;
          }
          printStream();
          ++printedln;
        }
        // adapt mesh if necessary
        if ((params->getIad() != 0) && (((i+1) % params->getIad()) == 0))
        {
          const int itad = pt.refine(true);
          const int ittan = pt.tangent(true);
          if (toprint || (i != 0  && ustab != ustabprev))
          {
            screenout << " " << itad << " " << ittan;
            printStream();
          }
        }
        if (toprint || stabchange)
        {
          screenout << "\n";
          printStream();
          storeCursor();
        }
        // file output
        pt.BinaryWrite(data(), bif, i);
        dataUpdated();

        // branch output
        if (branchFile)
        {
          for (int j = 0; j < par.size(); j++) ff << par(j) << "\t";
          ff << "\t" << norm << "\t" << pt.NormMX() << "\t" << ustab << "\n";
          ff.flush();
        }
        if ( (fabs(ds)*dsmul >= params->getDsMin()) && (fabs(ds)*dsmul <= params->getDsMax()) ) ds *= dsmul;
        else if ((itc >= params->getNItC()) && (fabs(ds)*dsmul < params->getDsMin()))
        {
          parValuePrint(screenout, par, params->getCp(), var, i, BifNoConvergence, norm, ustab, itc);
          screenout << '\n';
          printStream();
          break;
        }
        // stop continuation if CP has reached the bounds
        if ((par(params->getCp() - VarPAR0) < params->getCpMin())||
            (par(params->getCp() - VarPAR0) > params->getCpMax()))
        {
          parValuePrint(screenout, par, params->getCp(), var, i, BifMax, norm, ustab, itc);
          screenout << '\n';
          printStream();
          break;
        }
      }
      deleteData();
    }
    else
    {
      setData(new KNDataFile(params->getOutputFile(), params->getParNames(),
                           params->getSteps(), sys->ndim(), VarToIndex(VarEnd,sys->npar()),
                           params->getNInt1(), params->getNInt2(), params->getNDeg1(), params->getNDeg2()));

      screenout << "ENTERING THE TORUS CODE!\n";
      printStream();

      // construct the initial torus
      double alpha;
      KNVector Sol;
      KNVector TRe, TIm;
      if (params->getBranchSW() == TFTRSwitch)
      {
        if (dynamic_cast<KNDdePeriodicSolution*>(pt_ptr) == 0) P_MESSAGE1("Cannot switch to torus because it is not a delay equation.");
        KNDdePeriodicSolution& pt = *dynamic_cast<KNDdePeriodicSolution*>(pt_ptr);
        Sol.init(pt.getSol().size());
        TRe.init(pt.getSol().size());
        TIm.init(pt.getSol().size());

        // making the mesh for the conversion
        KNVector meshint(params->getNInt1() + 1), meshdeg(params->getNDeg1() + 1);
        for (int i = 0; i < meshint.size(); i++) meshint(i) = (double)i / (params->getNInt1());
        for (int i = 0; i < meshdeg.size(); i++) meshdeg(i) = (double)i / (params->getNDeg1());

        // getting the sol and tangents
        screenout << "\nSwitching to the torus.\n";
        printStream();
        pt.SwitchTRSol(Sol, meshint, meshdeg);
        pt.SwitchTFTRTan(TRe, TIm, alpha, meshint, meshdeg);

        // getting the parameters
        for (int j = 0; j < par.size(); j++) par(j) = pt.getPar()(j);

        // destroy point, construct KNDdeTorusSolution
        delete pt_ptr;
        pt_ptr = 0;
      }
      KNDdeTorusSolution pttr(this, *sys, eqn, var, params->getNDeg1(), params->getNDeg2(), params->getNInt1(), params->getNInt2());
      // construct the solution tangent from the eigenvectors
      // these next three functions could be only one
      if (params->getBranchSW() == TFTRSwitch)
      {
        pttr.importSolution(Sol);
        pttr.importTangent(TRe, TIm, alpha);
        pttr.startingPoint(params->getDsStart());
        pttr.setPar(par);
        pttr.setRho(alpha / (2.0*M_PI));
        pttr.setCont(params->getCp());
      }
      else if (params->getBranchSW() == NOSwitch)
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
        P_MESSAGE1("Invalid barnch switch.");
      }

      double ds = params->getDs();
      for (int i = 0; i < params->getSteps(); i++)
      {
        if (stopFlag)
        {
          deleteData();
          delete sys;
          return;
        }
        // same as for periodic orbits
        double angle;
        int it = pttr.nextStep(ds, angle, (i == 0) && (params->getBranchSW() != NOSwitch));
        for (int j = 0; j < par.size(); j++) par(j) = pttr.getPar()(j);
        double norm = pttr.norm();
        if (i % 24 == 0)
        {
          parNamePrint(screenout, npar, params->getCp(), var, parNames);
          screenout << "\n";
        }
        printStream();
        // console output
        parValuePrint(screenout, par, params->getCp(), var, i, BifNone, norm, 0, it);
        screenout << "\n";
        printStream();
        if (branchFile)
        {
          for (int j = 0; j < par.size(); j++) ff << par(j) << "\t";
          ff << norm << "\n";
          ff.flush();
        }
        pttr.savePoint(data(), i);
        dataUpdated();
      }
      deleteData();
    }
    // **********************************************************************************************************
    delete pt_ptr;
    delete sys;
  }
  catch (KNException ex)
  {
    delete pt_ptr;
    delete sys;
    deleteData();
    raiseException(ex);
  }
}
