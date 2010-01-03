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

  KNVector par(sys->npar() + ParEnd);
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

    std::ostringstream screenout;
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
//    print(screenout);

    // If it is a periodic solution OR switching from per. solotion to torus
    // Note that eqn_refine is _NEVER_ set to torus so we cannot rely on that.
    if (eqn_start(0) != EqnTORSol)
    {
      if (eqn_start(0) == EqnSol) pt_ptr = new KNDdePeriodicSolution(*sys, eqn_refine, var_refine, params->getNInt(), params->getNDeg(), params->getNMul(), params->getNMat());
      else if (eqn_start(0) == EqnODESol) pt_ptr = new KNOdePeriodicSolution(*sys, eqn_refine, var_refine, params->getNInt(), params->getNDeg());
      else P_MESSAGE3("There is no such solution type: ", eqn_start(0), ".");
      KNAbstractPeriodicSolution& pt = *pt_ptr;
  
      pt.setContIter(params->getNItC());
      pt.setRefIter(params->getNItR());
      pt.setKernIter(params->getNItK());
      pt.setRefEps(params->getEpsR());
      pt.setContEps(params->getEpsC());
      pt.setKernEps(params->getEpsK());
      pt.setCont(params->getCp() - VarPAR0);
  
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
        pt.BinaryRead(istr, params->getLabel()-1);
      }

      screenout   << "\n---      Refine supplied solution      ---\n";
      print(screenout);
      pt.refine(screenout);
      print(screenout);
      if (needFN)
      {
        screenout << "\n--- Finding the bifurcation point (TF) ---\n";
        pt.reset(eqn_start, var_start);
        pt.setCont(params->getCp() - VarPAR0);
        pt.StartTF(findangle, screenout);   // it only computes the characteristic multiplier refines the solution
        print(screenout);
      }
    }
    // start the continuation!
    if (eqn(0) != EqnTORSol)
    {
      KNAbstractPeriodicSolution& pt = *pt_ptr;
      setData(new KNDataFile(params->getOutputFile(), params->getParNames(),
                           params->getSteps(), sys->ndim(), sys->npar() + ParEnd,
                           params->getNInt(), params->getNDeg(), params->getNMul()));

      screenout   << "\n---     Starting the continuation      ---\n";

      for (int j = 0; j < par.size(); j++) par(j) = pt.getPar()(j);
      //
      parNamePrint(screenout, npar, params->getCp(), var, parNames);
      screenout << "\n";
      parValuePrint(screenout, par, params->getCp(), var, -1, BifNone, pt.norm(), 0, 0);
      screenout << "\n";
      print(screenout);

      // making tangents
      switch (params->getBranchSW())
      {
        case TFPDSwitch:
          screenout << "\nSwitching to the period two branch (TF).\n";
          print(screenout);
          pt.SwitchTFPD(params->getDsStart());
          break;
        case TFHBSwitch:
          screenout << "\nSwitching to the periodic solution branch at the HOPF point (TF).\n";
          print(screenout);
          pt.SwitchTFHB(params->getDsStart(), screenout);
          break;
        case TFBRSwitch:
        case TFBRAUTSwitch:
        case TFBRAUTROTSwitch:
          screenout << "\nSwitching to the other branch (TF).\n";
          print(screenout);
          pt.SwitchTFLP(static_cast<BranchSW>(params->getBranchSW()), params->getDsStart());
          break;
        default:
          screenout << "\nFinding the tangent.\n";
          pt.setCont(params->getCp() - VarPAR0);
          print(screenout);
          pt.tangent();
          break;
      }
      pt.reset(eqn, var);
      pt.setCont(params->getCp() - VarPAR0);

      // if no stability computation then clear the previously computed multipliers
      if (!(params->getStab())) pt.clearStability();

      screenout << '\n';
      int ustab = 0, ustabprev = 0;
      double norm = 0.0;
      const int ithist = ITSTEPS;
      KNArray1D<int> it(ithist);
      for (int i = 0; i < it.size(); i++) it(i) = 3;
      int itpos = 0;
      int printedln = 0;
      double ds = params->getDs();
      for (int i = 0; i < params->getSteps(); i++)  // 35
      {
        itpos = (itpos + 1) % ithist;
        //
        it(itpos) = pt.nextStep(ds, (i == 0) && (params->getBranchSW() != NOSwitch));
        if (stopFlag)
        {
          deleteData();
          delete pt_ptr;
          delete sys;
          return;
        }
        //
        if (params->getStab()) pt.Stability();
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
          if (printedln % 24 == 0)
          {
            parNamePrint(screenout, npar, params->getCp(), var, parNames);
            screenout << "\n";
          }
          parValuePrint(screenout, par, params->getCp(), var, i, bif, norm, ustab, it(itpos));
          print(screenout);
          ++printedln;
        }
        // adapt mesh if necessary
        if ((params->getIad() != 0) && (((i+1) % params->getIad()) == 0))
        {
          const int itad = pt.refine(screenout,true);
          const int ittan = pt.tangent(true);
          if (toprint || (i != 0  && ustab != ustabprev))
          {
            screenout << " " << itad << " " << ittan;
            print(screenout);
          }
        }
        if (toprint || stabchange)
        {
          screenout << "\n";
          print(screenout);
        }
        // file output
        pt.BinaryWrite(data(), i);
        dataUpdated();

        // branch output
        if (branchFile)
        {
          for (int j = 0; j < par.size(); j++) ff << par(j) << "\t";
          ff << "\t" << norm << "\t" << pt.NormMX() << "\t" << ustab << "\n";
          ff.flush();
        }
        // step size adaptation
        const int itc = it(itpos);
        if ((itc > ITLIM1) && 
            (fabs(ds) / 1.414 > fabs(params->getDsMin())) && 
            (fabs(ds) / 1.414 < fabs(params->getDsMax()))) ds /= 1.414;
        if ((itc > ITLIM2) && 
            (fabs(ds) / 2.0 > fabs(params->getDsMin())) && 
            (fabs(ds) / 2.0 < fabs(params->getDsMax()))) ds /= 2.0;
        bool decr = true;
        for (int l = 0; l < it.size(); l++) if (!(it(l) < ITLIM1)) decr = false;
        if ( decr && 
            (fabs(ds)*1.414 > fabs(params->getDsMin())) && 
            (fabs(ds)*1.414 < fabs(params->getDsMax()))) ds *= 1.414;
        if ((itc >= params->getNItC()) && (fabs(ds) / 2.0 < params->getDsMin()))
        {
          parValuePrint(screenout, par, params->getCp(), var, i, BifNoConvergence, norm, ustab, it(itpos));
          screenout << '\n';
          print(screenout);
          break;
        }
        // stop continuation if CP has reached the bounds
        if ((par(params->getCp() - VarPAR0) < params->getCpMin())||
            (par(params->getCp() - VarPAR0) > params->getCpMax()))
        {
          parValuePrint(screenout, par, params->getCp(), var, i, BifMax, norm, ustab, it(itpos));
          screenout << '\n';
          print(screenout);
          break;
        }
      }
      deleteData();
    }
    else
    {
      setData(new KNDataFile(params->getOutputFile(), params->getParNames(),
                           params->getSteps(), sys->ndim(), sys->npar() + ParEnd,
                           params->getNInt1(), params->getNInt2(), params->getNDeg1(), params->getNDeg2()));

      screenout << "ENTERING THE TORUS CODE!\n";
      print(screenout);

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
        screenout << "\nSwitching to the torus.\n"; print(screenout);
        pt.SwitchTRSol(Sol, meshint, meshdeg);
        pt.SwitchTFTRTan(TRe, TIm, alpha, meshint, meshdeg);

        // getting the parameters
        for (int j = 0; j < par.size(); j++) par(j) = pt.getPar()(j);

        // destroy point, construct KNDdeTorusSolution
        delete pt_ptr;
        pt_ptr = 0;
      }
      KNDdeTorusSolution pttr(*sys, eqn, var, params->getNDeg1(), params->getNDeg2(), params->getNInt1(), params->getNInt2());
      // construct the solution tangent from the eigenvectors
      // these next three functions could be only one
      if (params->getBranchSW() == TFTRSwitch)
      {
        pttr.importSolution(Sol);
        pttr.importTangent(TRe, TIm, alpha);
        pttr.startingPoint(params->getDsStart());
        pttr.setPar(par);
        pttr.setRho(alpha / (2.0*M_PI));
        pttr.setCont(params->getCp() - VarPAR0);
      }
      else if (params->getBranchSW() == NOSwitch)
      {
        if (params->getLabel() != 0)
        {
          KNDataFile istr(params->getInputFile());
          pttr.loadPoint(istr, params->getLabel()-1);
          screenout << "\nFinding the tangent.\n"; print(screenout);
          pttr.setCont(params->getCp() - VarPAR0);
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
        int it = pttr.nextStep(ds, (i == 0) && (params->getBranchSW() != NOSwitch));
        for (int j = 0; j < par.size(); j++) par(j) = pttr.getPar()(j);
        double norm = pttr.norm();
        if (i % 24 == 0)
        {
          parNamePrint(screenout, npar, params->getCp(), var, parNames);
          screenout << "\n";
        }
        print(screenout);
        // console output
        parValuePrint(screenout, par, params->getCp(), var, i, BifNone, norm, 0, it);
        screenout << "\n";
        print(screenout);
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
