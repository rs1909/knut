// ------------------------------------------------------------------------- //
//
// This is part of PDDE-CONT
// Copyright (c) 2006 by Robert Szalai
//
// For license, see the file COPYING in the root directory of the package
//
// ------------------------------------------------------------------------- //

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif

#include "compthread.h"
#include "pderror.h"
#include "point.h"
#include "system.h"
#include "torpoint.h"

#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <string>
#include <cmath>

#include <QErrorMessage>

void MThread::run()
{
  try
  {
    System sys(params->getSysName());
    if (sys.ndim() == 0) P_MESSAGE("zerodimensions");

    Vector par(params->getNPar() + ParEnd);

    //-----------------------------------------------------------------------------------------------------------
    //
    // Initialize the system definition for the three cases of
    // a) refinement b) branch switching c) continuation
    //
    //-----------------------------------------------------------------------------------------------------------

    Array1D<Eqn> eqn; // used for continuation
    Array1D<Var> var; // used for continuation
    Array1D<Eqn> eqn_refine;
    Array1D<Var> var_refine;
    Array1D<Eqn> eqn_start;
    Array1D<Var> var_start;
    Eqn          testFN;
    int trivial = params->toEqnVar(sys, eqn, var, eqn_refine, var_refine, eqn_start, var_start, testFN);
    const int npar = params->getNPar();
    Point* pt_ptr;
    try
    {
      pt_ptr = new Point(sys, eqn_refine, var_refine, params->getNInt(), params->getNDeg(), params->getNMul(), params->getNMat());
    }
    catch (pddeException ex)
    {
      emit exceptionOccured(ex);
      return;
    }
    Point& pt = *pt_ptr;

    pt.setContIter(params->getNItC());
    pt.setRefIter(params->getNItR());
    pt.setKernIter(params->getNItK());
    pt.setRefEps(params->getEpsR());
    pt.setContEps(params->getEpsC());
    pt.setKernEps(params->getEpsK());
    pt.setCont(params->getCp() - VarPAR0);

    // setting the symmetric components
    Array1D<int> sre(params->getNSym()), sim(params->getNSym());
    for (int i = 0; i < sre.Size(); ++i)
    {
      sre(i) = params->getSymRe(i);
      sim(i) = params->getSymIm(i);
    }
    if (params->getNSym() != 0) pt.setSym(sre, sim);

    std::ostringstream screenout;
    screenout << std::scientific;
    screenout.precision(6);

    // load the initial guess
    if (params->getLabel() != 0)
    {
      mat4Data istr(params->getInputFile());
      pt.BinaryRead(istr, params->getLabel());
    }

    pt.Refine(screenout);
    if (testFN != EqnNone)
    {
      screenout << "\n--- Finding the bifurcation point (TF) ---\n";
      pt.Reset(eqn_start, var_start);
      pt.setCont(params->getCp() - VarPAR0);
      pt.StartTF(testFN, screenout);   // it only computes the characteristic multiplier refines the solution
    }
#ifdef DEBUG
    std::cout << screenout << "\n";
#endif
    emit printToScreen(screenout.str());
    screenout.str("");
    // start the continuation!
    if (params->getBranchSW() != TFTRSwitch)
    {
      mat4Data out(params->getOutputFile(),
                   params->getSteps(), sys.ndim(), sys.npar() + ParEnd,
                   params->getNInt(), params->getNDeg(), params->getNMul());

      screenout << "\n--- Starting the continuation ---\n";

      for (int j = 0; j < par.Size(); j++) par(j) = pt.getPar()(j);
      //
      parNamePrint(screenout, npar, params->getCp(), var);
      screenout << "\n";
      parValuePrint(screenout, par, params->getCp(), var, 0, pt.Norm(), 0, 0);
#ifdef DEBUG
      std::cout << screenout << "\n";
#endif
      emit printToScreen(screenout.str());
      screenout.str("");

      // making tangents
      if (params->getBranchSW() == TFPDSwitch)
      {
        screenout << "\nSwitching to the period two branch (TF).\n";
        pt.SwitchTFPD(params->getDsStart());
      }
      else if (params->getBranchSW() == TFBRSwitch)
      {
        screenout << "\nSwitching to the other branch (TF).\n";
        pt.SwitchTFLP(params->getDsStart());
      }
      else if (params->getBranchSW() == TFHBSwitch)
      {
        screenout << "\nSwitching to the periodic solution branch at the HOPF point (TF).\n";
        pt.SwitchTFHB(params->getDsStart());
      }
      else
      {
        screenout << "\nFinding the tangent.\n";
        pt.setCont(params->getCp() - VarPAR0);
        pt.Tangent();
      }
      pt.Reset(eqn, var);
      pt.setCont(params->getCp() - VarPAR0);
#ifdef DEBUG
      std::cout << screenout << "\n";
#endif
      emit printToScreen(screenout.str());
      screenout.str("");

      // if no stability computation then clear the previously computed multipliers
      if (!(params->getStab())) pt.clearStability();

      int ustab = 0, ustabprev = 0;
      double norm = 0.0;
      const int ithist = 5;
      Array1D<int> it(ithist);
      for (int i = 0; i < it.Size(); i++) it(i) = 3;
      int itpos = 0;
      double ds = params->getDs();
      for (int i = 0; i < params->getSteps(); i++)  // 35
      {
        if (stopFlag)
        {
          delete pt_ptr;
          return;
        }
        if (i % 24 == 0)
        {
          parNamePrint(screenout, npar, params->getCp(), var);
          screenout << "\n";
        }
        itpos = (itpos + 1) % ithist;
        //
        it(itpos) = pt.Continue(ds, (i == 0) && (params->getBranchSW() == TFHBSwitch));
        //
        if (params->getStab()) pt.Stability();
        ustabprev = ustab;
        if (trivial == 0) ustab = pt.UStab();
        else if (trivial == 1) ustab = pt.UStabAUT();
        else ustab = pt.UStabAUTRot();
        for (int j = 0; j < par.Size(); j++) par(j) = pt.getPar()(j);
        norm = pt.Norm();

        // console output
        parValuePrint(screenout, par, params->getCp(), var, i, norm, ustab, it(itpos));
        // adapt mesh if necessary
        if((params->getIad() != 0) && (((i+1) % params->getIad()) == 0))
        {
          const int itad = pt.Refine(std::cout,true);
          const int ittan = pt.Tangent(true);
          screenout << " " << itad << " " << ittan;
        }
        if (i != 0  && ustab != ustabprev)
        {
          PtType bif = SolTF;
          if (trivial == 0) bif = pt.testBif();
          else if (trivial == 1) bif = pt.testBifAUT();
          else pt.testBifAUTRot();
          switch (bif)
          {
            case BifTFLP:
            case BifTFAUTLP:
              screenout << "  LP";
              break;
            case BifTFPD:
            case BifTFAUTPD:
              screenout << "  PD";
              break;
            case BifTFNS:
            case BifTFAUTNS:
              screenout << "  NS";
              break;
            default:
              screenout << "  ??";
              break;
          }
        }
#ifdef DEBUG
        std::cout << screenout << "\n";
#endif
        emit printToScreen(screenout.str());
        screenout.str("");

        // file output
        pt.BinaryWrite(out, i);

        int itc = it(itpos);
        if ((itc > 3) && (fabs(ds) / 1.414 > params->getDsMin()) && (fabs(ds) / 1.414 < params->getDsMax())) ds /= 1.414;
        if ((itc > 5) && (fabs(ds) / 2.0 > params->getDsMin()) && (fabs(ds) / 2.0 < params->getDsMax())) ds /= 2.0;
        bool decr = true;
        for (int l = 0; l < it.Size(); l++) if (it(l) > 3) decr = false;
        if (decr && (fabs(ds)*1.414 > params->getDsMin()) && (fabs(ds)*1.414 < params->getDsMax())) ds *= 1.414;
        if ((itc >= params->getNItC()) && (fabs(ds) / 2.0 < params->getDsMin()))
        {
          P_MESSAGE("reached minimum stepsize (DSMIN)");
        }
        std::cout.flush();
      }
    }
    else
    {
      mat4Data out(params->getOutputFile(),
                   params->getSteps(), sys.ndim(), sys.npar() + ParEnd,
                   params->getNInt1(), params->getNInt2(), params->getNDeg1(), params->getNDeg2());

      std::cout << "ENTERING THE TORUS CODE!\n";
      // construct the initial torus
      double alpha;
      Vector Sol(pt.getSol().Size());
      Vector TRe(pt.getSol().Size()), TIm(pt.getSol().Size());

      // making the mesh for the conversion
      Vector meshint(params->getNInt1() + 1), meshdeg(params->getNDeg1() + 1);
      for (int i = 0; i < meshint.Size(); i++) meshint(i) = (double)i / (params->getNInt1());
      for (int i = 0; i < meshdeg.Size(); i++) meshdeg(i) = (double)i / (params->getNDeg1());

      // getting the sol and tangents
      pt.SwitchTRSol(Sol, meshint, meshdeg);
      pt.SwitchTFTRTan(TRe, TIm, alpha, meshint, meshdeg);

      // getting the parameters
      for (int j = 0; j < par.Size(); j++) par(j) = pt.getPar()(j);

      // destroy point, construct PointTR
      delete pt_ptr;
      pt_ptr = 0;
      PointTR pttr(sys, eqn, var, params->getNDeg1(), params->getNDeg2(), params->getNInt1(), params->getNInt2());

      // construct the solution tangent from the eigenvectors
      // these next three functions could be only one
      pttr.ImportSol(Sol);
      pttr.ImportTan(TRe, TIm, alpha);
      pttr.Start(params->getDsStart());
      pttr.setPar(par);

      pttr.setRho(alpha / (2.0*M_PI));
      pttr.setCont(params->getCp() - VarPAR0);

      double ds = params->getDs();
      std::ostringstream fdata, fidx;
      for (int i = 0; i < params->getSteps(); i++)
      {
        if (stopFlag)
        {
          return;
        }
        pttr.Continue(ds, false);

        // write out the results
        for (int j = 0; j < npar; j++) par(j) = pttr.getPar()(j);
        screenout << i;
        for (int j = 0; j < npar; j++) screenout << "\t" << par(j);
#ifdef DEBUG
        std::cout << screenout << "\n";
#endif
        emit printToScreen(screenout.str());
        screenout.str("");
        pttr.WriteBinary(out, i);
      }
    }
    // **********************************************************************************************************
    delete pt_ptr;
  }
  catch (pddeException ex)
  {
    emit exceptionOccured(ex);
    return;
  }
}
