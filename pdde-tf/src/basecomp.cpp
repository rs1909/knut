#include "basecomp.h"
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

void BaseComp::run(const char* branchFile)
{
  System sys(params->getSysName());
  if (sys.ndim() == 0) P_MESSAGE("Number of dimensions are set to zero.");

  Vector par(sys.npar() + ParEnd);
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

  Array1D<Eqn> eqn; // used for continuation
  Array1D<Var> var; // used for continuation
  Array1D<Eqn> eqn_refine;
  Array1D<Var> var_refine;
  Array1D<Eqn> eqn_start;
  Array1D<Var> var_start;
  Eqn          testFN;


  //-----------------------------------------------------------------------------------------------------------
  //
  // END of initialization
  //
  //-----------------------------------------------------------------------------------------------------------


  // just a block to contain pt, which eats too much memory
  try
  {
    params->toEqnVar(sys, eqn, var, eqn_refine, var_refine, eqn_start, var_start, testFN);
    const int npar = sys.npar();

//     for( int i=0; i<eqn.Size(); i++ ) std::cout<<EqnToStr( eqn(i) )<<", ";
//     std::cout<<'\n';
//     for( int i=0; i<var.Size(); i++ ) std::cout<<VarToStr( var(i) )<<", ";
//     std::cout<<'\n';
    std::ostringstream screenout;
    screenout << std::scientific;
    screenout.precision(6);

    // If it is a periodic solution OR switching from per. solotion to torus
    // Note that eqn_refine is _NEVER_ set to torus so we cannot rely on that.
    Point* pt_ptr = 0;
    if (eqn_start(0) != EqnTORSol)
    {
      pt_ptr = new Point(sys, eqn_refine, var_refine, params->getNInt(), params->getNDeg(), params->getNMul(), params->getNMat());
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

      // load the initial guess
      if (params->getLabel() != 0)
      {
        mat4Data istr(params->getInputFile());
        pt.BinaryRead(istr, params->getLabel());
      }
  
      pt.Refine(screenout);
      print(screenout);
      if (testFN != EqnNone)
      {
        screenout << "\n--- Finding the bifurcation point (TF) ---\n";
        pt.Reset(eqn_start, var_start);
        pt.setCont(params->getCp() - VarPAR0);
        pt.StartTF(testFN, screenout);   // it only computes the characteristic multiplier refines the solution
        print(screenout);
      }
    }
    // start the continuation!
    if (eqn(0) != EqnTORSol)
    {
      Point& pt = *pt_ptr;
      mat4Data out(params->getOutputFile(),
                   params->getSteps(), sys.ndim(), sys.npar() + ParEnd,
                   params->getNInt(), params->getNDeg(), params->getNMul());

      screenout << "\n--- Starting the continuation ---\n";

      for (int j = 0; j < par.Size(); j++) par(j) = pt.getPar()(j);
      //
      parNamePrint(screenout, npar, params->getCp(), var);
      screenout << "\n";
      parValuePrint(screenout, par, params->getCp(), var, 0, pt.Norm(), 0, 0);
      screenout << "\n";
      print(screenout);

      // making tangents
      switch (params->getBranchSW())
      {
        case TFPDSwitch:
          screenout << "\nSwitching to the period two branch (TF).\n";
          pt.SwitchTFPD(params->getDsStart());
          break;
        case TFHBSwitch:
          screenout << "\nSwitching to the periodic solution branch at the HOPF point (TF).\n";
          pt.SwitchTFHB(params->getDsStart());
          break;
        case TFBRSwitch:
        case TFBRAUTSwitch:
        case TFBRAUTROTSwitch:
          screenout << "\nSwitching to the other branch (TF).\n";
          pt.SwitchTFLP(static_cast<BranchSW>(params->getBranchSW()), params->getDsStart());
          break;
        default:
          screenout << "\nFinding the tangent.\n";
          pt.setCont(params->getCp() - VarPAR0);
          pt.Tangent();
          break;
      }
      print(screenout);
      pt.Reset(eqn, var);
      pt.setCont(params->getCp() - VarPAR0);

      // if no stability computation then clear the previously computed multipliers
      if (!(params->getStab())) pt.clearStability();

      screenout << '\n';
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
        print(screenout);
        itpos = (itpos + 1) % ithist;
        //
        it(itpos) = pt.Continue(ds, (i == 0) && (params->getBranchSW() != NOSwitch));
        if (stopFlag)
        {
          delete pt_ptr;
          return;
        }
        //
        if (params->getStab()) pt.Stability();
        ustabprev = ustab;
        ustab = pt.UStab();
        for (int j = 0; j < par.Size(); j++) par(j) = pt.getPar()(j);
        norm = pt.Norm();
        // console output
        parValuePrint(screenout, par, params->getCp(), var, i, norm, ustab, it(itpos));
        print(screenout);
        // adapt mesh if necessary
        if((params->getIad() != 0) && (((i+1) % params->getIad()) == 0))
        {
          const int itad = pt.Refine(screenout,true);
          const int ittan = pt.Tangent(true);
          screenout << " " << itad << " " << ittan;
          print(screenout);
        }
        if (i != 0  && ustab != ustabprev)
        {
          PtType bif = SolTF;
          bif = pt.testBif();
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
        screenout << "\n";
        print(screenout);

        // file output
        pt.BinaryWrite(out, i);

        // branch output
        if (branchFile)
        {
          for (int j = 0; j < par.Size(); j++) ff << par(j) << "\t";
          ff << "\t" << norm << "\t" << pt.NormMX() << "\t" << ustab << "\n";
          ff.flush();
        }
        int itc = it(itpos);
        if ((itc > 3) && (fabs(ds) / 1.414 > params->getDsMin()) && (fabs(ds) / 1.414 < params->getDsMax())) ds /= 1.414;
        if ((itc > 5) && (fabs(ds) / 2.0 > params->getDsMin()) && (fabs(ds) / 2.0 < params->getDsMax())) ds /= 2.0;
        bool decr = true;
        for (int l = 0; l < it.Size(); l++) if (it(l) > 3) decr = false;
        if (decr && (fabs(ds)*1.414 > params->getDsMin()) && (fabs(ds)*1.414 < params->getDsMax())) ds *= 1.414;
        if ((itc >= params->getNItC()) && (fabs(ds) / 2.0 < params->getDsMin()))
        {
          P_MESSAGE("No convergence. The minimum arclength step size (DSMIN) has been reached.");
        }
        // stop continuation if CP has reached the bounds
        if (par(params->getCp() - VarPAR0) < params->getCpMin()) break;
        if (par(params->getCp() - VarPAR0) > params->getCpMax()) break;
      }
    }
    else
    {
      mat4Data out(params->getOutputFile(),
                   params->getSteps(), sys.ndim(), sys.npar() + ParEnd,
                   params->getNInt1(), params->getNInt2(), params->getNDeg1(), params->getNDeg2());

      screenout << "ENTERING THE TORUS CODE!\n";
      print(screenout);

      // construct the initial torus
      double alpha;
      Vector Sol;
      Vector TRe, TIm;
      if (params->getBranchSW() == TFTRSwitch)
      {
        Point& pt = *pt_ptr;
        Sol.Init(pt.getSol().Size());
        TRe.Init(pt.getSol().Size());
        TIm.Init(pt.getSol().Size());

        // making the mesh for the conversion
        Vector meshint(params->getNInt1() + 1), meshdeg(params->getNDeg1() + 1);
        for (int i = 0; i < meshint.Size(); i++) meshint(i) = (double)i / (params->getNInt1());
        for (int i = 0; i < meshdeg.Size(); i++) meshdeg(i) = (double)i / (params->getNDeg1());

        // getting the sol and tangents
        screenout << "\nSwitching to the torus.\n"; print(screenout);
        pt.SwitchTRSol(Sol, meshint, meshdeg);
        pt.SwitchTFTRTan(TRe, TIm, alpha, meshint, meshdeg);

        // getting the parameters
        for (int j = 0; j < par.Size(); j++) par(j) = pt.getPar()(j);

        // destroy point, construct PointTR
        delete pt_ptr;
        pt_ptr = 0;
      }
      PointTR pttr(sys, eqn, var, params->getNDeg1(), params->getNDeg2(), params->getNInt1(), params->getNInt2());
      // construct the solution tangent from the eigenvectors
      // these next three functions could be only one
      if (params->getBranchSW() == TFTRSwitch)
      {
        pttr.ImportSol(Sol);
        pttr.ImportTan(TRe, TIm, alpha);
        pttr.Start(params->getDsStart());
        pttr.setPar(par);
        pttr.setRho(alpha / (2.0*M_PI));
        pttr.setCont(params->getCp() - VarPAR0);
      }
      else if (params->getBranchSW() == NOSwitch)
      {
        if (params->getLabel() != 0)
        {
          mat4Data istr(params->getInputFile());
          pttr.ReadBinary(istr, params->getLabel());
          screenout << "\nFinding the tangent.\n"; print(screenout);
          pttr.setCont(params->getCp() - VarPAR0);
          pttr.Tangent();
        } else
        {
          P_MESSAGE("A torus cannot be started from scratch.");
        }
      } else
      {
        P_MESSAGE("Invalid barnch switch.");
      }

      double ds = params->getDs();
      for (int i = 0; i < params->getSteps(); i++)
      {
        if (stopFlag)
        {
          return;
        }
        // same as for periodic orbits
        int it = pttr.Continue(ds, (i == 0) && (params->getBranchSW() != NOSwitch));
        for (int j = 0; j < par.Size(); j++) par(j) = pttr.getPar()(j);
        double norm = pttr.Norm();
        if (i % 24 == 0)
        {
          parNamePrint(screenout, npar, params->getCp(), var);
          screenout << "\n";
        }
        print(screenout);
        // console output
        parValuePrint(screenout, par, params->getCp(), var, i, norm, 0, it);
        screenout << "\n";
        print(screenout);
        if (branchFile)
        {
          for (int j = 0; j < par.Size(); j++) ff << par(j) << "\t";
          ff << norm << "\n";
          ff.flush();
        }
        pttr.WriteBinary(out, i);
      }
    }
    // **********************************************************************************************************
    delete pt_ptr;
  }
  catch (pddeException ex)
  {
    raiseException(ex);
    return;
  }
}
