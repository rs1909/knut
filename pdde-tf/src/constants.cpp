// ------------------------------------------------------------------------- //
//
// This is part of KNUT
// Copyright (c) 2006 by Robert Szalai
//
// For license, see the file COPYING in the root directory of the package
//
// ------------------------------------------------------------------------- //

#include "constants.h"
#include <iostream>
#include <fstream>
#include <sstream>
// For loading XML constant files
#include <mxml.h>
static const char *knut_whitespace_cb(mxml_node_t *node, int where)	
{
  const unsigned int max_levels = 12;
  const unsigned int tab = 1;
  static char spaces[tab*max_levels+2];
  const char *name = node->value.element.name;

  bool chc = false;
  if (node->child) if (node->child->child) chc = true;
  if (strncmp(name,"?xml",4) == 0)
  {
  	if (where == MXML_WS_AFTER_OPEN) return "\n";
  	else return 0;
  }
  if ((where == MXML_WS_BEFORE_OPEN)||(chc&&(where == MXML_WS_BEFORE_CLOSE)))
  {
    unsigned int level = 0;
    mxml_node_t *parent = node->parent;
    while (parent) { ++level; parent = parent->parent; }

	level -= 1;
    if (level > max_levels) level = max_levels;
    else if (level < 0) level = 0;
    for (unsigned int i = 0; i < level*tab; ++i ) spaces[i] = ' ';
    spaces[level*tab] = '\0';
    return spaces;
  }
  else if ((where == MXML_WS_AFTER_CLOSE)||(chc&&(where == MXML_WS_AFTER_OPEN)))
  {
  	return "\n";
  } else if ((node->child == 0)&&(where == MXML_WS_AFTER_OPEN))
  {
  	return "\n";
  }
  return (0);
}


inline bool NConstants::inputAssert(std::istream& is)
{
  if (is.rdstate() & (std::istream::eofbit | std::istream::failbit | std::istream::badbit))
  {
    switch (is.rdstate())
    {
      case std::istream::eofbit:
        P_MESSAGE("Unexpected end of file");
        return true;
        break;
      case std::istream::failbit:
        P_MESSAGE("Input failed\n");
        return true;
        break;
      case std::istream::badbit:
        P_MESSAGE("Bad input\n");
        return true;
        break;
      default:
        P_MESSAGE("Unexpected error");
        return true;
        break;
    }
  }
  return false;
}

static inline const char *getNodeText(mxml_node_t* nd)
{
  if (nd != 0)
  {
    if (nd->child != 0)
    {
      if (nd->child == nd->last_child)
      {
        if (nd->child->type == MXML_TEXT)
        {
          return nd->child->value.text.string; 
        } else std::cout << "XML wrong element type.";
      } else std::cout << "XML node has too many children.";
    } else return "";
  } else std::cout << "XML node could not be found.";
  return "";
}

static inline const char getNodeChar(mxml_node_t* nd)
{
  const char* str = getNodeText(nd);
  return str[0];
}

static inline int getNodeInteger(mxml_node_t* nd)
{
  const char* str = getNodeText(nd);
  return strtol(str, 0, 10);
}

static inline double getNodeReal(mxml_node_t* nd)
{
  std::istringstream str(getNodeText(nd));
  double val;
  str >> val;
  return val;  
}

static inline const char* c2s(char *buf, char c)
{
  buf[0] = c;
  buf[1] = '\0';
  return buf;
}

void NConstants::loadXmlFile(const std::string &fileName)
{
  FILE *fp;
  mxml_node_t *tree;

  fp = fopen(fileName.c_str(), "r");
  tree = mxmlLoadFile(0, fp, MXML_TEXT_CALLBACK);
  fclose(fp);

  mxml_node_t* nd = 0, *nd_a = 0, *nd_b = 0;
  mxml_node_t* root_nd = mxmlFindElement(tree, tree, "knut", 0, 0, MXML_DESCEND_FIRST);
  
  nd = mxmlFindElement(root_nd, root_nd, "input", 0, 0, MXML_DESCEND_FIRST);
  setInputFileText(getNodeText(nd));
  
  nd = mxmlFindElement(root_nd, root_nd, "output", 0, 0, MXML_DESCEND_FIRST);
  setOutputFileText(getNodeText(nd));
  
  nd = mxmlFindElement(root_nd, root_nd, "sysname", 0, 0, MXML_DESCEND_FIRST);
  setSysNameText(getNodeText(nd));

  nd = mxmlFindElement(root_nd, root_nd, "label", 0, 0, MXML_DESCEND_FIRST);
  setLabel(getNodeInteger(nd));
  
  nd = mxmlFindElement(root_nd, root_nd, "pointtype", 0, 0, MXML_DESCEND_FIRST);
  setPointType(static_cast<PtType>(getNodeInteger(nd)));
    
  nd_a = mxmlFindElement(root_nd, root_nd, "cptype", 0, 0, MXML_DESCEND_FIRST);
  nd_b = mxmlFindElement(root_nd, root_nd, "cpnum", 0, 0, MXML_DESCEND_FIRST);
  setCp(getNodeChar(nd_a), getNodeInteger(nd_b));
  
  nd = mxmlFindElement(root_nd, root_nd, "switch", 0, 0, MXML_DESCEND_FIRST);
  setBranchSWIdx(getNodeInteger(nd));

  if (getPointType() != SolUser)
  {
    nd = mxmlFindElement(root_nd, root_nd, "nparx", 0, 0, MXML_DESCEND_FIRST);
    setNEqns(getNodeInteger(nd));
    if (getNEqns() != 0)
    {
      mxml_node_t* parx_nd = mxmlFindElement(root_nd, root_nd, "parx", 0, 0, MXML_DESCEND_FIRST);
      
      int it = 0;
      for (nd = mxmlFindElement(parx_nd, parx_nd, "par", 0, 0, MXML_DESCEND_FIRST);
      	   nd != 0;
      	   nd = mxmlFindElement(nd, parx_nd->child, "par", 0, 0, MXML_NO_DESCEND) )
      {
      	mxml_node_t* nd_type = mxmlFindElement(nd, nd, "type", 0, 0, MXML_DESCEND_FIRST);
      	mxml_node_t* nd_num = mxmlFindElement(nd, nd, "num", 0, 0, MXML_DESCEND_FIRST);
      	setParX(it, getNodeChar(nd_type), getNodeInteger(nd_num));
      	++it;
      }
    } 
  } else
  {
    nd = mxmlFindElement(root_nd, root_nd, "neqns", 0, 0, MXML_DESCEND_FIRST);
    setNEqns(getNodeInteger(nd));
    if (getNEqns() != 0)
    {
      mxml_node_t* eqns_nd = mxmlFindElement(root_nd, root_nd, "eqns", 0, 0, MXML_DESCEND_FIRST);
      int it = 0;
      for (nd = mxmlFindElement(eqns_nd, eqns_nd, "eqn", 0, 0, MXML_DESCEND_FIRST);
      	   nd != 0;
      	   nd = mxmlFindElement(nd, eqns_nd->child, "eqn", 0, 0, MXML_NO_DESCEND) )
      {
      	mxml_node_t* nd_type = mxmlFindElement(nd, nd, "type", 0, 0, MXML_DESCEND_FIRST);
      	mxml_node_t* nd_num = mxmlFindElement(nd, nd, "num", 0, 0, MXML_DESCEND_FIRST);
      	setEqns(it, getNodeChar(nd_type), getNodeInteger(nd_num));
      	++it;
      }
      mxml_node_t* vars_nd = mxmlFindElement(root_nd, root_nd, "vars", 0, 0, MXML_DESCEND_FIRST);
      it = 0;
      for (nd = mxmlFindElement(vars_nd, vars_nd, "var", 0, 0, MXML_DESCEND_FIRST);
      	   nd != 0;
      	   nd = mxmlFindElement(nd, vars_nd->child, "var", 0, 0, MXML_NO_DESCEND) )
      {
      	mxml_node_t* nd_type = mxmlFindElement(nd, nd, "type", 0, 0, MXML_DESCEND_FIRST);
      	mxml_node_t* nd_num = mxmlFindElement(nd, nd, "num", 0, 0, MXML_DESCEND_FIRST);
      	setVars(it, getNodeChar(nd_type), getNodeInteger(nd_num));
      	++it;
      }
    }
  }
  
  nd = mxmlFindElement(root_nd, root_nd, "nint", 0, 0, MXML_DESCEND_FIRST);
  setNInt(getNodeInteger(nd));

  nd = mxmlFindElement(root_nd, root_nd, "ndeg", 0, 0, MXML_DESCEND_FIRST);
  setNDeg(getNodeInteger(nd));
  
  nd = mxmlFindElement(root_nd, root_nd, "nmul", 0, 0, MXML_DESCEND_FIRST);
  setNMul(getNodeInteger(nd));

  nd = mxmlFindElement(root_nd, root_nd, "stab", 0, 0, MXML_DESCEND_FIRST);
  setStab(getNodeInteger(nd) != 0);

  nd = mxmlFindElement(root_nd, root_nd, "nmat", 0, 0, MXML_DESCEND_FIRST);
  setNMat(getNodeInteger(nd));

  nd = mxmlFindElement(root_nd, root_nd, "nint1", 0, 0, MXML_DESCEND_FIRST);
  setNInt1(getNodeInteger(nd));
  
  nd = mxmlFindElement(root_nd, root_nd, "nint2", 0, 0, MXML_DESCEND_FIRST);
  setNInt2(getNodeInteger(nd));
  
  nd = mxmlFindElement(root_nd, root_nd, "ndeg1", 0, 0, MXML_DESCEND_FIRST);
  setNDeg1(getNodeInteger(nd));
  
  nd = mxmlFindElement(root_nd, root_nd, "ndeg2", 0, 0, MXML_DESCEND_FIRST);
  setNDeg2(getNodeInteger(nd));
  
  nd = mxmlFindElement(root_nd, root_nd, "steps", 0, 0, MXML_DESCEND_FIRST);
  setSteps(getNodeInteger(nd));
  
  nd = mxmlFindElement(root_nd, root_nd, "iad", 0, 0, MXML_DESCEND_FIRST);
  setIad(getNodeInteger(nd));

  nd = mxmlFindElement(root_nd, root_nd, "cpmin", 0, 0, MXML_DESCEND_FIRST);
  setCpMin(getNodeReal(nd));

  nd = mxmlFindElement(root_nd, root_nd, "cpmax", 0, 0, MXML_DESCEND_FIRST);
  setCpMax(getNodeReal(nd));
  
  nd = mxmlFindElement(root_nd, root_nd, "ds", 0, 0, MXML_DESCEND_FIRST);
  setDs(getNodeReal(nd));
  
  nd = mxmlFindElement(root_nd, root_nd, "dsmin", 0, 0, MXML_DESCEND_FIRST);
  setDsMin(getNodeReal(nd));
  
  nd = mxmlFindElement(root_nd, root_nd, "dsmax", 0, 0, MXML_DESCEND_FIRST);
  setDsMax(getNodeReal(nd));
  
  nd = mxmlFindElement(root_nd, root_nd, "dsstart", 0, 0, MXML_DESCEND_FIRST);
  setDsStart(getNodeReal(nd));

  nd = mxmlFindElement(root_nd, root_nd, "epsc", 0, 0, MXML_DESCEND_FIRST);
  setEpsC(getNodeReal(nd));
  
  nd = mxmlFindElement(root_nd, root_nd, "epsr", 0, 0, MXML_DESCEND_FIRST);
  setEpsR(getNodeReal(nd));
  
  nd = mxmlFindElement(root_nd, root_nd, "epsk", 0, 0, MXML_DESCEND_FIRST);
  setEpsK(getNodeReal(nd));

  nd = mxmlFindElement(root_nd, root_nd, "nitc", 0, 0, MXML_DESCEND_FIRST);
  setNItC(getNodeInteger(nd));
  
  nd = mxmlFindElement(root_nd, root_nd, "nitr", 0, 0, MXML_DESCEND_FIRST);
  setNItR(getNodeInteger(nd));
  
  nd = mxmlFindElement(root_nd, root_nd, "nitk", 0, 0, MXML_DESCEND_FIRST);
  setNItK(getNodeInteger(nd));

  nd = mxmlFindElement(root_nd, root_nd, "nsym", 0, 0, MXML_DESCEND_FIRST);
  setNSym(getNodeInteger(nd));

  mxml_node_t* sym_nd = mxmlFindElement(root_nd, root_nd, "sym", 0, 0, MXML_DESCEND_FIRST);
  int it = 0;
  for (nd = mxmlFindElement(sym_nd, sym_nd, "dim", 0, 0, MXML_DESCEND_FIRST);
  	   nd != 0;
  	   nd = mxmlFindElement(nd, sym_nd->child, "dim", 0, 0, MXML_NO_DESCEND))
  {
  	mxml_node_t* nd_real = mxmlFindElement(nd, nd, "real", 0, 0, MXML_DESCEND_FIRST);
  	mxml_node_t* nd_imag = mxmlFindElement(nd, nd, "imag", 0, 0, MXML_DESCEND_FIRST);
  	setSymRe(it, getNodeInteger(nd_real));
    setSymIm(it, getNodeInteger(nd_imag));
    ++it; 
  }
}

void NConstants::saveXmlFile(const std::string &fileName)
{
  std::ofstream file(fileName.c_str());
  printXmlFile(file);
}

void NConstants::printXmlFile(std::ostream& file)
{
  char cbuf[2];
  mxml_node_t *node = 0;
  
  mxml_node_t *xml = mxmlNewXML("cfile 1.0");
  mxml_node_t *data = mxmlNewElement(xml, "knut");
  
  node = mxmlNewElement(data, "input");
  mxmlNewText(node, 0, getInputFile().c_str());
  
  node = mxmlNewElement(data, "output");
  mxmlNewText(node, 0, getOutputFile().c_str());
  
  node = mxmlNewElement(data, "sysname");
  mxmlNewText(node, 0, getSysName().c_str());
  
  node = mxmlNewElement(data, "label");
  mxmlNewInteger(node, getLabel());
  
  node = mxmlNewElement(data, "pointtype");
  mxmlNewInteger(node, getPointType());
  
  node = mxmlNewElement(data, "cptype");
  mxmlNewText(node, 0, c2s(cbuf,getCpType()));

  node = mxmlNewElement(data, "cpnum");
  mxmlNewInteger(node, getCpNum());
  
  node = mxmlNewElement(data, "switch");
  mxmlNewInteger(node, getBranchSW());
 
  if (getPointType() != SolUser)
  {
    node = mxmlNewElement(data, "nparx");
    mxmlNewInteger(node, getNEqns());

    if (getNEqns() != 0)
    {
      mxml_node_t *group_parx = mxmlNewElement(data, "parx");
      for (int i = 0; i < getNEqns(); ++i)
      {
        mxml_node_t *group_par = mxmlNewElement(group_parx, "par");
        
        node = mxmlNewElement(group_par, "type");
        mxmlNewText(node, 0, c2s(cbuf,getParXType(i)));
        
        node = mxmlNewElement(group_par, "num");
        mxmlNewInteger(node, getParXNum(i));
      }
    }
  }
  else
  {
    node = mxmlNewElement(data, "neqns");
    mxmlNewInteger(node, getNEqns());

    if (getNEqns() != 0)
    {
      mxml_node_t *group_eqns = mxmlNewElement(data, "eqns");
      mxml_node_t *group_vars = mxmlNewElement(data, "vars");

      for (int i = 0; i < getNEqns(); ++i)
      {
        mxml_node_t *group_eqn = mxmlNewElement(group_eqns, "eqn");
        
        node = mxmlNewElement(group_eqn, "type");
        mxmlNewText(node, 0, c2s(cbuf,getEqnsType(i)));
        
        node = mxmlNewElement(group_eqn, "num");
        mxmlNewInteger(node, getEqnsNum(i));
      }
      for (int i = 0; i < getNEqns(); ++i)
      {        
        mxml_node_t *group_var = mxmlNewElement(group_vars, "var");

        node = mxmlNewElement(group_var, "type");
        mxmlNewText(node, 0, c2s(cbuf,getVarsType(i)));

        node = mxmlNewElement(group_var, "num");
        mxmlNewInteger(node, getVarsNum(i));
      }
    }
  }
  node = mxmlNewElement(data, "nint");
  mxmlNewInteger(node, getNInt());

  node = mxmlNewElement(data, "ndeg");
  mxmlNewInteger(node, getNDeg());
  
  node = mxmlNewElement(data, "nmul");
  mxmlNewInteger(node, getNMul());
  
  node = mxmlNewElement(data, "stab");
  mxmlNewInteger(node, getStab());

  node = mxmlNewElement(data, "nmat");
  mxmlNewInteger(node, getNMat());
  
  node = mxmlNewElement(data, "nint1");
  mxmlNewInteger(node, getNInt1());

  node = mxmlNewElement(data, "nint2");
  mxmlNewInteger(node, getNInt2());
  
  node = mxmlNewElement(data, "ndeg1");
  mxmlNewInteger(node, getNDeg1());
  
  node = mxmlNewElement(data, "ndeg2");
  mxmlNewInteger(node, getNDeg2());
  
  node = mxmlNewElement(data, "steps");
  mxmlNewInteger(node, getSteps());

  node = mxmlNewElement(data, "iad");
  mxmlNewInteger(node, getIad());
  
  node = mxmlNewElement(data, "cpmin");
  mxmlNewReal(node, getCpMin());
  
  node = mxmlNewElement(data, "cpmax");
  mxmlNewReal(node, getCpMax());
  
  node = mxmlNewElement(data, "ds");
  mxmlNewReal(node, getDs());
  
  node = mxmlNewElement(data, "dsmin");
  mxmlNewReal(node, getDsMin());
  
  node = mxmlNewElement(data, "dsmax");
  mxmlNewReal(node, getDsMax());
  
  node = mxmlNewElement(data, "dsstart");
  mxmlNewReal(node, getDsStart());

  node = mxmlNewElement(data, "epsc");
  mxmlNewReal(node, getEpsC());
  
  node = mxmlNewElement(data, "epsr");
  mxmlNewReal(node, getEpsR());
  
  node = mxmlNewElement(data, "epsk");
  mxmlNewReal(node, getEpsK());

  node = mxmlNewElement(data, "nitc");
  mxmlNewInteger(node, getNItC());
  
  node = mxmlNewElement(data, "nitr");
  mxmlNewInteger(node, getNItR());
  
  node = mxmlNewElement(data, "nitk");
  mxmlNewInteger(node, getNItK());

  node = mxmlNewElement(data, "nsym");
  mxmlNewInteger(node, getNSym());

  if (getNSym() != 0)
  {
    mxml_node_t *group_sym = mxmlNewElement(data, "sym");
    for (int i = 0; i < getNSym(); ++i)
    {
      mxml_node_t *group_dim = mxmlNewElement(group_sym, "dim");
      
      node = mxmlNewElement(group_dim, "real");
      mxmlNewInteger(node, getSymRe(i));
      
      node = mxmlNewElement(group_dim, "imag");
      mxmlNewInteger(node, getSymIm(i));
    }
  }
  
  char *xmlString = mxmlSaveAllocString (xml, knut_whitespace_cb);
  if (xmlString)
  {
   file << xmlString;
   free(xmlString);
  }
  mxmlDelete(xml);
}

int NConstants::toEqnVar(System& sys,
                         Array1D<Eqn>& eqn, Array1D<Var>& var,                 // input
                         Array1D<Eqn>& eqn_refine, Array1D<Var>& var_refine,   // output
                         Array1D<Eqn>& eqn_start, Array1D<Var>& var_start, Eqn& testFN)
{
  // initializing the equations and variables
  if (getPointType() == SolUser)
  {
    eqn.Init(getNEqns());
    var.Init(getNEqns());
    for (int i = 0; i < getNEqns(); i++)
    {
      eqn(i) = getEqns(i);
      var(i) = getVars(i);
    }
  }
  else
  {
    Array1D<Var> L_PARX(getNEqns());
    for (int i = 0; i < getNEqns(); i++)
    {
      L_PARX(i) = getParX(i);
    }
    setBranchSW(PtToEqnVar(eqn, var, getPointType(), L_PARX, sys.npar()));
  }

  // checking whether it is an autonomous problem or not
  bool aut = false;
  bool phaseRot = false;
  testFN = EqnNone;
  int  testFN_idx = -1;
  for (int i = 0; i < eqn.Size(); i++)
  {
    if (eqn(i) == EqnPhase) aut = true;
    if (eqn(i) == EqnPhaseRot) phaseRot = true;
    if ((eqn(i) == EqnTFPD) ||
        (eqn(i) == EqnTFLP) ||
        (eqn(i) == EqnTFLPAUT) ||
        (eqn(i) == EqnTFLPAUTROT) ||
        (eqn(i) == EqnTFCPLX_RE))
    {
      if (testFN_idx == -1)
      {
        testFN = eqn(i);
        testFN_idx = i;
      }
      else P_MESSAGE("Too many test functionals.");
    }
  }

  // setting up for refinement
  if (aut)
  {
    if (phaseRot)
    {
      // std::cout<<"Phase and PhaseRot\n";
      eqn_refine.Init(3);
      var_refine.Init(3);
      eqn_refine(0) = EqnSol;
      eqn_refine(1) = EqnPhase;
      eqn_refine(2) = EqnPhaseRot;
      var_refine(0) = VarSol;
      var_refine(1) = var(var.Size() - 2);
      var_refine(2) = var(var.Size() - 1);
    }
    else
    {
      // std::cout<<"Phase\n";
      eqn_refine.Init(2);
      var_refine.Init(2);
      eqn_refine(0) = EqnSol;
      eqn_refine(1) = EqnPhase;
      var_refine(0) = VarSol;
      var_refine(1) = var(var.Size() - 1);
    }
  }
  else
  {
    if (phaseRot)
    {
      // this happens when a steady state solution of a laser is computed
      eqn_refine.Init(2);
      var_refine.Init(2);
      eqn_refine(0) = EqnSol;
      eqn_refine(1) = EqnPhaseRot;
      var_refine(0) = VarSol;
      var_refine(1) = var(var.Size() - 1);
    }
    else
    {
      eqn_refine.Init(1);
      var_refine.Init(1);
      eqn_refine(0) = EqnSol;
      var_refine(0) = VarSol;
    }
  }

  if (getBranchSW() == TFHBSwitch)
  {
      Array1D<Eqn> ee(eqn_refine);
      Array1D<Var> vv(var_refine);
      eqn_refine.Init(ee.Size()-1);
      var_refine.Init(ee.Size()-1);
      for (int i = 0, j = 0; i < ee.Size(); ++i)
      {
        if (ee(i) != EqnPhase) { eqn_refine(j) = ee(i); var_refine(j) = vv(j); ++j; }
      }
  }

  // Here, we set up the branch switching.
  // We suppose that if there is a switch we use one parameter continuation afterwards
  // without using characteristic matrices. This means that we can switch on the characteristic matrix,
  // include the equation for the eigenvector norm before the other equations and
  // add CP to the variables as a normal parameter.
  Eqn eqn_temp;
  switch (getBranchSW())
  {
/// with TEST FUNCTIONALS
    case TFBRSwitch:
      if (aut) eqn_temp = EqnTFLPAUT;
      else eqn_temp = EqnTFLP;
      goto tfskip;
    case TFPDSwitch:
      eqn_temp = EqnTFPD;
      goto tfskip;
tfskip:
      eqn_start.Init(eqn_refine.Size() + 1);
      var_start.Init(var_refine.Size() + 1);
      eqn_start(0) = eqn_refine(0);
      var_start(0) = var_refine(0);
      eqn_start(1) = eqn_temp;
      var_start(var_refine.Size()) = getCp();
      for (int i = 1; i < eqn_refine.Size(); i++)
      {
        eqn_start(i + 1) = eqn_refine(i);
        var_start(i) = var_refine(i);
      }
      testFN = eqn_temp;
      break;
    case TFTRSwitch:
    case TFHBSwitch:
      eqn_start.Init(eqn_refine.Size() + 2);
      var_start.Init(var_refine.Size() + 2);
      eqn_start(0) = eqn_refine(0);
      var_start(0) = var_refine(0);
      eqn_start(1) = EqnTFCPLX_RE;
      eqn_start(2) = EqnTFCPLX_IM;
      var_start(1) = (Var)(VarPAR0 + sys.npar() + ParAngle); // CH
      var_start(var_refine.Size() + 1) = getCp();
      for (int i = 1; i < eqn_refine.Size(); i++)
      {
        eqn_start(i + 2) = eqn_refine(i);
        var_start(i + 1) = var_refine(i);
      }
      testFN = EqnTFCPLX_RE;
      break;
    default:
      eqn_start.Init(eqn.Size());
      var_start.Init(var.Size());
      eqn_start = eqn;
      var_start = var;
      break;
  }
  if (!aut) return 0;
  else if (!phaseRot) return 1;
  else return 2;
}
