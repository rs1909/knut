// ------------------------------------------------------------------------- //
//
// This is part of KNUT
// Copyright (c) 2012 by Robert Szalai
//
// For license, see the file COPYING in the root directory of the package
//
// ------------------------------------------------------------------------- //

#include "config.h"
#include "constants.h"
#include <iostream>
#include <fstream>
#include <sstream>

// For loading XML constant files
#include "mxml.h"

// from constants.cpp
const char *knut_whitespace_cb(mxml_node_t *node, int where);

static inline void mxmlNewDouble(mxml_node_t *node, double dbl)
{
  mxmlNewTextf(node, 0, "%.15E", dbl);
}

static inline mxml_node_t *setNodeIndex(mxml_node_t *parent, size_t idx)
{
  std::ostringstream strstr;
  strstr << idx;
  return mxmlNewOpaque (parent, strstr.str().c_str());
}

static inline const char *getNodeText(mxml_node_t* nd)
{
  if (nd != 0)
  {
    if (nd->child != 0)
    {
      if (nd->child == nd->last_child)
      {
        if (nd->child->type == MXML_OPAQUE)
        {
          return nd->child->value.opaque; //text.string;
        } else P_MESSAGE1("MXML: wrong element type.");
      } else P_MESSAGE1("MXML: node has too many children.");
    } else return 0;
  } else return 0;
}

static inline const char *getNodeText(mxml_node_t* nd, const char* def_value)
{
  const char* str = getNodeText(nd);
//#ifdef DEBUG
//  P_ERROR_X1(str != 0, "MXML: node could not be found at line.");
//#endif
  if (str) return str;
  else return def_value;
}

static inline char getNodeChar(mxml_node_t* nd, const char *line = "")
{
  const char* str = getNodeText(nd);
  P_ERROR_X3(str != 0, "MXML: node could not be found at line ", line, ".");
  return str[0];
}

static inline int getNodeInteger(mxml_node_t* nd, const char *line = "")
{
  const char* str = getNodeText(nd);
  P_ERROR_X3(str != 0, "MXML: node could not be found at line ", line, ".");
  return toInt(strtol(str, 0, 10));
}

static inline int getNodeInteger(mxml_node_t* nd, int def_value)
{
  const char* str = getNodeText(nd);
//#ifdef DEBUG
//  P_ERROR_X1(str != 0, "MXML: node could not be found at line.");
//#endif
  if (str != 0) return toInt(strtol(str, 0, 10));
  else return def_value;
}

static inline size_t getNodeIndex(mxml_node_t* nd, const char *line = "")
{
  const char* str = getNodeText(nd);
  P_ERROR_X3(str != 0, "MXML: node could not be found at line ", line, ".");
  return toSizeT(strtol(str, 0, 10));
}

static inline size_t getNodeIndex(mxml_node_t* nd, size_t def_value)
{
  const char* str = getNodeText(nd);
//#ifdef DEBUG
//  P_ERROR_X1(str != 0, "MXML: node could not be found at line.");
//#endif
  if (str != 0) return toSizeT(strtol(str, 0, 10));
  else return def_value;
}

static inline double getNodeReal(mxml_node_t* nd, const char *line = "")
{
  const char* str = getNodeText(nd);
  P_ERROR_X3(str != 0, "MXML: node could not be found at line ", line, ".");
  std::istringstream strstr(str);
  double val;
  strstr >> val;
  return val;
}

static inline double getNodeReal(mxml_node_t* nd, double def_value)
{
  const char* str = getNodeText(nd);
//#ifdef DEBUG
//  P_ERROR_X1(str != 0, "MXML: node could not be found at line.");
//#endif
  if (str != 0)
  {
    std::istringstream strstr(str);
    double val;
    strstr >> val;
    return val;
  } else return def_value;
}

#define STRX(x) #x
#define STR(x) STRX(x)
#define getNodeCharM(nd) getNodeChar(nd,STR(__LINE__))
#define getNodeIntegerM(nd) getNodeInteger(nd,STR(__LINE__))
#define getNodeIndexM(nd) getNodeIndex(nd,STR(__LINE__))
#define getNodeRealM(nd) getNodeReal(nd,STR(__LINE__))

void KNConstantsBase::loadXmlFileV4(const std::string &fileName)
{
  FILE *fp;
  mxml_node_t *tree;

  fp = fopen(fileName.c_str(), "r");
  P_ERROR_X3(fp != 0, "Cannot open file '", fileName, "'." );
  tree = mxmlLoadFile(0, fp, MXML_OPAQUE_CALLBACK);
  fclose(fp);

  mxml_node_t* nd = 0;
  mxml_node_t* root_nd = mxmlFindElement(tree, tree, "knut", 0, 0, MXML_DESCEND_FIRST);
  if (!root_nd)
  {
    root_nd = mxmlFindElement(tree, tree, "pdde", 0, 0, MXML_DESCEND_FIRST);
  }
  
  nd = mxmlFindElement(root_nd, root_nd, "input", 0, 0, MXML_DESCEND_FIRST);
  setInputFile(getNodeText(nd, ""));
  
  nd = mxmlFindElement(root_nd, root_nd, "output", 0, 0, MXML_DESCEND_FIRST);
  setOutputFile(getNodeText(nd, ""));
  
  nd = mxmlFindElement(root_nd, root_nd, "sysname", 0, 0, MXML_DESCEND_FIRST);
  // This has to load the system definition
  const char *s_type = mxmlElementGetAttr(nd, "type");
  if (s_type) setSysType(s_type);
  setSysNameText(getNodeText(nd, ""));

  nd = mxmlFindElement(root_nd, root_nd, "fromtype", 0, 0, MXML_DESCEND_FIRST);
  setFromType(BifTypeTable.CodeToType(getNodeText(nd))); // -> default

  nd = mxmlFindElement(root_nd, root_nd, "label", 0, 0, MXML_DESCEND_FIRST);
  setLabel(getNodeIndexM(nd));
  
  nd = mxmlFindElement(root_nd, root_nd, "pointtype", 0, 0, MXML_DESCEND_FIRST);
  setPointType(PtTypeTable.CodeToType(getNodeText(nd)));
    
  nd = mxmlFindElement(root_nd, root_nd, "cp", 0, 0, MXML_DESCEND_FIRST);
  setCp(VarTable.CodeToType(getNodeText(nd)));
  
  nd = mxmlFindElement(root_nd, root_nd, "switch", 0, 0, MXML_DESCEND_FIRST);
  setBranchSW(BranchSWTable.CodeToType(getNodeText(nd))); // -> default

  if (getPointType() != SolUser)
  {
    nd = mxmlFindElement(root_nd, root_nd, "nparx", 0, 0, MXML_DESCEND_FIRST);
    setParxSize(getNodeIndexM(nd));
    if (getParxSize() != 0)
    {
      mxml_node_t* parx_nd = mxmlFindElement(root_nd, root_nd, "parx", 0, 0, MXML_DESCEND_FIRST);
      
      size_t it = 0;
      for (nd = mxmlFindElement(parx_nd, parx_nd, "par", 0, 0, MXML_DESCEND_FIRST);
           nd != 0;
           nd = mxmlFindElement(nd, parx_nd->child, "par", 0, 0, MXML_NO_DESCEND) )
      {
        setParx(it, VarTable.CodeToType(getNodeText(nd)));
        ++it;
      }
    }
  } else
  {
    nd = mxmlFindElement(root_nd, root_nd, "neqns", 0, 0, MXML_DESCEND_FIRST);
    setEqnsSize(getNodeIndexM(nd));
    setVarsSize(getNodeIndexM(nd));
    if (getEqnsSize() != 0)
    {
      mxml_node_t* eqns_nd = mxmlFindElement(root_nd, root_nd, "eqns", 0, 0, MXML_DESCEND_FIRST);
      size_t it = 0;
      for (nd = mxmlFindElement(eqns_nd, eqns_nd, "eqn", 0, 0, MXML_DESCEND_FIRST);
           nd != 0;
           nd = mxmlFindElement(nd, eqns_nd->child, "eqn", 0, 0, MXML_NO_DESCEND) )
      {
        setEqns(it, EqnTable.CodeToType(getNodeText(nd)));
        ++it;
      }
      mxml_node_t* vars_nd = mxmlFindElement(root_nd, root_nd, "vars", 0, 0, MXML_DESCEND_FIRST);
      it = 0;
      for (nd = mxmlFindElement(vars_nd, vars_nd, "var", 0, 0, MXML_DESCEND_FIRST);
           nd != 0;
           nd = mxmlFindElement(nd, vars_nd->child, "var", 0, 0, MXML_NO_DESCEND) )
      {
        setVars(it, VarTable.CodeToType(getNodeText(nd)));
        ++it;
      }
    }
  }
  nd = mxmlFindElement(root_nd, root_nd, "nint", 0, 0, MXML_DESCEND_FIRST);
  setNInt(getNodeIndex(nd, 20));

  nd = mxmlFindElement(root_nd, root_nd, "ndeg", 0, 0, MXML_DESCEND_FIRST);
  setNDeg(getNodeIndex(nd, 5));
  
  nd = mxmlFindElement(root_nd, root_nd, "nmul", 0, 0, MXML_DESCEND_FIRST);
  setNMul(getNodeIndex(nd, 5));

  nd = mxmlFindElement(root_nd, root_nd, "stab", 0, 0, MXML_DESCEND_FIRST);
  setStab(getNodeIndex(nd, (size_t)0) != 0);
  
  nd = mxmlFindElement(root_nd, root_nd, "curvature", 0, 0, MXML_DESCEND_FIRST);
  setCAngle(getNodeReal(nd, getCAngle()));
  
  nd = mxmlFindElement(root_nd, root_nd, "nint1", 0, 0, MXML_DESCEND_FIRST);
  setNInt1(getNodeIndex(nd, 12));
  
  nd = mxmlFindElement(root_nd, root_nd, "nint2", 0, 0, MXML_DESCEND_FIRST);
  setNInt2(getNodeIndex(nd, 12));
  
  nd = mxmlFindElement(root_nd, root_nd, "ndeg1", 0, 0, MXML_DESCEND_FIRST);
  setNDeg1(getNodeIndex(nd, 4));
  
  nd = mxmlFindElement(root_nd, root_nd, "ndeg2", 0, 0, MXML_DESCEND_FIRST);
  setNDeg2(getNodeIndex(nd, 4));
  
  nd = mxmlFindElement(root_nd, root_nd, "steps", 0, 0, MXML_DESCEND_FIRST);
  setSteps(getNodeIndex(nd, 100));
  
  nd = mxmlFindElement(root_nd, root_nd, "iad", 0, 0, MXML_DESCEND_FIRST);
  setIad(getNodeIndex(nd, 3));
  
  nd = mxmlFindElement(root_nd, root_nd, "npr", 0, 0, MXML_DESCEND_FIRST);
  setNPr(getNodeIndex(nd, 10));
  P_ERROR_X1(getNPr() > 0, "NPR must be greater than zero." );
  
  nd = mxmlFindElement(root_nd, root_nd, "cpmin", 0, 0, MXML_DESCEND_FIRST);
  setCpMin(getNodeRealM(nd));

  nd = mxmlFindElement(root_nd, root_nd, "cpmax", 0, 0, MXML_DESCEND_FIRST);
  setCpMax(getNodeReal(nd));
  
  nd = mxmlFindElement(root_nd, root_nd, "ds", 0, 0, MXML_DESCEND_FIRST);
  setDs(getNodeRealM(nd));
  
  nd = mxmlFindElement(root_nd, root_nd, "dsmin", 0, 0, MXML_DESCEND_FIRST);
  setDsMin(getNodeRealM(nd));
  
  nd = mxmlFindElement(root_nd, root_nd, "dsmax", 0, 0, MXML_DESCEND_FIRST);
  setDsMax(getNodeRealM(nd));
  
  nd = mxmlFindElement(root_nd, root_nd, "dsstart", 0, 0, MXML_DESCEND_FIRST);
  setDsStart(getNodeRealM(nd));

  nd = mxmlFindElement(root_nd, root_nd, "epsc", 0, 0, MXML_DESCEND_FIRST);
  setEpsC(getNodeRealM(nd));
  
  nd = mxmlFindElement(root_nd, root_nd, "epsr", 0, 0, MXML_DESCEND_FIRST);
  setEpsR(getNodeRealM(nd));
  
  nd = mxmlFindElement(root_nd, root_nd, "epsk", 0, 0, MXML_DESCEND_FIRST);
  setEpsK(getNodeRealM(nd));

  nd = mxmlFindElement(root_nd, root_nd, "nitc", 0, 0, MXML_DESCEND_FIRST);
  setNItC(getNodeIndex(nd, 5));
  
  nd = mxmlFindElement(root_nd, root_nd, "nitr", 0, 0, MXML_DESCEND_FIRST);
  setNItR(getNodeIndex(nd, 12));
  
  nd = mxmlFindElement(root_nd, root_nd, "nitk", 0, 0, MXML_DESCEND_FIRST);
  setNItK(getNodeIndex(nd, 12));

  nd = mxmlFindElement(root_nd, root_nd, "nderi", 0, 0, MXML_DESCEND_FIRST);
  setNDeri(getNodeIndex(nd, 8)); // this is just a number, higher then sys_nderi() provides

  nd = mxmlFindElement(root_nd, root_nd, "nsym", 0, 0, MXML_DESCEND_FIRST);
  setSymReSize(getNodeIndex(nd, (size_t)0));
  setSymImSize(getNodeIndex(nd, (size_t)0));

  mxml_node_t* sym_nd = mxmlFindElement(root_nd, root_nd, "sym", 0, 0, MXML_DESCEND_FIRST);
  size_t it = 0;
  for (nd = mxmlFindElement(sym_nd, sym_nd, "dim", 0, 0, MXML_DESCEND_FIRST);
       nd != 0;
       nd = mxmlFindElement(nd, sym_nd->child, "dim", 0, 0, MXML_NO_DESCEND))
  {
    mxml_node_t* nd_real = mxmlFindElement(nd, nd, "real", 0, 0, MXML_DESCEND_FIRST);
    mxml_node_t* nd_imag = mxmlFindElement(nd, nd, "imag", 0, 0, MXML_DESCEND_FIRST);
    setSymRe(it, getNodeIndexM(nd_real));
    setSymIm(it, getNodeIndexM(nd_imag));
    ++it; 
  }
  mxmlDelete(tree);
}

void KNConstantsBase::printXmlFileV4(std::ostream& file)
{
  mxml_node_t *node = 0;
  
  mxml_node_t *xml = mxmlNewXML("cfile 1.0");
  mxml_node_t *data = mxmlNewElement(xml, "knut");
  mxmlElementSetAttr( data, "version", PACKAGE_VERSION);
  
  node = mxmlNewElement(data, "input");
  mxmlNewText(node, 0, getInputFile().c_str());
  
  node = mxmlNewElement(data, "output");
  mxmlNewText(node, 0, getOutputFile().c_str());
  
  node = mxmlNewElement(data, "sysname");
  mxmlNewText(node, 0, getSysName().c_str());
  mxmlElementSetAttr(node, "type", getSysType().c_str());
  
  node = mxmlNewElement(data, "fromtype");
  mxmlNewText(node, 0, BifTypeTable.TypeToCode(getFromType()).c_str());
  
  node = mxmlNewElement(data, "label");
  setNodeIndex(node, getLabel());
  
  node = mxmlNewElement(data, "pointtype");
  mxmlNewText(node, 0, PtTypeTable.TypeToCode(getPointType()).c_str());
  
  node = mxmlNewElement(data, "cp");
  mxmlNewText(node, 0, VarTable.TypeToCode(getCp()).c_str());

  node = mxmlNewElement(data, "switch");
  mxmlNewText(node, 0, BranchSWTable.TypeToCode(getBranchSW()).c_str());
 
  if (getPointType() != SolUser)
  {
    node = mxmlNewElement(data, "nparx");
    setNodeIndex(node, getParxSize());

    if (getParxSize() != 0)
    {
      mxml_node_t *group_parx = mxmlNewElement(data, "parx");
      for (size_t i = 0; i < getParxSize(); ++i)
      {
        node = mxmlNewElement(group_parx, "par");
        mxmlNewText(node, 0, VarTable.TypeToCode(getParx(i)).c_str());
      }
    }
  }
  else
  {
    node = mxmlNewElement(data, "neqns");
    setNodeIndex(node, getEqnsSize());

    if (getEqnsSize() != 0)
    {
      mxml_node_t *group_eqns = mxmlNewElement(data, "eqns");
      mxml_node_t *group_vars = mxmlNewElement(data, "vars");

      for (size_t i = 0; i < getEqnsSize(); ++i)
      {
        node = mxmlNewElement(group_eqns, "eqn");
        mxmlNewText(node, 0, EqnTable.TypeToCode(getEqns(i)).c_str());
      }
      for (size_t i = 0; i < getVarsSize(); ++i)
      {        
        node = mxmlNewElement(group_vars, "var");
        mxmlNewText(node, 0, VarTable.TypeToCode(getVars(i)).c_str());
      }
    }
  }
  node = mxmlNewElement(data, "nint");
  setNodeIndex(node, getNInt());

  node = mxmlNewElement(data, "ndeg");
  setNodeIndex(node, getNDeg());
  
  node = mxmlNewElement(data, "nmul");
  setNodeIndex(node, getNMul());
  
  node = mxmlNewElement(data, "stab");
  setNodeIndex(node, getStab());
  
  node = mxmlNewElement(data, "curvature");
  mxmlNewDouble(node, getCAngle());
  
  node = mxmlNewElement(data, "nint1");
  setNodeIndex(node, getNInt1());

  node = mxmlNewElement(data, "nint2");
  setNodeIndex(node, getNInt2());
  
  node = mxmlNewElement(data, "ndeg1");
  setNodeIndex(node, getNDeg1());
  
  node = mxmlNewElement(data, "ndeg2");
  setNodeIndex(node, getNDeg2());
  
  node = mxmlNewElement(data, "steps");
  setNodeIndex(node, getSteps());

  node = mxmlNewElement(data, "iad");
  setNodeIndex(node, getIad());
  
  node = mxmlNewElement(data, "npr");
  setNodeIndex(node, getNPr());
  
  node = mxmlNewElement(data, "cpmin");
  mxmlNewDouble(node, getCpMin());
  
  node = mxmlNewElement(data, "cpmax");
  mxmlNewDouble(node, getCpMax());
  
  node = mxmlNewElement(data, "ds");
  mxmlNewDouble(node, getDs());
  
  node = mxmlNewElement(data, "dsmin");
  mxmlNewDouble(node, getDsMin());
  
  node = mxmlNewElement(data, "dsmax");
  mxmlNewDouble(node, getDsMax());
  
  node = mxmlNewElement(data, "dsstart");
  mxmlNewDouble(node, getDsStart());

  node = mxmlNewElement(data, "epsc");
  mxmlNewDouble(node, getEpsC());
  
  node = mxmlNewElement(data, "epsr");
  mxmlNewDouble(node, getEpsR());
  
  node = mxmlNewElement(data, "epsk");
  mxmlNewDouble(node, getEpsK());

  node = mxmlNewElement(data, "nitc");
  setNodeIndex(node, getNItC());
  
  node = mxmlNewElement(data, "nitr");
  setNodeIndex(node, getNItR());
  
  node = mxmlNewElement(data, "nitk");
  setNodeIndex(node, getNItK());
  
  node = mxmlNewElement(data, "nderi");
  setNodeIndex(node, getNDeri());

  node = mxmlNewElement(data, "nsym");
  setNodeIndex(node, getSymReSize());

  if (getSymReSize() != 0)
  {
    mxml_node_t *group_sym = mxmlNewElement(data, "sym");
    for (size_t i = 0; i < getSymReSize(); ++i)
    {
      mxml_node_t *group_dim = mxmlNewElement(group_sym, "dim");
      
      node = mxmlNewElement(group_dim, "real");
      setNodeIndex(node, getSymRe(i));
      
      node = mxmlNewElement(group_dim, "imag");
      setNodeIndex(node, getSymIm(i));
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

void KNConstantsBase::saveXmlFileV4(const std::string &fileName)
{
  std::ofstream file(fileName.c_str());
  printXmlFileV4(file);
}

void KNConstantsBase::loadXmlFileV2(const std::string &fileName)
{
  FILE *fp;
  mxml_node_t *tree;

  fp = fopen(fileName.c_str(), "r");
  P_ERROR_X3(fp != 0, "Cannot open file '", fileName, "'." );
  tree = mxmlLoadFile(0, fp, MXML_OPAQUE_CALLBACK);
  fclose(fp);

  mxml_node_t* nd = 0, *nd_a = 0, *nd_b = 0;
  mxml_node_t* root_nd = mxmlFindElement(tree, tree, "knut", 0, 0, MXML_DESCEND_FIRST);
  if (!root_nd)
  {
    root_nd = mxmlFindElement(tree, tree, "pdde", 0, 0, MXML_DESCEND_FIRST);
  }
  
  nd = mxmlFindElement(root_nd, root_nd, "input", 0, 0, MXML_DESCEND_FIRST);
  setInputFile(getNodeText(nd, ""));
  
  nd = mxmlFindElement(root_nd, root_nd, "output", 0, 0, MXML_DESCEND_FIRST);
  setOutputFile(getNodeText(nd, ""));
  
  nd = mxmlFindElement(root_nd, root_nd, "sysname", 0, 0, MXML_DESCEND_FIRST);
  // This has to load the system definition
  setSysNameText(getNodeText(nd, ""));

  nd = mxmlFindElement(root_nd, root_nd, "fromtype", 0, 0, MXML_DESCEND_FIRST);
  setFromType(static_cast<BifType>(getNodeIndex(nd, BifNone)));

  nd = mxmlFindElement(root_nd, root_nd, "label", 0, 0, MXML_DESCEND_FIRST);
  setLabel(getNodeIndexM(nd));
  
  nd = mxmlFindElement(root_nd, root_nd, "pointtype", 0, 0, MXML_DESCEND_FIRST);
  setPointType(static_cast<PtType>(getNodeIndexM(nd)));
    
  nd_a = mxmlFindElement(root_nd, root_nd, "cptype", 0, 0, MXML_DESCEND_FIRST);
  const char cp_type = getNodeCharM(nd_a); 
  
  nd_b = mxmlFindElement(root_nd, root_nd, "cpnum", 0, 0, MXML_DESCEND_FIRST);
  const size_t cp_num = getNodeIndexM(nd_b);
  
  setCp(varFromTypeNum(cp_type,cp_num));
  
  nd = mxmlFindElement(root_nd, root_nd, "switch", 0, 0, MXML_DESCEND_FIRST);
  setBranchSW(static_cast<BranchSW>(getNodeIndex(nd, NOSwitch)));

  if (getPointType() != SolUser)
  {
    nd = mxmlFindElement(root_nd, root_nd, "nparx", 0, 0, MXML_DESCEND_FIRST);
    setParxSize(getNodeIndexM(nd));
    if (getParxSize() != 0)
    {
      mxml_node_t* parx_nd = mxmlFindElement(root_nd, root_nd, "parx", 0, 0, MXML_DESCEND_FIRST);
      
      size_t it = 0;
      for (nd = mxmlFindElement(parx_nd, parx_nd, "par", 0, 0, MXML_DESCEND_FIRST);
           nd != 0;
           nd = mxmlFindElement(nd, parx_nd->child, "par", 0, 0, MXML_NO_DESCEND) )
      {
        mxml_node_t* nd_type = mxmlFindElement(nd, nd, "type", 0, 0, MXML_DESCEND_FIRST);
        const char type = getNodeCharM(nd_type);
        mxml_node_t* nd_num = mxmlFindElement(nd, nd, "num", 0, 0, MXML_DESCEND_FIRST);
        const size_t num = getNodeIndexM(nd_num);
        setParx(it, varFromTypeNum(type,num));
        ++it;
      }
    } 
  } else
  {
    nd = mxmlFindElement(root_nd, root_nd, "neqns", 0, 0, MXML_DESCEND_FIRST);
    setEqnsSize(getNodeIndexM(nd));
    setVarsSize(getNodeIndexM(nd));
    if (getEqnsSize() != 0)
    {
      mxml_node_t* eqns_nd = mxmlFindElement(root_nd, root_nd, "eqns", 0, 0, MXML_DESCEND_FIRST);
      size_t it = 0;
      for (nd = mxmlFindElement(eqns_nd, eqns_nd, "eqn", 0, 0, MXML_DESCEND_FIRST);
           nd != 0;
           nd = mxmlFindElement(nd, eqns_nd->child, "eqn", 0, 0, MXML_NO_DESCEND) )
      {
        mxml_node_t* nd_type = mxmlFindElement(nd, nd, "type", 0, 0, MXML_DESCEND_FIRST);
        const char type = getNodeCharM(nd_type);
        mxml_node_t* nd_num = mxmlFindElement(nd, nd, "num", 0, 0, MXML_DESCEND_FIRST);
        const size_t num = getNodeIndexM(nd_num);
        setEqns(it, eqnFromTypeNum(type,num));
        ++it;
      }
      mxml_node_t* vars_nd = mxmlFindElement(root_nd, root_nd, "vars", 0, 0, MXML_DESCEND_FIRST);
      it = 0;
      for (nd = mxmlFindElement(vars_nd, vars_nd, "var", 0, 0, MXML_DESCEND_FIRST);
           nd != 0;
           nd = mxmlFindElement(nd, vars_nd->child, "var", 0, 0, MXML_NO_DESCEND) )
      {
        mxml_node_t* nd_type = mxmlFindElement(nd, nd, "type", 0, 0, MXML_DESCEND_FIRST);
        const char type = getNodeCharM(nd_type);
        mxml_node_t* nd_num = mxmlFindElement(nd, nd, "num", 0, 0, MXML_DESCEND_FIRST);
        const size_t num = getNodeIndexM(nd_num);
        setVars(it, varFromTypeNum(type,num));
        ++it;
      }
    }
  }
  
  nd = mxmlFindElement(root_nd, root_nd, "nint", 0, 0, MXML_DESCEND_FIRST);
  setNInt(getNodeIndex(nd, 20));

  nd = mxmlFindElement(root_nd, root_nd, "ndeg", 0, 0, MXML_DESCEND_FIRST);
  setNDeg(getNodeIndex(nd, 5));
  
  nd = mxmlFindElement(root_nd, root_nd, "nmul", 0, 0, MXML_DESCEND_FIRST);
  setNMul(getNodeIndex(nd, 5));

  nd = mxmlFindElement(root_nd, root_nd, "stab", 0, 0, MXML_DESCEND_FIRST);
  setStab(getNodeIndex(nd, (size_t)0) != 0);
  
  nd = mxmlFindElement(root_nd, root_nd, "curvature", 0, 0, MXML_DESCEND_FIRST);
  setCAngle(getNodeReal(nd, getCAngle()));
  
  nd = mxmlFindElement(root_nd, root_nd, "nint1", 0, 0, MXML_DESCEND_FIRST);
  setNInt1(getNodeIndex(nd, 12));
  
  nd = mxmlFindElement(root_nd, root_nd, "nint2", 0, 0, MXML_DESCEND_FIRST);
  setNInt2(getNodeIndex(nd, 12));
  
  nd = mxmlFindElement(root_nd, root_nd, "ndeg1", 0, 0, MXML_DESCEND_FIRST);
  setNDeg1(getNodeIndex(nd, 4));
  
  nd = mxmlFindElement(root_nd, root_nd, "ndeg2", 0, 0, MXML_DESCEND_FIRST);
  setNDeg2(getNodeIndex(nd, 4));
  
  nd = mxmlFindElement(root_nd, root_nd, "steps", 0, 0, MXML_DESCEND_FIRST);
  setSteps(getNodeIndex(nd, 100));
  
  nd = mxmlFindElement(root_nd, root_nd, "iad", 0, 0, MXML_DESCEND_FIRST);
  setIad(getNodeIndex(nd, 3));
  
  nd = mxmlFindElement(root_nd, root_nd, "npr", 0, 0, MXML_DESCEND_FIRST);
  setNPr(getNodeIndex(nd, 10));
  P_ERROR_X1(getNPr() > 0, "NPR must be greater than zero." );
  
  nd = mxmlFindElement(root_nd, root_nd, "cpmin", 0, 0, MXML_DESCEND_FIRST);
  setCpMin(getNodeRealM(nd));

  nd = mxmlFindElement(root_nd, root_nd, "cpmax", 0, 0, MXML_DESCEND_FIRST);
  setCpMax(getNodeRealM(nd));
  
  nd = mxmlFindElement(root_nd, root_nd, "ds", 0, 0, MXML_DESCEND_FIRST);
  setDs(getNodeRealM(nd));
  
  nd = mxmlFindElement(root_nd, root_nd, "dsmin", 0, 0, MXML_DESCEND_FIRST);
  setDsMin(getNodeRealM(nd));
  
  nd = mxmlFindElement(root_nd, root_nd, "dsmax", 0, 0, MXML_DESCEND_FIRST);
  setDsMax(getNodeRealM(nd));
  
  nd = mxmlFindElement(root_nd, root_nd, "dsstart", 0, 0, MXML_DESCEND_FIRST);
  setDsStart(getNodeRealM(nd));

  nd = mxmlFindElement(root_nd, root_nd, "epsc", 0, 0, MXML_DESCEND_FIRST);
  setEpsC(getNodeRealM(nd));
  
  nd = mxmlFindElement(root_nd, root_nd, "epsr", 0, 0, MXML_DESCEND_FIRST);
  setEpsR(getNodeRealM(nd));
  
  nd = mxmlFindElement(root_nd, root_nd, "epsk", 0, 0, MXML_DESCEND_FIRST);
  setEpsK(getNodeRealM(nd));

  nd = mxmlFindElement(root_nd, root_nd, "nitc", 0, 0, MXML_DESCEND_FIRST);
  setNItC(getNodeIndex(nd, 5));
  
  nd = mxmlFindElement(root_nd, root_nd, "nitr", 0, 0, MXML_DESCEND_FIRST);
  setNItR(getNodeIndex(nd, 12));
  
  nd = mxmlFindElement(root_nd, root_nd, "nitk", 0, 0, MXML_DESCEND_FIRST);
  setNItK(getNodeIndex(nd, 12));

  nd = mxmlFindElement(root_nd, root_nd, "nderi", 0, 0, MXML_DESCEND_FIRST);
  setNDeri(getNodeIndex(nd, 8)); // this is just a number, higher then sys_nderi() provides

  nd = mxmlFindElement(root_nd, root_nd, "nsym", 0, 0, MXML_DESCEND_FIRST);
  setSymReSize(getNodeIndex(nd, (size_t)0));
  setSymImSize(getNodeIndex(nd, (size_t)0));

  mxml_node_t* sym_nd = mxmlFindElement(root_nd, root_nd, "sym", 0, 0, MXML_DESCEND_FIRST);
  size_t it = 0;
  for (nd = mxmlFindElement(sym_nd, sym_nd, "dim", 0, 0, MXML_DESCEND_FIRST);
       nd != 0;
       nd = mxmlFindElement(nd, sym_nd->child, "dim", 0, 0, MXML_NO_DESCEND))
  {
    mxml_node_t* nd_real = mxmlFindElement(nd, nd, "real", 0, 0, MXML_DESCEND_FIRST);
    mxml_node_t* nd_imag = mxmlFindElement(nd, nd, "imag", 0, 0, MXML_DESCEND_FIRST);
    setSymRe(it, getNodeIndexM(nd_real));
    setSymIm(it, getNodeIndexM(nd_imag));
    ++it; 
  }
  mxmlDelete(tree);
}
