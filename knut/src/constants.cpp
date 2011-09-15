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
#include "mxml.h"

std::map<std::string, KNConstantNames::tuple_t> KNConstantNames::nlist;

void KNConstantNames::addConstantName(const char * type, const char * name, void * var, void(*nsetFun)())
{ 
  if (nlist.find (name) == nlist.end())
  {
    tuple_t t = { type, var, nsetFun };
    nlist[name] = t;
//     std::cout << "CNAME " << nlist[name].type << " " << name << "\n";
  }
}

const std::string& KNConstantNames::findType(const char * name)
{
  return nlist[name].type;
}

const void* KNConstantNames::findValue(const char * name)
{
  return nlist[name].var;
}

KNConstantNames::FPTR KNConstantNames::findFun(const char * name)
{
  std::map<std::string, tuple_t>::iterator it = nlist.find (name);
  if (it != nlist.end())
  {
    return nlist[name].nsetFun;
  } else return 0;
}

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

static inline bool inputAssert(std::istream& is)
{
  if (is.rdstate() & (std::istream::eofbit | std::istream::failbit | std::istream::badbit))
  {
    switch (is.rdstate())
    {
      case std::istream::eofbit:
        P_MESSAGE1("std::istream: unexpected end of file.");
        return true;
        break;
      case std::istream::failbit:
        P_MESSAGE1("std::istream: Input failed.");
        return true;
        break;
      case std::istream::badbit:
        P_MESSAGE1("std::istream: bad input.");
        return true;
        break;
      default:
        P_MESSAGE1("std::istream: unexpected error.");
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
  if (str) return str;
  else return def_value;
}
static inline char getNodeChar(mxml_node_t* nd)
{
  const char* str = getNodeText(nd);
  P_ERROR_X1(str != 0, "MXML: node could not be found.");
  return str[0];
}

static inline int getNodeInteger(mxml_node_t* nd)
{
  const char* str = getNodeText(nd);
  P_ERROR_X1(str != 0, "MXML: node could not be found.");
  return toInt(strtol(str, 0, 10));
}

static inline int getNodeInteger(mxml_node_t* nd, int def_value)
{
  const char* str = getNodeText(nd);
  if (str != 0) return toInt(strtol(str, 0, 10));
  else return def_value;
}

static inline double getNodeReal(mxml_node_t* nd)
{
  const char* str = getNodeText(nd);
  P_ERROR_X1(str != 0, "MXML: node could not be found.");
  std::istringstream strstr(str);
  double val;
  strstr >> val;
  return val;
}

static inline double getNodeReal(mxml_node_t* nd, double def_value)
{
  const char* str = getNodeText(nd);
  if (str != 0)
  {
    std::istringstream strstr(str);
    double val;
    strstr >> val;
    return val;
  } else return def_value;
}

static inline const char* c2s(char *buf, char c)
{
  buf[0] = c;
  buf[1] = '\0';
  return buf;
}

void KNConstants::loadXmlFile(const std::string &fileName)
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
  setFromType(BifTypeTable.CodeToType(getNodeText(nd))); // -> default

  nd = mxmlFindElement(root_nd, root_nd, "label", 0, 0, MXML_DESCEND_FIRST);
  setLabel(getNodeInteger(nd));
  
  nd = mxmlFindElement(root_nd, root_nd, "pointtype", 0, 0, MXML_DESCEND_FIRST);
  setPointType(PtTypeTable.CodeToType(getNodeText(nd)));
    
  nd = mxmlFindElement(root_nd, root_nd, "cp", 0, 0, MXML_DESCEND_FIRST);
  setCp(VarTable.CodeToType(getNodeText(nd)));
  
  nd = mxmlFindElement(root_nd, root_nd, "switch", 0, 0, MXML_DESCEND_FIRST);
  setBranchSW(BranchSWTable.CodeToType(getNodeText(nd))); // -> default

  if (getPointType() != SolUser)
  {
    nd = mxmlFindElement(root_nd, root_nd, "nparx", 0, 0, MXML_DESCEND_FIRST);
    setParxSize(getNodeInteger(nd));
    if (getParxSize() != 0)
    {
      mxml_node_t* parx_nd = mxmlFindElement(root_nd, root_nd, "parx", 0, 0, MXML_DESCEND_FIRST);
      
      int it = 0;
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
    setEqnsSize(getNodeInteger(nd));
    setVarsSize(getNodeInteger(nd));
    if (getEqnsSize() != 0)
    {
      mxml_node_t* eqns_nd = mxmlFindElement(root_nd, root_nd, "eqns", 0, 0, MXML_DESCEND_FIRST);
      int it = 0;
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
  setNInt(getNodeInteger(nd, 20));

  nd = mxmlFindElement(root_nd, root_nd, "ndeg", 0, 0, MXML_DESCEND_FIRST);
  setNDeg(getNodeInteger(nd, 5));
  
  nd = mxmlFindElement(root_nd, root_nd, "nmul", 0, 0, MXML_DESCEND_FIRST);
  setNMul(getNodeInteger(nd, 5));

  nd = mxmlFindElement(root_nd, root_nd, "stab", 0, 0, MXML_DESCEND_FIRST);
  setStab(getNodeInteger(nd, 0) != 0);
  
  nd = mxmlFindElement(root_nd, root_nd, "curvature", 0, 0, MXML_DESCEND_FIRST);
  setCAngle(getNodeReal(nd, getCAngle()));
  
  nd = mxmlFindElement(root_nd, root_nd, "nint1", 0, 0, MXML_DESCEND_FIRST);
  setNInt1(getNodeInteger(nd, 12));
  
  nd = mxmlFindElement(root_nd, root_nd, "nint2", 0, 0, MXML_DESCEND_FIRST);
  setNInt2(getNodeInteger(nd, 12));
  
  nd = mxmlFindElement(root_nd, root_nd, "ndeg1", 0, 0, MXML_DESCEND_FIRST);
  setNDeg1(getNodeInteger(nd, 4));
  
  nd = mxmlFindElement(root_nd, root_nd, "ndeg2", 0, 0, MXML_DESCEND_FIRST);
  setNDeg2(getNodeInteger(nd, 4));
  
  nd = mxmlFindElement(root_nd, root_nd, "steps", 0, 0, MXML_DESCEND_FIRST);
  setSteps(getNodeInteger(nd, 100));
  
  nd = mxmlFindElement(root_nd, root_nd, "iad", 0, 0, MXML_DESCEND_FIRST);
  setIad(getNodeInteger(nd, 3));
  
  nd = mxmlFindElement(root_nd, root_nd, "npr", 0, 0, MXML_DESCEND_FIRST);
  setNPr(getNodeInteger(nd, 10));
  P_ERROR_X1(getNPr() > 0, "NPR must be greater than zero." );
  
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
  setNItC(getNodeInteger(nd, 5));
  
  nd = mxmlFindElement(root_nd, root_nd, "nitr", 0, 0, MXML_DESCEND_FIRST);
  setNItR(getNodeInteger(nd, 12));
  
  nd = mxmlFindElement(root_nd, root_nd, "nitk", 0, 0, MXML_DESCEND_FIRST);
  setNItK(getNodeInteger(nd, 12));

  nd = mxmlFindElement(root_nd, root_nd, "nderi", 0, 0, MXML_DESCEND_FIRST);
  setNDeri(getNodeInteger(nd, 8)); // this is just a number, higher then sys_nderi() provides

  nd = mxmlFindElement(root_nd, root_nd, "nsym", 0, 0, MXML_DESCEND_FIRST);
  setSymReSize(getNodeInteger(nd, 0));
  setSymImSize(getNodeInteger(nd, 0));

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
  mxmlDelete(tree);
}

void KNConstants::loadXmlFileOLD(const std::string &fileName)
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
  setFromType(static_cast<BifType>(getNodeInteger(nd, BifNone)));

  nd = mxmlFindElement(root_nd, root_nd, "label", 0, 0, MXML_DESCEND_FIRST);
  setLabel(getNodeInteger(nd));
  
  nd = mxmlFindElement(root_nd, root_nd, "pointtype", 0, 0, MXML_DESCEND_FIRST);
  setPointType(static_cast<PtType>(getNodeInteger(nd)));
    
  nd_a = mxmlFindElement(root_nd, root_nd, "cptype", 0, 0, MXML_DESCEND_FIRST);
  const char cp_type = getNodeChar(nd_a); 
  
  nd_b = mxmlFindElement(root_nd, root_nd, "cpnum", 0, 0, MXML_DESCEND_FIRST);
  const char cp_num = getNodeInteger(nd_b);
  
  setCp(varFromTypeNum(cp_type,cp_num));
  
  nd = mxmlFindElement(root_nd, root_nd, "switch", 0, 0, MXML_DESCEND_FIRST);
  setBranchSW(static_cast<BranchSW>(getNodeInteger(nd, NOSwitch)));

  if (getPointType() != SolUser)
  {
    nd = mxmlFindElement(root_nd, root_nd, "nparx", 0, 0, MXML_DESCEND_FIRST);
    setParxSize(getNodeInteger(nd));
    if (getParxSize() != 0)
    {
      mxml_node_t* parx_nd = mxmlFindElement(root_nd, root_nd, "parx", 0, 0, MXML_DESCEND_FIRST);
      
      int it = 0;
      for (nd = mxmlFindElement(parx_nd, parx_nd, "par", 0, 0, MXML_DESCEND_FIRST);
           nd != 0;
           nd = mxmlFindElement(nd, parx_nd->child, "par", 0, 0, MXML_NO_DESCEND) )
      {
        mxml_node_t* nd_type = mxmlFindElement(nd, nd, "type", 0, 0, MXML_DESCEND_FIRST);
        const char type = getNodeChar(nd_type);
        mxml_node_t* nd_num = mxmlFindElement(nd, nd, "num", 0, 0, MXML_DESCEND_FIRST);
        const char num = getNodeInteger(nd_num);
        setParx(it, varFromTypeNum(type,num));
        ++it;
      }
    } 
  } else
  {
    nd = mxmlFindElement(root_nd, root_nd, "neqns", 0, 0, MXML_DESCEND_FIRST);
    setEqnsSize(getNodeInteger(nd));
    setVarsSize(getNodeInteger(nd));
    if (getEqnsSize() != 0)
    {
      mxml_node_t* eqns_nd = mxmlFindElement(root_nd, root_nd, "eqns", 0, 0, MXML_DESCEND_FIRST);
      int it = 0;
      for (nd = mxmlFindElement(eqns_nd, eqns_nd, "eqn", 0, 0, MXML_DESCEND_FIRST);
           nd != 0;
           nd = mxmlFindElement(nd, eqns_nd->child, "eqn", 0, 0, MXML_NO_DESCEND) )
      {
        mxml_node_t* nd_type = mxmlFindElement(nd, nd, "type", 0, 0, MXML_DESCEND_FIRST);
        const char type = getNodeChar(nd_type);
        mxml_node_t* nd_num = mxmlFindElement(nd, nd, "num", 0, 0, MXML_DESCEND_FIRST);
        const char num = getNodeInteger(nd_num);
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
        const char type = getNodeChar(nd_type);
        mxml_node_t* nd_num = mxmlFindElement(nd, nd, "num", 0, 0, MXML_DESCEND_FIRST);
        const char num = getNodeInteger(nd_num);
        setVars(it, varFromTypeNum(type,num));
        ++it;
      }
    }
  }
  
  nd = mxmlFindElement(root_nd, root_nd, "nint", 0, 0, MXML_DESCEND_FIRST);
  setNInt(getNodeInteger(nd, 20));

  nd = mxmlFindElement(root_nd, root_nd, "ndeg", 0, 0, MXML_DESCEND_FIRST);
  setNDeg(getNodeInteger(nd, 5));
  
  nd = mxmlFindElement(root_nd, root_nd, "nmul", 0, 0, MXML_DESCEND_FIRST);
  setNMul(getNodeInteger(nd, 5));

  nd = mxmlFindElement(root_nd, root_nd, "stab", 0, 0, MXML_DESCEND_FIRST);
  setStab(getNodeInteger(nd, 0) != 0);
  
  nd = mxmlFindElement(root_nd, root_nd, "curvature", 0, 0, MXML_DESCEND_FIRST);
  setCAngle(getNodeReal(nd, getCAngle()));
  
  nd = mxmlFindElement(root_nd, root_nd, "nint1", 0, 0, MXML_DESCEND_FIRST);
  setNInt1(getNodeInteger(nd, 12));
  
  nd = mxmlFindElement(root_nd, root_nd, "nint2", 0, 0, MXML_DESCEND_FIRST);
  setNInt2(getNodeInteger(nd, 12));
  
  nd = mxmlFindElement(root_nd, root_nd, "ndeg1", 0, 0, MXML_DESCEND_FIRST);
  setNDeg1(getNodeInteger(nd, 4));
  
  nd = mxmlFindElement(root_nd, root_nd, "ndeg2", 0, 0, MXML_DESCEND_FIRST);
  setNDeg2(getNodeInteger(nd, 4));
  
  nd = mxmlFindElement(root_nd, root_nd, "steps", 0, 0, MXML_DESCEND_FIRST);
  setSteps(getNodeInteger(nd, 100));
  
  nd = mxmlFindElement(root_nd, root_nd, "iad", 0, 0, MXML_DESCEND_FIRST);
  setIad(getNodeInteger(nd, 3));
  
  nd = mxmlFindElement(root_nd, root_nd, "npr", 0, 0, MXML_DESCEND_FIRST);
  setNPr(getNodeInteger(nd, 10));
  P_ERROR_X1(getNPr() > 0, "NPR must be greater than zero." );
  
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
  setNItC(getNodeInteger(nd, 5));
  
  nd = mxmlFindElement(root_nd, root_nd, "nitr", 0, 0, MXML_DESCEND_FIRST);
  setNItR(getNodeInteger(nd, 12));
  
  nd = mxmlFindElement(root_nd, root_nd, "nitk", 0, 0, MXML_DESCEND_FIRST);
  setNItK(getNodeInteger(nd, 12));

  nd = mxmlFindElement(root_nd, root_nd, "nderi", 0, 0, MXML_DESCEND_FIRST);
  setNDeri(getNodeInteger(nd, 8)); // this is just a number, higher then sys_nderi() provides

  nd = mxmlFindElement(root_nd, root_nd, "nsym", 0, 0, MXML_DESCEND_FIRST);
  setSymReSize(getNodeInteger(nd, 0));
  setSymImSize(getNodeInteger(nd, 0));

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
  mxmlDelete(tree);
}

void KNConstants::saveXmlFile(const std::string &fileName)
{
  std::ofstream file(fileName.c_str());
  printXmlFile(file);
}

static inline void mxmlNewDouble(mxml_node_t *node, double dbl)
{
  mxmlNewTextf(node, 0, "%.15E", dbl);
}

void KNConstants::printXmlFile(std::ostream& file)
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
  
  node = mxmlNewElement(data, "fromtype");
  mxmlNewText(node, 0, BifTypeTable.TypeToCode(getFromType()).c_str());
  
  node = mxmlNewElement(data, "label");
  mxmlNewInteger(node, getLabel());
  
  node = mxmlNewElement(data, "pointtype");
  mxmlNewText(node, 0, PtTypeTable.TypeToCode(getPointType()).c_str());
  
  node = mxmlNewElement(data, "cp");
  mxmlNewText(node, 0, VarTable.TypeToCode(getCp()).c_str());

  node = mxmlNewElement(data, "switch");
  mxmlNewText(node, 0, BranchSWTable.TypeToCode(getBranchSW()).c_str());
 
  if (getPointType() != SolUser)
  {
    node = mxmlNewElement(data, "nparx");
    mxmlNewInteger(node, getParxSize());

    if (getParxSize() != 0)
    {
      mxml_node_t *group_parx = mxmlNewElement(data, "parx");
      for (int i = 0; i < getParxSize(); ++i)
      {
        node = mxmlNewElement(group_parx, "par");
        mxmlNewText(node, 0, VarTable.TypeToCode(getParx(i)).c_str());
      }
    }
  }
  else
  {
    node = mxmlNewElement(data, "neqns");
    mxmlNewInteger(node, getEqnsSize());

    if (getEqnsSize() != 0)
    {
      mxml_node_t *group_eqns = mxmlNewElement(data, "eqns");
      mxml_node_t *group_vars = mxmlNewElement(data, "vars");

      for (int i = 0; i < getEqnsSize(); ++i)
      {
        node = mxmlNewElement(group_eqns, "eqn");
        mxmlNewText(node, 0, EqnTable.TypeToCode(getEqns(i)).c_str());
      }
      for (int i = 0; i < getVarsSize(); ++i)
      {        
        node = mxmlNewElement(group_vars, "var");
        mxmlNewText(node, 0, VarTable.TypeToCode(getVars(i)).c_str());
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
  
  node = mxmlNewElement(data, "curvature");
  mxmlNewDouble(node, getCAngle());
  
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
  
  node = mxmlNewElement(data, "npr");
  mxmlNewInteger(node, getNPr());
  
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
  mxmlNewInteger(node, getNItC());
  
  node = mxmlNewElement(data, "nitr");
  mxmlNewInteger(node, getNItR());
  
  node = mxmlNewElement(data, "nitk");
  mxmlNewInteger(node, getNItK());
  
  node = mxmlNewElement(data, "nderi");
  mxmlNewInteger(node, getNDeri());

  node = mxmlNewElement(data, "nsym");
  mxmlNewInteger(node, getSymReSize());

  if (getSymReSize() != 0)
  {
    mxml_node_t *group_sym = mxmlNewElement(data, "sym");
    for (int i = 0; i < getSymReSize(); ++i)
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

bool KNConstants::toEqnVar(KNSystem& sys,
                          KNArray1D<Eqn>& eqn, KNArray1D<Var>& var,                 // input
                          KNArray1D<Eqn>& eqn_refine, KNArray1D<Var>& var_refine,   // output
                          KNArray1D<Eqn>& eqn_start, KNArray1D<Var>& var_start, bool& findangle)
{
  // initializing the equations and variables
  if (getPointType() == SolUser)
  {
    eqn.init(getEqnsSize());
    var.init(getVarsSize());
    for (int i = 0; i < getEqnsSize(); i++)
    {
      eqn(i) = getEqns(i);
      var(i) = getVars(i);
    }
  }
  else
  {
    KNArray1D<Var> L_PARX(getParxSize());
    for (int i = 0; i < getParxSize(); i++)
    {
      L_PARX(i) = getParx(i);
    }
    setBranchSW(PtToEqnVar(eqn, var, getPointType(), L_PARX, sys.npar()));
  }

  // checking whether it is an autonomous problem or not
  bool aut = false;
  bool phaseRot = false;
  bool needTF = false;
  findangle = false;
  for (int i = 0; i < eqn.size(); i++)
  {
    if ((eqn(i) == EqnPhase) || (eqn(i) == EqnTORPhase0)) aut = true;
    if (eqn(i) == EqnPhaseRot) phaseRot = true;
    if ((eqn(i) == EqnTFPD) ||
        (eqn(i) == EqnTFLP) ||
        (eqn(i) == EqnTFLPAUT) ||
        (eqn(i) == EqnTFLPAUTROT) ||
        (eqn(i) == EqnTFCPLX_RE))
    {
      needTF = true;
    }
    if (eqn(i) == EqnTFCPLX_RE) findangle = true;
  }

  // setting up for refinement
  if (aut)
  {
    if (phaseRot)
    {
      // std::cout<<"Phase and PhaseRot\n";
      eqn_refine.init(3);
      var_refine.init(3);
      eqn_refine(0) = eqn(0);
      eqn_refine(1) = EqnPhase;
      eqn_refine(2) = EqnPhaseRot;
      var_refine(0) = var(0);
      var_refine(1) = var(var.size() - 2);
      var_refine(2) = var(var.size() - 1);
    }
    else
    {
      // std::cout<<"Phase\n";
      eqn_refine.init(2);
      var_refine.init(2);
      eqn_refine(0) = eqn(0);
      eqn_refine(1) = EqnPhase;
      var_refine(0) = var(0);
      var_refine(1) = var(var.size() - 1);
    }
  }
  else
  {
    if (phaseRot)
    {
      // this happens when a steady state solution of a laser is computed
      eqn_refine.init(2);
      var_refine.init(2);
      eqn_refine(0) = eqn(0);
      eqn_refine(1) = EqnPhaseRot;
      var_refine(0) = var(0);
      var_refine(1) = var(var.size() - 1);
    }
    else
    {
      eqn_refine.init(1);
      var_refine.init(1);
      eqn_refine(0) = eqn(0);
      var_refine(0) = var(0);
    }
  }

  if (getBranchSW() == TFHBSwitch)
  {
      KNArray1D<Eqn> ee(eqn_refine);
      KNArray1D<Var> vv(var_refine);
      eqn_refine.init(ee.size()-1);
      var_refine.init(ee.size()-1);
      for (int i = 0, j = 0; i < ee.size(); ++i)
      {
        if (ee(i) != EqnPhase) { eqn_refine(j) = ee(i); var_refine(j) = vv(j); ++j; }
      }
  }
  
  // Swithing from a Neimark-Sacker to a Torus (Delay Equation)
  if (getBranchSW() == TFTRSwitch)
  {
    eqn_refine(0) = EqnSol;
    var_refine(0) = VarSol;
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
      eqn_start.init(eqn_refine.size() + 1);
      var_start.init(var_refine.size() + 1);
      eqn_start(0) = eqn_refine(0);
      var_start(0) = var_refine(0);
      eqn_start(1) = eqn_temp;
      var_start(var_refine.size()) = getCp();
      for (int i = 1; i < eqn_refine.size(); i++)
      {
        eqn_start(i + 1) = eqn_refine(i);
        var_start(i) = var_refine(i);
      }
      needTF = true;
      break;
    case TFTRSwitch:
    case TFHBSwitch:
      eqn_start.init(eqn_refine.size() + 2);
      var_start.init(var_refine.size() + 2);
      eqn_start(0) = eqn_refine(0);
      var_start(0) = var_refine(0);
      eqn_start(1) = EqnTFCPLX_RE;
      eqn_start(2) = EqnTFCPLX_IM;
      var_start(1) = VarAngle; // VarToIndex(VarAngle,sys.npar()); // CH
      var_start(var_refine.size() + 1) = getCp();
      for (int i = 1; i < eqn_refine.size(); i++)
      {
        eqn_start(i + 2) = eqn_refine(i);
        var_start(i + 1) = var_refine(i);
      }
      needTF = true;
      findangle = true;
      break;
    default:
      eqn_start.init(eqn.size());
      var_start.init(var.size());
      eqn_start = eqn;
      var_start = var;
      break;
  }
  return needTF;
}

void KNConstants::initDimensions(const KNSystem* sys)
{
  initParNames(sys->npar());
  const char **names = new const char *[sys->npar()+1];
  for (int i=0; i<sys->npar()+1; ++i) names[i] = 0;
  sys->parnames(names);
  for (int i=0; i<sys->npar(); ++i) if (names[i] != 0) parNames[i] = names[i];
  delete[] names;
  // this will send a signal so parNames must be ready
  setNPar(sys->npar());
  setNDim(sys->ndim());
  VarTable.update();
  // just to send signals about updated VarTable
  setNPar(sys->npar());
}

void KNConstants::setSysNameText(const std::string& str, bool testing)
{
  setSysName(str);
  KNSystem* sys = 0;
  try
  {
    sys = new KNSystem(getSysName(), getNDeri());
    initDimensions(sys);
    delete sys;
  }
  catch (KNException ex)
  {
    delete sys;
    setNPar(0);
    setNDim(0);
    throw(ex);
  }
}

template<> const TypeTuple<PtType> KNConstants::TypeTupleTab<PtType>::tabStatic[] = {
  {0,  SolUser,      "USER",             "User defined"}, // User
  {1,  SolODE,       "ODE",              "ODE limit cycle"}, // ODE limit cycle
  {2,  SolAUTODE,    "ODE_AUT",          "ODE limit cycle (aut)" }, // ODE limit cycle (aut)
  {3,  SolTF,        "DDE_LC",           "Limit cycle"}, // "Limit cycle
  {4,  BifTFLP,      "DDE_LP",           "Limit point"}, // Limit point
  {5,  BifTFPD,      "DDE_PD",           "Period doubling"}, // Period doubling
  {6,  BifTFNS,      "DDE_NS",           "Neimark-Sacker"}, // Neimark-Sacker
  {7,  SolTFBRSW,    "DDE_BP_SW",        "Branch switch"}, // Branch switch
  {8,  SolTFPDSW,    "DDE_PD_SW",        "Period doubling switch"}, // Period doubling switch
  {9,  SolTFAUT,     "DDE_AUT",          "Limit cycle (aut)"}, // Limit cycle (aut)
  {10, BifTFAUTLP,   "DDE_AUT_LP",       "Limit point (aut)"}, //
  {11, BifTFAUTPD,   "DDE_AUT_PD",       "Period doubling (aut)"},
  {12, BifTFAUTNS,   "DDE_AUT_NS",       "Neimark-Sacker (aut)"},
  {13, SolTFAUTBRSW, "DDE_AUT_BP_SW",    "Branch switch (aut)"},
  {14, SolTFAUTPDSW, "DDE_AUT_PD_SW",    "Period doubling switch (aut)"},
  {15, SolTFAUTHBSW, "DDE_AUT_HOPF_SW",  "Hopf switch (aut)"},
  {16, SolTor,       "DDE_TORUS",        "Torus"},
  {17, SolTorNS,     "DDE_TORUS_SW",     "Torus from NS"},
  {18, SolAUTTor,    "DDE_AUT_TORUS",    "Torus (aut)"},
  {19, SolAUTTorNS,  "DDE_AUT_TORUS_SW", "Torus from NS (aut)"},
  {-1, SolUser, "", ""}};

template<> const TypeTuple<Eqn> KNConstants::TypeTupleTab<Eqn>::tabStatic[] = {
  {0,  EqnNone,       "NX",                "None" },
  {1,  EqnSol,        "DDE_LC",            "Limit cycle"},
  {2,  EqnODESol,     "ODE_LC",            "ODE Limit cycle"},
  {3,  EqnPhase,      "PHASE",             "PHCND (transl)"},
  {4,  EqnPhaseRot,   "PHASE_ROT",         "PHCND (rot)"},
  {5,  EqnTORSol,     "DDE_TORUS",         "Torus"},
  {6,  EqnTORPhase0,  "TORUS_PHS1",        "PHCND 1 (torus)"},
  {7,  EqnTORPhase1,  "TORUS_PHS2",        "PHCND 2 (torus)"},
  {8,  EqnTFLP,       "DDE_TF_LP",         "LP TF"},
  {9,  EqnTFPD,       "DDE_TF_PD",         "PD TF"},
  {10, EqnTFLPAUT,    "DDE_AUT_TF_LP",     "LP TF (aut)"},
  {11, EqnTFLPAUTROT, "DDE_AUT_ROT_TF_LP", "LP TF (sym)"},
  {12, EqnTFCPLX_RE,  "DDE_TF_NS_RE",      "NS TF (real)"},
  {13, EqnTFCPLX_IM,  "DDE_TF_NS_IM",      "NS TF (imag)"},
  {-1, EqnNone, "", ""}};
  
template<> const TypeTuple<BranchSW> KNConstants::TypeTupleTab<BranchSW>::tabStatic[] = {
  {0, NOSwitch,         "NX",         "No switch"},
  {1, TFBRSwitch,       "BP",         "Branch"},
  {2, TFPDSwitch,       "PD",         "Period doubling"},
  {3, TFHBSwitch,       "HOPF",       "Hopf (aut)"},
  {4, TFTRSwitch,       "TORUS",      "Torus"},
  {5, TFBRAUTSwitch,    "AUT_BP",     "Branch (aut)"},
  {6, TFBRAUTROTSwitch, "AUT_ROT_BP", "Branch (sym)"},
  {-1, NOSwitch, "", ""}};

template<> const TypeTuple<Var> KNConstants::TypeTupleTab<Var>::tabStatic[] = {
  {0, VarNone,   "NX",          "None"},
  {1, VarSol,    "DDE_LC",      "Limit cycle"},
  {2, VarODESol, "ODE_LC",      "ODE limit cycle"},
  {3, VarTORSol, "DDE_TORUS",   "Torus"},
  {4, VarPAR0,   "PAR0",        "Par 0"},
  {5, VarAngle,  "PAR_ANGLE",   "Angle"},
  {6, VarPeriod, "PAR_PERIOD",  "Period"},
  {7, VarRot,    "PAR_ROT_NUM", "Rot. num."},
  {-1, VarEnd, "", ""}};

template<> const TypeTuple<BifType> KNConstants::TypeTupleTab<BifType>::tabStatic[] = {
  {0, BifNone,          "NX", "None"},
  {1, BifLP,            "LP", "LP"},
  {2, BifPD,            "PD", "PD"},
  {3, BifNS,            "NS", "NS"},
  {4, BifNS,            "UN", "UN"},
  {5, BifMax,           "MX", "MX"},
  {6, BifEndPoint,      "EP", "EP"},
  {7, BifNoConvergence, "NC", "NC"},
  {-1, BifNone, "", ""}};
  
// specialization to Var

template<> void KNConstants::TypeTupleTab<Var>::update()
{
  size_t sz = 0;
  while (tabStatic[sz].index >= 0) ++sz;
  tab.resize(sz + parent->getNPar());

  size_t k = 0;
  for (; tabStatic[k].type != VarPAR0; ++k)
  {
    tab[k].index = tabStatic[k].index;
    tab[k].type = tabStatic[k].type;
    tab[k].code = tabStatic[k].code;
    tab[k].name = tabStatic[k].name;
  }
  for (int j = 0; j < parent->getNPar(); ++j)
  {
    tab[k+j].index = k+j;
    tab[k+j].type = static_cast<Var>(VarPAR0 + j);
    tab[k+j].code = parent->parNames.at(j);
    tab[k+j].name = parent->parNames.at(j);
  }
  k++;
  for (; k<sz; ++k)
  {
    tab[k + parent->getNPar() - 1].index = k + parent->getNPar() - 1;
    tab[k + parent->getNPar() - 1].type = tabStatic[k].type;
    tab[k + parent->getNPar() - 1].code = tabStatic[k].code;
    tab[k + parent->getNPar() - 1].name = tabStatic[k].name;
  }
}
