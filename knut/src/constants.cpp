// ------------------------------------------------------------------------- //
//
// This is part of KNUT
// Copyright (c) 2006 by Robert Szalai
//
// For license, see the file COPYING in the root directory of the package
//
// ------------------------------------------------------------------------- //

#include "config.h"
#include "constants.h"
#include "matrix.h"
#include <iostream>
#include <fstream>
#include <sstream>

// For loading XML constant files
#include "mxml.h"

class KNXmlConstantsReader : public KNAbstractConstantsReader
{
  public:
    KNXmlConstantsReader();
    KNXmlConstantsReader(mxml_node_t * tree);
    ~KNXmlConstantsReader();
    bool getSystem(std::string& strout) override;
    bool getTextField(const char* field, std::string& strout) override;
    bool getIndexField(const char* field, size_t& out) override;
    bool getDoubleField(const char* field, double& out) override;
    bool getStringListField(const char* field, std::vector<std::string>& out) override;
    bool getIndexListField(const char* field, std::vector<size_t>& out) override;
  private:
    const std::string filename;
    mxml_node_t* tree;
    mxml_node_t* knut_node;
    mxml_node_t* branch_node;
};

static int getXmlTree(const std::string& filename, mxml_node_t ** ptree)
{
  FILE *fp;
  mxml_node_t *tree;

  fp = fopen(filename.c_str(), "r");
  if (fp == nullptr) return -1;
  tree = mxmlLoadFile(nullptr, fp, MXML_OPAQUE_CALLBACK);
  fclose(fp);

  mxml_node_t* root_nd = mxmlFindElement(tree, tree, "knut", nullptr, nullptr, MXML_DESCEND_FIRST);
  if (!root_nd)
  {
    root_nd = mxmlFindElement(tree, tree, "pdde", nullptr, nullptr, MXML_DESCEND_FIRST);
  }
  int ver = -1;
  const char *attr = mxmlElementGetAttr(root_nd, "version");
  if (attr) ver = strtol(attr, nullptr, 10);
  *ptree = tree;
  return ver;  
}

KNXmlConstantsReader::KNXmlConstantsReader(mxml_node_t * xtree) : tree(xtree), knut_node(nullptr), branch_node(nullptr)
{
  knut_node = mxmlFindElement(xtree, xtree, "knut", nullptr, nullptr, MXML_DESCEND_FIRST);
  if (knut_node == nullptr)
  {
    knut_node = mxmlFindElement(xtree, xtree, "pdde", nullptr, nullptr, MXML_DESCEND_FIRST);
  }
  branch_node = mxmlFindElement(knut_node, knut_node, "branch", nullptr, nullptr, MXML_DESCEND_FIRST);
}

KNXmlConstantsReader::~KNXmlConstantsReader()
{
  mxmlDelete(tree);
}

class KNXmlConstantsWriter : public KNAbstractConstantsWriter
{
  public:
    KNXmlConstantsWriter();
    KNXmlConstantsWriter(const std::string& filename);
    ~KNXmlConstantsWriter();
    void setSystem(const std::string& str) override;
    void setTextField(const char* fieldname, const std::string& str) override;
    void setDoubleField(const char* fieldname, double val) override;
    void setIndexField(const char* fieldname, size_t val) override;
    void setStringListField(const char* field, const std::vector<std::string>& vec) override;
    void setIndexListField(const char* field, const std::vector<size_t>& vec) override;
  private:
    const std::string filename;
    mxml_node_t *xml;
    mxml_node_t *knut_node;
    mxml_node_t *system_node;
    mxml_node_t *branch_node;
};

void KNConstantNames::addConstantName(const char * type, const char * name, void * var, void(*nsetFun)())
{ 
  if (nlist.find (name) == nlist.end())
  {
    tuple_t t = { type, var, nsetFun };
    nlist[name] = t;
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
  } else return nullptr;
}

const char *knut_whitespace_cb(mxml_node_t *node, int where)
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
    else return nullptr;
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
  } else if ((node->child == nullptr)&&(where == MXML_WS_AFTER_OPEN))
  {
    return "\n";
  }
  return (nullptr);
}

// static inline bool inputAssert(std::istream& is)
// {
//   if (is.rdstate() & (std::istream::eofbit | std::istream::failbit | std::istream::badbit))
//   {
//     switch (is.rdstate())
//     {
//       case std::istream::eofbit:
//         P_MESSAGE1("std::istream: unexpected end of file.");
//         return true;
//         break;
//       case std::istream::failbit:
//         P_MESSAGE1("std::istream: Input failed.");
//         return true;
//         break;
//       case std::istream::badbit:
//         P_MESSAGE1("std::istream: bad input.");
//         return true;
//         break;
//       default:
//         P_MESSAGE1("std::istream: unexpected error.");
//         return true;
//         break;
//     }
//   }
//   return false;
// }
// 
// static inline const char* c2s(char *buf, char c)
// {
//   buf[0] = c;
//   buf[1] = '\0';
//   return buf;
// }

static void commaToVector(const std::string& str, std::vector<std::string>& list)
{
  list.clear();
  size_t pos = 0, npos;
  do {
    npos = str.find(',', pos);
    if (npos == std::string::npos) npos = str.size();
    if (!str.substr(pos,npos-pos).empty())
    { 
      list.push_back(str.substr(pos,npos-pos));
//      std::cout << list[list.size()-1] << "\n";
    }
    pos = npos+1;
  } while (npos != str.size());
}

static void commaToVector(const std::string& str, std::vector<size_t>& list)
{
  list.clear();
  size_t pos = 0, npos;
  do {
    npos = str.find(',', pos);
    if (npos == std::string::npos) npos = str.size();
    if (!str.substr(pos,npos-pos).empty())
    { 
      list.push_back(toSizeT(strtol(str.substr(pos,npos-pos).c_str(), nullptr, 10)));
//      std::cout << list[list.size()-1] << "\n";
    }
    pos = npos+1;
  } while (npos != str.size());
}

void KNConstantsBase::loadXmlFileV5(const std::string &fileName)
{
  mxml_node_t *tree;
  int ver = getXmlTree(fileName, &tree);

  // compatibility options
  if ((ver < 5) && (ver > 2))
  {
    loadXmlFileV4(fileName);
    return;
  } else if (ver < 3)
  { 
    loadXmlFileV2(fileName);
    return;
  }
  // continuing normal
  KNXmlConstantsReader reader(tree);
  read(reader);
}

void KNConstantsBase::saveXmlFileV5(const std::string &fileName)
{
  KNXmlConstantsWriter writer(fileName);
  write(writer);
}

void KNConstantsBase::write(KNAbstractConstantsWriter& writer)
{
  writer.setSystem (getSystem());
  writer.setTextField ("input", getInputFile());
  writer.setTextField ("output", getOutputFile());
  if (!getSysName().empty()) writer.setTextField ("sysname", getSysName());
  if (!getSysType().empty()) writer.setTextField ("systype", getSysType());
  writer.setTextField ("fromtype", BifTypeTable.TypeToCode(getFromType()));
  writer.setIndexField ("label", getLabel());
  writer.setTextField ("pointtype", PtTypeTable.TypeToCode(getPointType()));
  writer.setTextField ("cp", VarTable.TypeToCode(getCp()));
  writer.setTextField ("switch", BranchSWTable.TypeToCode(getBranchSW()));
  
  std::vector<std::string> parx_list;
  VarTable.TypeToCodeVector (parx_list, getParxVector());
  writer.setStringListField("parx", parx_list);

  std::vector<std::string> eqns_list;
  EqnTable.TypeToCodeVector (eqns_list, getEqnsVector());
  writer.setStringListField("eqns", eqns_list);

  std::vector<std::string> vars_list;
  VarTable.TypeToCodeVector (vars_list, getVarsVector());
  writer.setStringListField("vars", vars_list);
    
  writer.setIndexField ("nint", getNInt());
  writer.setIndexField ("ndeg", getNDeg());
  writer.setIndexField ("nmul", getNMul());
  writer.setIndexField ("stab", getStab());
  writer.setDoubleField ("curvature", getCAngle());
  writer.setIndexField ("nint1", getNInt1());
  writer.setIndexField ("nint2", getNInt2());
  writer.setIndexField ("steps", getSteps());
  writer.setIndexField ("iad", getIad());
  writer.setIndexField ("npr", getNPr());
  writer.setIndexField ("ssm", getSSM());
  writer.setDoubleField ("cpmin", getCpMin());
  writer.setDoubleField ("cpmax", getCpMax());
  writer.setDoubleField ("ds", getDs());
  writer.setDoubleField ("dsmin", getDsMin());
  writer.setDoubleField ("dsmax", getDsMax());
  writer.setDoubleField ("dsstart", getDsStart());
  writer.setDoubleField ("epsc", getEpsC());
  writer.setDoubleField ("epsr", getEpsR());
  writer.setDoubleField ("epsk", getEpsK());
  writer.setIndexField ("nitc", getNItC());
  writer.setIndexField ("nitr", getNItR());
  writer.setIndexField ("nitk", getNItK());
  writer.setIndexField ("nderi", getNDeri());
  
  writer.setIndexListField("symre", getSymReVector());
  writer.setIndexListField("symim", getSymImVector());
}

void KNConstantsBase::read(KNAbstractConstantsReader& reader)
{
  std::string buf;
  size_t idx;
  double dbl;

  if (reader.getSystem(buf)) if (!buf.empty()) setSystemText (buf);
  if (reader.getTextField("input", buf)) setInputFile(buf);
  if (reader.getTextField("output", buf)) setOutputFile(buf);
  if (reader.getTextField("sysname", buf)) if (!buf.empty()) setSysNameText(buf);
  if (reader.getTextField("systype", buf)) setSysType(buf);
  if (reader.getTextField("fromtype", buf)) setFromType(BifTypeTable.CodeToType(buf.c_str()));
  if (reader.getIndexField("label", idx)) setLabel(idx);
  if (reader.getTextField("pointtype", buf)) setPointType(PtTypeTable.CodeToType(buf.c_str()));
  if (reader.getTextField("cp", buf)) setCp(VarTable.CodeToType(buf.c_str()));
  if (reader.getTextField("switch", buf)) setBranchSW(BranchSWTable.CodeToType(buf.c_str()));

  std::vector<std::string> parx_list;
  std::vector<Var> parx_tp;
  reader.getStringListField("parx", parx_list);
  VarTable.CodeToTypeVector(parx_tp, parx_list);
  setParxVector(parx_tp);

  std::vector<std::string> eqns_list;
  std::vector<Eqn> eqns_tp;
  reader.getStringListField("eqns", eqns_list);
  EqnTable.CodeToTypeVector(eqns_tp, eqns_list);
  setEqnsVector(eqns_tp);
  
  std::vector<std::string> vars_list;
  std::vector<Var> vars_tp;
  reader.getStringListField("vars", vars_list);
  VarTable.CodeToTypeVector(vars_tp, vars_list);
  setVarsVector(vars_tp);

  if (reader.getIndexField("nint", idx)) setNInt(idx);
  if (reader.getIndexField("ndeg", idx)) setNDeg(idx);
  if (reader.getIndexField("nmul", idx)) setNMul(idx);
  if (reader.getIndexField("stab", idx)) setStab(idx);
  if (reader.getDoubleField("curvature", dbl)) setCAngle(dbl);
  if (reader.getIndexField("nint1", idx)) setNInt1(idx);
  if (reader.getIndexField("nint2", idx)) setNInt2(idx);
  if (reader.getIndexField("ndeg1", idx)) setNDeg1(idx);
  if (reader.getIndexField("ndeg2", idx)) setNDeg2(idx);
  if (reader.getIndexField("steps", idx)) setSteps(idx);
  if (reader.getIndexField("iad", idx)) setIad(idx);
  if (reader.getIndexField("npr", idx)) setNPr(idx);
  if (reader.getIndexField("ssm", idx)) setSSM(idx);
  if (reader.getDoubleField("cpmin", dbl)) setCpMin(dbl);
  if (reader.getDoubleField("cpmax", dbl)) setCpMax(dbl);
  if (reader.getDoubleField("ds", dbl)) setDs(dbl);
  if (reader.getDoubleField("dsmin", dbl)) setDsMin(dbl);
  if (reader.getDoubleField("dsmax", dbl)) setDsMax(dbl);
  if (reader.getDoubleField("dsstart", dbl)) setDsStart(dbl);
  if (reader.getDoubleField("epsc", dbl)) setEpsC(dbl);
  if (reader.getDoubleField("epsr", dbl)) setEpsR(dbl);
  if (reader.getDoubleField("epsk", dbl)) setEpsK(dbl);
  if (reader.getIndexField("nitc", idx)) setNItC(idx);
  if (reader.getIndexField("nitr", idx)) setNItR(idx);
  if (reader.getIndexField("nitk", idx)) setNItK(idx);
  if (reader.getIndexField("nderi", idx)) setNDeri(idx);

  std::vector<size_t> symre_list, symim_list;
  reader.getIndexListField("symre", symre_list);
  reader.getIndexListField("symim", symim_list);
  setSymReVector(symre_list);
  setSymImVector(symim_list);
}

bool KNXmlConstantsReader::getSystem(std::string& strout)
{
  mxml_node_t* system_node = mxmlFindElement(knut_node, knut_node, "system", nullptr, nullptr, MXML_DESCEND_FIRST);
  if (system_node != nullptr)
  {
    const char* str = mxmlGetOpaque (system_node);
    if (str)
    {
      strout = str;
      return true;
    }
  }
  return false;
}

bool KNXmlConstantsReader::getTextField(const char* field, std::string& strout)
{
  const char* str = mxmlElementGetAttr(branch_node, field);
  if (str)
  {
    strout = str;
    return true;
  }
  return false;
}

bool KNXmlConstantsReader::getIndexField(const char* field, size_t& out)
{
  const char* str = mxmlElementGetAttr(branch_node, field);
  if (str)
  {
    out = toSizeT(strtol(str, nullptr, 10));
    return true;
  }
  return false;
}

bool KNXmlConstantsReader::getDoubleField(const char* field, double& out)
{
  const char* str = mxmlElementGetAttr(branch_node, field);
  if (str)
  {
    std::istringstream strstr(str);
    strstr >> out;
    return true;
  }
  return false;
}

bool KNXmlConstantsReader::getStringListField(const char* field, std::vector<std::string>& out)
{
  std::string str;
  if ( getTextField(field, str) )
  {
    commaToVector(str, out);
    return true;
  }
  return false;
}

bool KNXmlConstantsReader::getIndexListField(const char* field, std::vector<size_t>& out)
{
  std::string str;
  if ( getTextField(field, str) )
  {
    commaToVector(str, out);
    return true;
  }
  return false;
}

KNXmlConstantsWriter::KNXmlConstantsWriter(const std::string& fname) : filename(fname)
{
  xml = mxmlNewXML ("cfile 1.0");
  knut_node = mxmlNewElement (xml, "knut");
  mxmlElementSetAttr (knut_node, "version", PACKAGE_VERSION);
  system_node = mxmlNewElement (knut_node, "system");
  branch_node = mxmlNewElement (knut_node, "branch");  
}

KNXmlConstantsWriter::~KNXmlConstantsWriter()
{
  char *xmlString = mxmlSaveAllocString (xml, knut_whitespace_cb);
  if (xmlString)
  {
    if (!filename.empty()) 
    {
      std::ofstream file(filename.c_str());
      file << xmlString;
    } else
    {
      std::cout << xmlString;
    }
    free(xmlString);
  }
  mxmlDelete(xml);
}

void KNXmlConstantsWriter::setSystem (const std::string& str)
{
  mxmlNewOpaque (system_node, str.c_str());
}

void KNXmlConstantsWriter::setTextField(const char* fieldname, const std::string& str)
{
  mxmlElementSetAttr(branch_node, fieldname, str.c_str());
}

void KNXmlConstantsWriter::setDoubleField(const char* fieldname, double val)
{
  std::ostringstream strstr;
  strstr.precision(15);
  strstr << val;
  mxmlElementSetAttr (branch_node, fieldname, strstr.str().c_str());
}

void KNXmlConstantsWriter::setIndexField(const char* fieldname, size_t val)
{
  std::ostringstream strstr;
  strstr << val;
  mxmlElementSetAttr (branch_node, fieldname, strstr.str().c_str());
}

void KNXmlConstantsWriter::setStringListField(const char* field, const std::vector<std::string>& vec)
{
  if (vec.size() != 0)
  {
    std::string field_str;
    for (size_t i = 0; i < vec.size(); ++i)
    {
      if (i != 0) field_str.append(",");
      field_str.append(vec[i]);
    }
    setTextField(field, field_str); // --- this is a comma separated array
  }
}

void KNXmlConstantsWriter::setIndexListField(const char* field, const std::vector<size_t>& vec)
{
  if (vec.size() != 0)
  {
    std::ostringstream field_str;
    for (size_t i = 0; i < vec.size(); ++i)
    {
      if (i != 0) field_str << ",";
      field_str << vec[i];
    }
    setTextField(field, field_str.str().c_str()); // --- this is a comma separated array
  }
}

bool KNConstantsBase::toEqnVar(KNExprSystem& sys,
                          KNArray1D<Eqn>& eqn, KNArray1D<Var>& var,                 // input
                          KNArray1D<Eqn>& eqn_refine, KNArray1D<Var>& var_refine,   // output
                          KNArray1D<Eqn>& eqn_start, KNArray1D<Var>& var_start, bool& findangle)
{
  // initializing the equations and variables
  if (getPointType() == SolUser)
  {
    eqn.init(getEqnsSize());
    var.init(getVarsSize());
    for (size_t i = 0; i < getEqnsSize(); i++)
    {
      eqn(i) = getEqns(i);
      var(i) = getVars(i);
    }
  }
  else
  {
    KNArray1D<Var> L_PARX(getParxSize());
    for (size_t i = 0; i < getParxSize(); i++)
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
  for (size_t i = 0; i < eqn.size(); i++)
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

  if ((getBranchSW() == TFHBSwitch) || (getBranchSW() == TFSSMSwitch))
  {
      KNArray1D<Eqn> ee(eqn_refine);
      KNArray1D<Var> vv(var_refine);
      eqn_refine.init(ee.size()-1);
      var_refine.init(ee.size()-1);
      for (size_t i = 0, j = 0; i < ee.size(); ++i)
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
    case TFBRAUTSwitch:
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
      for (size_t i = 1; i < eqn_refine.size(); i++)
      {
        eqn_start(i + 1) = eqn_refine(i);
        var_start(i) = var_refine(i);
      }
      needTF = true;
      break;
    case TFTRSwitch:
    case TFHBSwitch:
    case TFSSMSwitch:
      eqn_start.init(eqn_refine.size() + 2);
      var_start.init(var_refine.size() + 2);
      eqn_start(0) = eqn_refine(0);
      var_start(0) = var_refine(0);
      eqn_start(1) = EqnTFCPLX_RE;
      eqn_start(2) = EqnTFCPLX_IM;
      var_start(1) = VarAngle; // VarToIndex(VarAngle,sys.npar()); // CH
      var_start(var_refine.size() + 1) = getCp();
      for (size_t i = 1; i < eqn_refine.size(); i++)
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

void KNConstantsBase::initDimensions(const KNExprSystem* sys)
{
  initParNames(sys->npar());
  std::vector<const char*> names(sys->npar()+1, nullptr);

  sys->parnames(names.data());
  for (size_t i=0; i<sys->npar(); ++i) if (names[i] != nullptr) parNames[i] = names[i];

  // this will send a signal so parNames must be ready
  setNPar(sys->npar());
  setNDim(sys->ndim());
  VarTable.update();
  // just to send signals about updated VarTable
  setNPar(sys->npar());
}

void KNConstantsBase::setSystemText (const std::string& str)
{
  setSystem (str);
  KNExprSystem* sys = nullptr;
  try
  {
    sys = new KNExprSystem (getSystem (), false);
    initDimensions (sys);
    delete sys;
    sys = nullptr;
  }
  catch (KNException& ex)
  {
    delete sys;
    sys = nullptr;
    setNPar (0);
    setNDim (0);
    ex.exprStr (getSystem ());
    throw (ex);
  }
}

void KNConstantsBase::setSysNameText(const std::string& str, bool /*testing*/)
{
  setSysName(str);
  KNSystem* sys = nullptr;
  try
  {
    sys = new KNSystem(getSysName(), false);
    initDimensions(sys);
    delete sys;
    sys = nullptr;
  }
  catch (KNException& ex)
  {
    delete sys;
    sys = nullptr;
    setNPar(0);
    setNDim(0);
    throw(ex);
  }
}
  
// specialization to Var
template<> void TypeTupleTab<Var>::update()
{
  size_t sz = 0;
  while (TypeTupleTabBase<Var>::tabStatic[sz].index < ~(size_t)0) ++sz;
  tab.resize(sz + parent->getNPar());

  size_t k = 0;
  for (; TypeTupleTabBase<Var>::tabStatic[k].type != VarPAR0; ++k)
  {
    tab[k].index = TypeTupleTabBase<Var>::tabStatic[k].index;
    tab[k].type = TypeTupleTabBase<Var>::tabStatic[k].type;
    tab[k].code = TypeTupleTabBase<Var>::tabStatic[k].code;
    tab[k].name = TypeTupleTabBase<Var>::tabStatic[k].name;
  }
  for (size_t j = 0; j < parent->getNPar(); ++j)
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
    tab[k + parent->getNPar() - 1].type = TypeTupleTabBase<Var>::tabStatic[k].type;
    tab[k + parent->getNPar() - 1].code = TypeTupleTabBase<Var>::tabStatic[k].code;
    tab[k + parent->getNPar() - 1].name = TypeTupleTabBase<Var>::tabStatic[k].name;
  }
}
