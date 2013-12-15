// ------------------------------------------------------------------------- //
//
// This is part of KNUT
// Copyright (c) 2013 by Robert Szalai
//
// For license, see the file COPYING in the root directory of the package
//
// ------------------------------------------------------------------------- //

#include "mat4data.h"
#include "basecomp.h"
#include "mex.h"

class KNConstants;

class fieldWriter
{
public:
  void addField(const char *field, const mxArray*  src)
  {
    if (src)
    {
      fieldname.push_back(field);
      source.push_back(const_cast<mxArray*>(src));
    }
  }
  mxArray* write()
  {
    if (fieldname.size() > 0)
    {
      const char * fields[fieldname.size()];
      for (size_t k = 0; k < fieldname.size(); k++) fields[k] = fieldname[k];
      mxArray* res = mxCreateStructMatrix(1, 1, fieldname.size(), fields);
      for (size_t k = 0; k < fieldname.size(); k++) mxSetField(res, 0, fieldname[k], source[k]);
      return res;
    }
    return 0;
  }
private:
  std::vector<const char *> fieldname;
  std::vector<mxArray*> source;
};

class KNMatlabData : public KNAbstractData
{
  public:
    KNMatlabData(const std::vector<std::string>& parNames, size_t steps_, size_t ndim_, size_t npar_, size_t nint_, size_t ndeg_, size_t nmul_);
    // Constructor for quasi-periodic orbits
    KNMatlabData(const std::vector<std::string>& parNames, size_t steps_, size_t ndim_, size_t npar_, size_t nint1_, size_t nint2_, size_t ndeg1_, size_t ndeg2_);
    // Constructor for opening an existing file
    KNMatlabData(const mxArray* data);
    // unmaps the memory, truncates the file if necessary, closes the file
    virtual ~KNMatlabData();
    void lockRead() const;
    void lockWrite() const;
    void unlock() const;
    void   setNTrivMul(const size_t j, size_t i);
    void   setMagic(size_t n, int32_t magic);
    void   setPar(size_t n, const KNVector& par);
    void   setMul(size_t n, const KNVector& real, const KNVector& imag);
    void   setElem(size_t n, const KNVector& el);
    void   setMesh(size_t n, const KNVector& mesh);
    void   setProfile(size_t n, const KNVector& profile);
    void   getPar(size_t n, KNVector& par) const;
    void   getMul(size_t n, KNVector& real, KNVector& imag) const;
    void   getMesh(size_t n, KNVector& mesh) const;
    size_t  getNPar() const;
    size_t  getNDim() const;
    size_t  getNInt() const;
    size_t  getNDeg() const;
    size_t  getNCols() const;
    int32_t getMagic(size_t n) const;
    void   getProfile(size_t n, KNVector& profile) const;
    // the extra bit for tori
    void   setMesh1(size_t n, size_t j, double d);
    void   setMesh2(size_t n, size_t j, double d);
    void   getBlanket(size_t n, KNVector& blanket) const;
    void   setBlanket(size_t n, const KNVector& blanket);
    
    mxArray* getData() { return mx_alldata; }
    
    static KNMatlabData* createDataStatic (DataType t, size_t ndim, size_t npar, KNConstants* prms);

  private:
    size_t size;
    size_t ncols;
    // TODO implement write permission
//     const bool wperm;
    size_t ndim;
    size_t npar;
    size_t nint;
    size_t ndeg;
    size_t nmul;

    size_t nint1;
    size_t nint2;
    size_t ndeg1;
    size_t ndeg2;
    
    bool torus;

    const mxArray* mx_npoints;
    const mxArray* mx_par;
    const mxArray* mx_parnames;
    const mxArray* mx_mul;
    const mxArray* mx_ntrivmul;
    const mxArray* mx_magic;
    const mxArray* mx_ndim;     //T
    const mxArray* mx_elem;
    const mxArray* mx_mesh;
    const mxArray* mx_prof;
    // for the torus
    const mxArray* mx_nint1;     //T
    const mxArray* mx_nint2;     //T
    const mxArray* mx_ndeg1;     //T
    const mxArray* mx_ndeg2;     //T
    const mxArray* mx_mesh1;     //T
    const mxArray* mx_mesh2;     //T
    const mxArray* mx_blanket;   //T
    mxArray* mx_alldata;
};
