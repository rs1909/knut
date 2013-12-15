// ------------------------------------------------------------------------- //
//
// This is part of KNUT
// Copyright (c) 2013 by Robert Szalai
//
// For license, see the file COPYING in the root directory of the package
//
// ------------------------------------------------------------------------- //

#include "mxdata.h"
#include "constants.h"

KNMatlabData::KNMatlabData(const std::vector<std::string>& parNames, size_t steps_, size_t ndim_, size_t npar_, size_t nint_, size_t ndeg_, size_t nmul_)
{
  torus = false;

  ncols = steps_;
  ndim = ndim_;
  npar = npar_;
  nint = nint_;
  ndeg = ndeg_;
  nmul = nmul_;

  size_t max_namesize = 0;
  for (size_t i=0; i<parNames.size(); ++i)
    if (parNames[i].size() > max_namesize) max_namesize = parNames[i].size();

  // create a structure
  mx_npoints = mxCreateDoubleMatrix (1, 1, mxREAL);
  mx_par = mxCreateDoubleMatrix (npar, ncols, mxREAL);
  mx_parnames = mxCreateNumericMatrix (parNames.size(), max_namesize + 1, mxCHAR_CLASS, mxREAL);
  mx_mul = mxCreateDoubleMatrix(nmul, ncols, mxCOMPLEX);
  mx_ntrivmul = mxCreateDoubleMatrix(1, 3, mxREAL);
  mx_magic = mxCreateDoubleMatrix(1, ncols, mxREAL);
  mx_ndim = mxCreateDoubleMatrix(1, 1, mxREAL);
  mx_elem = mxCreateDoubleMatrix(ndeg+1, ncols, mxREAL);
  mx_mesh = mxCreateDoubleMatrix(nint+1, ncols, mxREAL);
  mx_prof = mxCreateDoubleMatrix(ndim * (ndeg * nint + 1), ncols, mxREAL);
  
  fieldWriter writer;
  writer.addField ("npoints", mx_npoints);
  writer.addField ("par", mx_par);
  writer.addField ("parnames", mx_parnames);
  writer.addField ("mul", mx_mul);
  writer.addField ("ntrivmul", mx_ntrivmul);
  writer.addField ("magic", mx_magic);
  writer.addField ("ndim", mx_ndim);
  writer.addField ("elem", mx_elem);
  writer.addField ("mesh", mx_mesh);
  writer.addField ("prof", mx_prof);
  mx_alldata = writer.write();
  
  double* pt = mxGetPr (mx_npoints);
  pt[0] = (double)ncols;
  pt = mxGetPr (mx_ndim);
  pt[0] = (double)ndim;
}

// Constructor for quasi-periodic orbits
KNMatlabData::KNMatlabData(const std::vector<std::string>& parNames, size_t steps_, size_t ndim_, size_t npar_, size_t nint1_, size_t nint2_, size_t ndeg1_, size_t ndeg2_)
{
  torus = true;

  ncols = steps_;
  ndim = ndim_;
  npar = npar_;
  nint1 = nint1_;
  nint2 = nint2_;
  ndeg1 = ndeg1_;
  ndeg2 = ndeg2_;

  size_t max_namesize = 0;
  for (size_t i=0; i<parNames.size(); ++i)
    if (parNames[i].size() > max_namesize) max_namesize = parNames[i].size();

  mx_npoints = mxCreateDoubleMatrix (1, 1, mxREAL);
  mx_par = mxCreateDoubleMatrix (npar, ncols, mxREAL);
  mx_parnames = mxCreateNumericMatrix (parNames.size(), max_namesize + 1, mxINT8_CLASS, mxREAL);
  mx_ndim = mxCreateDoubleMatrix(1, 1, mxREAL);
  mx_nint1 = mxCreateDoubleMatrix(1, 1, mxREAL);
  mx_nint2 = mxCreateDoubleMatrix(1, 1, mxREAL);
  mx_ndeg1 = mxCreateDoubleMatrix(1, 1, mxREAL);
  mx_ndeg2 = mxCreateDoubleMatrix(1, 1, mxREAL);
  mx_mesh1 = mxCreateDoubleMatrix(nint1 * ndeg1, ncols, mxREAL);
  mx_mesh2 = mxCreateDoubleMatrix(nint2 * ndeg2, ncols, mxREAL);
  mx_blanket = mxCreateDoubleMatrix(ndim * nint1 * ndeg1 * nint2 * ndeg2, ncols, mxREAL);

  fieldWriter writer;
  writer.addField ("npoints", mx_npoints);
  writer.addField ("par", mx_par);
  writer.addField ("parnames", mx_parnames);
  writer.addField ("ndim", mx_ndim);
  writer.addField ("nint1", mx_nint1);
  writer.addField ("nint2", mx_nint2);
  writer.addField ("ndeg1", mx_ndeg1);
  writer.addField ("ndeg2", mx_ndeg2);
  writer.addField ("mesh1", mx_mesh1);
  writer.addField ("mesh2", mx_mesh2);
  writer.addField ("blanket", mx_blanket);
  mx_alldata = writer.write();
  
  double* pt = mxGetPr (mx_npoints);
  pt[0] = (double)ncols;
  pt = mxGetPr (mx_nint1);
  pt[0] = (double)nint1;
  pt = mxGetPr (mx_nint2);
  pt[0] = (double)nint2;
  pt = mxGetPr (mx_ndeg1);
  pt[0] = (double)ndeg1;
  pt = mxGetPr (mx_ndeg2);
  pt[0] = (double)ndeg2;

}

static inline const mxArray * mxStructField (const mxArray * data, const char * field)
{
  mxArray * elem = mxGetField (data, 0, field);
  P_ERROR_X3(elem != nullptr, "Cannot find element '", field, "'." );
  return elem;
}

// Constructor for opening an existing file
KNMatlabData::KNMatlabData (const mxArray* data)
{
  if (mxGetField (data, 0, "blanket") != nullptr) torus = true;
  else torus = false;
  
  mx_npoints = mxStructField (data, "npoints");
  mx_par = mxStructField (data, "par");
  mx_parnames = mxStructField (data, "parnames");
  mx_ndim = mxStructField (data, "ndim");

  ncols = mxGetScalar (mx_npoints);
  ndim = mxGetScalar (mx_ndim);
  npar = mxGetM (mx_par);

  if (torus)
  {
    mx_nint1 = mxStructField (data, "nint1");
    mx_nint2 = mxStructField (data, "nint2");
    mx_ndeg1 = mxStructField (data, "ndeg1");
    mx_ndeg2 = mxStructField (data, "ndeg2");
    mx_mesh1 = mxStructField (data, "mesh1");
    mx_mesh2 = mxStructField (data, "mesh2");
    mx_blanket = mxStructField (data, "blanket");
    
    ndeg1 = mxGetScalar (mx_ndeg1);
    ndeg2 = mxGetScalar (mx_ndeg2);
    nint1 = mxGetScalar (mx_nint1);
    nint2 = mxGetScalar (mx_nint2);
  } else
  {
    mx_mul = mxStructField (data, "mul");
    mx_ntrivmul = mxStructField (data, "ntrivmul");
    mx_magic = mxStructField (data, "magic");
    mx_elem = mxStructField (data, "elem");
    mx_mesh = mxStructField (data, "mesh");
    mx_prof = mxStructField (data, "prof");
    
    ndeg = mxGetM (mx_elem) - 1;
    nint = mxGetM (mx_mesh) - 1;
    nmul = mxGetM (mx_mul);
  }
}

// unmaps the memory, truncates the file if necessary, closes the file
KNMatlabData::~KNMatlabData()
{
}

void KNMatlabData::lockRead() const { }
void KNMatlabData::lockWrite() const { }
void KNMatlabData::unlock() const { }

void KNMatlabData::setNTrivMul(const size_t j, size_t i)
{
  double* pt = mxGetPr (mx_ntrivmul);
  pt[j] = (double)i;
}

void KNMatlabData::setMagic(size_t n, int32_t magic)
{
  double* pt = mxGetPr (mx_magic);
  pt[n] = magic;
}

void KNMatlabData::setPar(size_t n, const KNVector& par)
{
  double* pt = mxGetPr (mx_par);
  for (size_t k=0; k < npar; k++) pt[n*npar + k] = par(k);
}

void KNMatlabData::setMul(size_t n, const KNVector& real, const KNVector& imag)
{
  double* pt = mxGetPr (mx_mul);
  for (size_t k=0; k < nmul; k++) pt[n*nmul + k] = real(k);

  double* pi = mxGetPi (mx_mul);
  for (size_t k=0; k < nmul; k++) pi[n*nmul + k] = imag(k);
}

void KNMatlabData::setElem(size_t n, const KNVector& el)
{
  double* pt = mxGetPr (mx_elem);
  for (size_t k=0; k < ndeg+1; k++) pt[n*(ndeg+1) + k] = el(k);
}

void KNMatlabData::setMesh(size_t n, const KNVector& mesh)
{
  double* pt = mxGetPr (mx_mesh);
  for (size_t k=0; k < nint+1; k++) pt[n*(nint+1) + k] = mesh(k);
}

void KNMatlabData::setProfile(size_t n, const KNVector& profile)
{
  double* pt = mxGetPr (mx_prof);
  for (size_t k=0; k < ndim * (ndeg * nint + 1); k++) 
    pt[n*ndim * (ndeg * nint + 1) + k] = profile(k);
}

// the reading functions
void KNMatlabData::getPar(size_t n, KNVector& par) const
{
  const double* pt = mxGetPr (mx_par);
  for (size_t k=0; k < npar; k++) par(k) = pt[n*npar + k];
}

void KNMatlabData::getMul(size_t n, KNVector& real, KNVector& imag) const
{
  const double* pt = mxGetPr (mx_mul);
  for (size_t k=0; k < nmul; k++) real(k) = pt[n*nmul + k];

  const double* pi = mxGetPi (mx_mul);
  for (size_t k=0; k < nmul; k++) imag(k) = pi[n*nmul + k];
}

void KNMatlabData::getMesh(size_t n, KNVector& mesh) const
{
  const double* pt = mxGetPr (mx_mesh);
  for (size_t k=0; k < nint+1; k++) mesh(k) = pt[n*(nint+1) + k];
}

size_t KNMatlabData::getNPar() const { return npar; }
size_t KNMatlabData::getNDim() const { return ndim; }
size_t KNMatlabData::getNInt() const { return nint; }
size_t KNMatlabData::getNDeg() const { return ndeg; }
size_t KNMatlabData::getNCols() const { return ncols; }
int32_t KNMatlabData::getMagic(size_t n) const
{
  const double* pt = mxGetPr (mx_magic);
  return pt[n];
}

void KNMatlabData::getProfile(size_t n, KNVector& profile) const
{
  const double* pt = mxGetPr (mx_prof);
  for (size_t k=0; k < ndim * (ndeg * nint + 1); k++) 
    profile(k) = pt[n*ndim * (ndeg * nint + 1) + k];
}

// the extra bit for tori
void KNMatlabData::setMesh1(size_t n, size_t j, double d)
{
  double* pt = mxGetPr (mx_mesh1);
  pt[n*(nint+1) + j] = d;
}

void KNMatlabData::setMesh2(size_t n, size_t j, double d)
{
  double* pt = mxGetPr (mx_mesh2);
  pt[n*(nint+1) + j] = d;
}

void KNMatlabData::getBlanket(size_t n, KNVector& blanket) const
{
  const double* pt = mxGetPr (mx_blanket);
  for (size_t k=0; k < ndim*(ndeg1*nint1*ndeg2*nint2); k++) 
    blanket(k) = pt[n*ndim*(ndeg1*nint1*ndeg2*nint2) + k];
}

void KNMatlabData::setBlanket(size_t n, const KNVector& blanket)
{
  double* pt = mxGetPr (mx_blanket);
  for (size_t k=0; k < ndim*(ndeg1*nint1*ndeg2*nint2); k++) 
    pt[n*ndim*(ndeg1*nint1*ndeg2*nint2) + k] = blanket(k);
}

KNMatlabData* KNMatlabData::createDataStatic (DataType t, size_t ndim, size_t npar, KNConstants* prms)
{
  if (t == DataType::ST) {
    return new KNMatlabData(prms->getParNames(), prms->getSteps(), 
                            ndim, VarToIndex(VarEnd,npar),
                            0, 0, prms->getNMul());
  } else if (t == DataType::LC) {
    return new KNMatlabData(prms->getParNames(), prms->getSteps(), 
                            ndim, VarToIndex(VarEnd,npar),
                            prms->getNInt(), prms->getNDeg(), prms->getNMul());
  } else if (t == DataType::TR) {
    return new KNMatlabData(prms->getParNames(), prms->getSteps(), 
                            ndim, VarToIndex(VarEnd,npar), prms->getNInt1(), prms->getNInt2(), 
                            prms->getNDeg1(), prms->getNDeg2());
  }
  return nullptr;
}

