// ------------------------------------------------------------------------- //
//
// This is part of KNUT
// Copyright (c) 2006 by Robert Szalai
//
// For license, see the file COPYING in the root directory of the package
//
// ------------------------------------------------------------------------- //

#include <string>
#include <cmath>
#include <cfloat>
#include <cstring>
#include <fstream>

#ifndef _WIN32
#  include <sys/types.h>
#  include <sys/stat.h>
#  include <fcntl.h>
#  include <sys/mman.h>
#  include <cerrno>
#  include <unistd.h>
#  include <sys/file.h>
#endif

#include "matrix.h"
#include "knerror.h"
#include "mat4data.h"
#include "multipliers.h"

#define MAT_BIG_ENDIAN 1000
#define MAT_LITTLE_ENDIAN 0
#define MAT_NUMERIC 0
#define MAT_TEXT 1

static int32_t byte_order()
{
  int32_t tmp = 0x01020304;
  if ( ((char*)(&tmp))[0] == 1 ) return MAT_BIG_ENDIAN;
  else if ( ((char*)(&tmp))[0] == 4 ) return MAT_LITTLE_ENDIAN;
  else P_MESSAGE1( "Fatal error. Unrecognized byte order." );
}

// this is one big global mutex.
std::mutex KNDataFile::fileLock;

void KNDataFile::lockRead() const
{
  fileLock.lock();
}

void KNDataFile::lockWrite() const
{
  fileLock.lock();
}

void KNDataFile::unlock() const
{
  fileLock.unlock();
}

size_t KNDataFile::findMatrix(const char* name, KNDataFile::header* found, bool test, size_t r, size_t c, uint32_t imag, const char* fileName, int32_t type)
{
  struct header hd;
  size_t cur_off = 0;
  size_t cur_size;
  do
  {
    memcpy(&hd, (char*)address + cur_off, sizeof(struct header));
    if (hd.imagf == 0)
      cur_size = sizeof(struct header) + hd.namelen * sizeof(char) + hd.mrows * hd.ncols * sizeof(double);
    else
      cur_size = sizeof(struct header) + hd.namelen * sizeof(char) + 2 * hd.mrows * hd.ncols * sizeof(double);
    if (strncmp(name, (char*)address + cur_off + sizeof(struct header), 20) == 0)
    {
      memcpy(found, &hd, sizeof(struct header));
      // Checking the parameters
      if (r != 0) P_ERROR_X5(hd.mrows == r, "Wrong number of rows of '", name, "' in file '", fileName, "'.");
      if (c != 0) P_ERROR_X5(hd.ncols == c, "Wrong number of columns of '", name, "' in file '", fileName, "'.");
      if ((imag == 0)&&(hd.imagf != 0)) P_MESSAGE5("KNMatrix '", name, "' in file '", fileName, "' has complex elements, but real was expected.");
      if ((imag != 0)&&(hd.imagf == 0)) P_MESSAGE5("KNMatrix '", name, "' in file '", fileName, "' has real elements, but complex was expected.");
      // check type
      P_ERROR_X5(hd.type == (byte_order()+type), "'", name, " in file '", fileName, "' is not a matrix of required type.");
      return cur_off;
    }
    cur_off += cur_size;
  }
  while (cur_off + sizeof(struct header) < size);
  if (test) P_MESSAGE5("Could not find '", name, "' in file '", fileName, "'.");
  return size;
}

void KNDataFile::mmapFileWrite(const std::string& fileName, size_t size)
{
  // Remove the existing file
  file_mapping::remove(fileName.c_str());
  std::filebuf fbuf;
  fbuf.open(fileName, std::ios_base::in | std::ios_base::out
                      | std::ios_base::trunc | std::ios_base::binary);
  // Set the size
  fbuf.pubseekoff(size-1, std::ios_base::beg);
  fbuf.sputc(0);
  fbuf.close();
//   std::cout << " made a file " << fileName << " of size " << size << "\n";
  try {
    // Create a file mapping
    file_mapping m_file(fileName.c_str(), read_write);
    matFileMapping.swap(m_file);
  
    //Map the whole file with read-write permissions in this process
    mapped_region region(matFileMapping, read_write);
    matFileMemory.swap(region);
  
    //Get the address of the mapped region
    address = matFileMemory.get_address();
  //  std::size_t size  = matFileMemory.get_size();
  }
  catch (interprocess_exception& exc)
  {
    unlock();
    P_MESSAGE4("Boost::interprocess_exception. MAT=", fileName, " Problem: ", exc.what());
  }
}

void KNDataFile::mmapFileRead(const std::string& fileName)
{
  try {
    // Create a file mapping
    file_mapping m_file(fileName.c_str(), read_only);
    matFileMapping.swap(m_file);
  
    //Map the whole file with read-write permissions in this process
    mapped_region region(matFileMapping, read_only);
    matFileMemory.swap(region);
  
    //Get the address of the mapped region
    address = matFileMemory.get_address();
    size = matFileMemory.get_size();
//     std::cout << " read a file " << fileName << " of size " << size << "\n";
  }
  catch (interprocess_exception& exc)
  {
    P_MESSAGE4("Boost::interprocess_exception. MAT=", fileName, " Problem: ", exc.what());
  }  
}

// returns the size of the whole data
inline size_t KNDataFile::createMatrixHeader(KNDataFile::header* hd, const char* name, size_t rows, size_t cols, int32_t type)
{
  hd->type = byte_order() + type;
  hd->mrows = static_cast<uint32_t>(rows);
  hd->ncols = static_cast<uint32_t>(cols);
  hd->imagf = 0;
  hd->namelen = static_cast<uint32_t>(((strlen(name) + sizeof(KNDataFile::header) + 1) / sizeof(double) + 1) * sizeof(double) - sizeof(KNDataFile::header));
  return sizeof(KNDataFile::header) + hd->namelen + rows * cols * sizeof(double);
}

// returns the size of the whole data
inline size_t KNDataFile::createComplexMatrixHeader(KNDataFile::header* hd, const char* name, size_t rows, size_t cols, int32_t type)
{
  hd->type = byte_order() + type;
  hd->mrows = static_cast<uint32_t>(rows);
  hd->ncols = static_cast<uint32_t>(cols);
  hd->imagf = 1;
  hd->namelen = static_cast<uint32_t>(((strlen(name) + sizeof(KNDataFile::header) + 1) / sizeof(double) + 1) * sizeof(double) - sizeof(KNDataFile::header));
  return sizeof(KNDataFile::header) + hd->namelen + 2 * rows * cols * sizeof(double);
}

inline void KNDataFile::writeMatrixHeader(void* address, size_t offset, KNDataFile::header* hd, const char* name)
{
  strncpy((char*)address + offset + sizeof(KNDataFile::header), name, hd->namelen);
  *((KNDataFile::header*)((char*)address + offset)) = *hd;
}

KNDataFile::KNDataFile(const std::string& fileName, const std::vector<std::string>& parNames, size_t steps_, size_t ndim_, size_t npar_, size_t nint_, size_t ndeg_, size_t nmul_)
  : matFileName(fileName), wperm(true)
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

  // creating the matrices
  char npoints_string[] = "knut_npoints";
  npoints_offset = 0;
  size_t npoints_size = createMatrixHeader(&npoints_header, npoints_string, 1, 1);

  char par_string[] = "knut_par";
  par_offset = npoints_offset + npoints_size;
  size_t par_size = createMatrixHeader(&par_header, par_string, npar, ncols);
  
  char parnames_string[] = "knut_parnames";
  parnames_offset = par_offset + par_size;
  size_t parnames_size = createMatrixHeader(&parnames_header, parnames_string,
                           parNames.size(), max_namesize, MAT_TEXT);

  char mul_string[] = "knut_mul";
  mul_offset = parnames_offset + parnames_size;
  size_t mul_size = createComplexMatrixHeader(&mul_header, mul_string, nmul, ncols);

  char ntrivmul_string[] = "knut_ntrivmul";
  ntrivmul_offset = mul_offset + mul_size;
  size_t ntrivmul_size = createMatrixHeader(&ntrivmul_header, ntrivmul_string, 3, 1);

  char magic_string[] = "knut_magic";
  magic_offset = ntrivmul_offset + ntrivmul_size;
  size_t magic_size = createMatrixHeader(&magic_header, magic_string, 1, ncols);
  
  char ndim_string[] = "knut_ndim";
  ndim_offset = magic_offset + magic_size;
  size_t ndim_size = createMatrixHeader(&ndim_header, ndim_string, 1, ncols);

  char elem_string[] = "knut_elem";
  elem_offset = ndim_offset + ndim_size;
  size_t elem_size = createMatrixHeader(&elem_header, elem_string, ndeg + 1, ncols);

  char mesh_string[] = "knut_mesh";
  mesh_offset = elem_offset + elem_size;
  size_t mesh_size = createMatrixHeader(&mesh_header, mesh_string, nint + 1, ncols);

  char prof_string[] = "knut_prof";
  prof_offset = mesh_offset + mesh_size;
  size_t prof_size = createMatrixHeader(&prof_header, prof_string, ndim * (ndeg * nint + 1), ncols);
  size = prof_offset + prof_size;

//  const size_t approxSize = 8 * (sizeof(header) + 20) + sizeof(double) * (1 + ncols * (npar + 2 * nmul + 1 + (ndeg + 1) + (ndim + 1) * (ndeg * nint + 1)));
  lockWrite();
  mmapFileWrite(fileName, size);

  writeMatrixHeader(address, npoints_offset, &npoints_header, npoints_string);
  writeMatrixHeader(address, par_offset, &par_header, par_string);
  writeMatrixHeader(address, parnames_offset, &parnames_header, parnames_string);
  writeMatrixHeader(address, mul_offset, &mul_header, mul_string);
  writeMatrixHeader(address, ntrivmul_offset, &ntrivmul_header, ntrivmul_string);
  writeMatrixHeader(address, magic_offset, &magic_header, magic_string);
  writeMatrixHeader(address, ndim_offset, &ndim_header, ndim_string);
  writeMatrixHeader(address, elem_offset, &elem_header, elem_string);
  writeMatrixHeader(address, mesh_offset, &mesh_header, mesh_string);
  writeMatrixHeader(address, prof_offset, &prof_header, prof_string);
  setParNames(parNames);
//   std::cout << "mat4data address WR1 " << address << "\n";
  unlock();  
}

void KNDataFile::setParNames(const std::vector<std::string>& parNames)
{
  for (size_t i=0; i<parNames.size(); ++i)
  {
    for (size_t j=0; j<parNames[i].size(); ++j)
    {
      elem(parnames_offset, i, j) = parNames[i][j];
    }
    for (size_t j=parNames[i].size(); j<getHeader(parnames_offset)->getCols(); ++j)
    {
      elem(parnames_offset, i, j) = 0;
    }
  }
}

void KNDataFile::getParNames(std::vector<std::string>& parNames) const
{
  parNames.resize(getRows(parnames_offset));
  char* buf = new char[getCols(parnames_offset)+1];
  for (size_t i = 0; i < getRows(parnames_offset); ++i)
  {
    for (size_t j = 0; j < getCols(parnames_offset); ++j)
    {
      buf[j] = (char)elem(parnames_offset, i, j);
    }
    buf[getCols(parnames_offset)] = '\0';
    parNames[i] = buf;
  }
  delete[] buf;
}

KNDataFile::KNDataFile(const std::string& fileName, const std::vector<std::string>& parNames, size_t steps_, size_t ndim_, size_t npar_, size_t nint1_, size_t nint2_, size_t ndeg1_, size_t ndeg2_)
  : matFileName(fileName), wperm(true)
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
  
  // creating the matrices
  char npoints_string[] = "knut_npoints";
  npoints_offset = 0;
  size_t npoints_size = createMatrixHeader(&npoints_header, npoints_string, 1, 1);

  char par_string[] = "knut_par";
  par_offset = npoints_offset + npoints_size;
  size_t par_size = createMatrixHeader(&par_header, par_string, npar, ncols);
  
  char parnames_string[] = "knut_parnames";
  parnames_offset = par_offset + par_size;
  size_t parnames_size = createMatrixHeader(&parnames_header, parnames_string,
                           parNames.size(), max_namesize, MAT_TEXT);

  char ndim_string[] = "knut_ndim";
  ndim_offset = parnames_offset + parnames_size;
  size_t ndim_size = createMatrixHeader(&ndim_header, ndim_string, 1, ncols);

  char nint1_string[] = "knut_nint1";
  nint1_offset = ndim_offset + ndim_size;
  size_t nint1_size = createMatrixHeader(&nint1_header, nint1_string, 1, ncols);

  char nint2_string[] = "knut_nint2";
  nint2_offset = nint1_offset + nint1_size;
  size_t nint2_size = createMatrixHeader(&nint2_header, nint2_string, 1, ncols);

  char ndeg1_string[] = "knut_ndeg1";
  ndeg1_offset = nint2_offset + nint2_size;
  size_t ndeg1_size = createMatrixHeader(&ndeg1_header, ndeg1_string, 1, ncols);

  char ndeg2_string[] = "knut_ndeg2";
  ndeg2_offset = ndeg1_offset + ndeg1_size;
  size_t ndeg2_size = createMatrixHeader(&ndeg2_header, ndeg2_string, 1, ncols);

  char mesh1_string[] = "knut_mesh1";
  mesh1_offset = ndeg2_offset + ndeg2_size;
  size_t mesh1_size = createMatrixHeader(&mesh1_header, mesh1_string, nint1 * ndeg1, ncols);

  char mesh2_string[] = "knut_mesh2";
  mesh2_offset = mesh1_offset + mesh1_size;
  size_t mesh2_size = createMatrixHeader(&mesh2_header, mesh2_string, nint2 * ndeg2, ncols);

  char blanket_string[] = "knut_blanket";
  blanket_offset = mesh2_offset + mesh2_size;
  size_t blanket_size = createMatrixHeader(&blanket_header, blanket_string, ndim * nint1 * ndeg1 * nint2 * ndeg2, ncols);
  size = blanket_offset + blanket_size;

  lockWrite();
  mmapFileWrite(fileName, size);

  writeMatrixHeader(address, npoints_offset, &npoints_header, npoints_string);
  writeMatrixHeader(address, par_offset, &par_header, par_string);
  writeMatrixHeader(address, parnames_offset, &parnames_header, parnames_string);
  writeMatrixHeader(address, ndim_offset, &ndim_header, ndim_string);
  writeMatrixHeader(address, nint1_offset, &nint1_header, nint1_string);
  writeMatrixHeader(address, nint2_offset, &nint2_header, nint2_string);
  writeMatrixHeader(address, ndeg1_offset, &ndeg1_header, ndeg1_string);
  writeMatrixHeader(address, ndeg2_offset, &ndeg2_header, ndeg2_string);
  writeMatrixHeader(address, mesh1_offset, &mesh1_header, mesh1_string);
  writeMatrixHeader(address, mesh2_offset, &mesh2_header, mesh2_string);
  writeMatrixHeader(address, blanket_offset, &blanket_header, blanket_string);
  setParNames(parNames);
//   std::cout << "mat4data address WR2 " << address << "\n";
  unlock();
}

KNDataFile::KNDataFile(const std::string& fileName)
  : matFileName(fileName), wperm(false)
{
  openReadOnly(fileName);
}

void KNDataFile::resizeMatrix(const char* name, size_t newcol, uint32_t imag, int32_t type)
{
  header mathead;
  size_t matoffset;
  if ((matoffset = findMatrix(name, &mathead, false, 0, 0, imag, "", type)) >= size) P_MESSAGE3("KNMatrix '", name, "' cannot be found.");
  P_ASSERT_X1(mathead.ncols >= newcol, "Cannot increase the size of the matrix.");
  char* to = (char*)address + matoffset + mathead.enddata(newcol);
  char* from = (char*)address + matoffset + mathead.enddata(mathead.ncols);
  size_t count = size - matoffset - mathead.enddata(mathead.ncols);
  if (mathead.imagf!= 0)
  {
    // copy the imaginary part
    char* to_im = (char*)address + matoffset + mathead.col_off(newcol);
    char* from_im = (char*)address + matoffset + mathead.col_off(mathead.ncols);
    memmove(to_im, from_im, mathead.col_off(newcol));
  }
  memmove(to, from, count);
  size -= mathead.enddata(mathead.ncols) - mathead.enddata(newcol);
  mathead.ncols = static_cast<uint32_t>(newcol);
  writeMatrixHeader(address, matoffset, &mathead, name);
}

void KNDataFile::condenseData()
{
//   std::cout << "Writable " << wperm << "\n";
  // condense the fields so that it is only the size of npoints, where it counts
  if (wperm)
  {
    lockWrite();
    const size_t npoints = getNPoints();
    if (!torus)
    {
      resizeMatrix("knut_prof", npoints, 0);
      resizeMatrix("knut_mesh", npoints, 0);
      resizeMatrix("knut_elem", npoints, 0);
      resizeMatrix("knut_ndim", npoints, 0);
      resizeMatrix("knut_magic", npoints, 0);
      resizeMatrix("knut_mul", npoints, 1);
      resizeMatrix("knut_par", npoints, 0);
    }
    unlock();
  }
}

void KNDataFile::initHeaders()
{
  npoints_offset = findMatrix("knut_npoints", &npoints_header, true, 1, 1, 0, matFileName.c_str());
  
  par_offset = findMatrix("knut_par", &par_header, true, 0, 0, 0, matFileName.c_str());
  npar = par_header.mrows;
  ncols = par_header.ncols;
  
  parnames_offset = findMatrix("knut_parnames", &parnames_header, true, 0, 0, 0, matFileName.c_str(), MAT_TEXT);

  ndim_offset = findMatrix("knut_ndim", &ndim_header, true, 1, ncols, 0, matFileName.c_str());
  ndim = static_cast<size_t>(*((double*)((char*)address + ndim_offset + ndim_header.col_off(0))));

  if ((mul_offset = findMatrix("knut_mul", &mul_header, false, 0, ncols, 1, matFileName.c_str())) < size)
  {
    torus = false;
    nmul = mul_header.mrows;
    
    // periodic solutions
    magic_offset = findMatrix("knut_magic", &magic_header, true, 1, ncols, 0, matFileName.c_str());
    
    ntrivmul_offset = findMatrix("knut_ntrivmul", &ntrivmul_header, true, 3, 1, 0, matFileName.c_str());

    elem_offset = findMatrix("knut_elem", &elem_header, true, 0, ncols, 0, matFileName.c_str());
    ndeg = elem_header.mrows - 1;

    mesh_offset = findMatrix("knut_mesh", &mesh_header, true, 0, ncols, 0, matFileName.c_str());
    nint = mesh_header.mrows - 1;

    prof_offset = findMatrix("knut_prof", &prof_header, true, ndim*(ndeg*nint + 1), ncols, 0, matFileName.c_str());
  }
  else
  {
    torus = true;
    // quasiperiodic solutions
    nint1_offset = findMatrix("knut_nint1", &nint1_header, true, 1, ncols, 0, matFileName.c_str());
    nint1 = static_cast<size_t>(*((double*)((char*)address + nint1_offset + nint1_header.col_off(0))));

    nint2_offset = findMatrix("knut_nint2", &nint2_header, true, 1, ncols, 0, matFileName.c_str());
    nint2 = static_cast<size_t>(*((double*)((char*)address + nint2_offset + nint2_header.col_off(0))));

    ndeg1_offset = findMatrix("knut_ndeg1", &ndeg1_header, true, 1, ncols, 0, matFileName.c_str());
    ndeg1 = static_cast<size_t>(*((double*)((char*)address + ndeg1_offset + ndeg1_header.col_off(0))));

    ndeg2_offset = findMatrix("knut_ndeg2", &ndeg2_header, true, 1, ncols, 0, matFileName.c_str());
    ndeg2 = static_cast<size_t>(*((double*)((char*)address + ndeg2_offset + ndeg2_header.col_off(0))));

    mesh1_offset = findMatrix("knut_mesh1", &mesh1_header, true, nint1*ndeg1, ncols, 0, matFileName.c_str());

    mesh2_offset = findMatrix("knut_mesh2", &mesh2_header, true, nint2*ndeg2, ncols, 0, matFileName.c_str());

    blanket_offset = findMatrix("knut_blanket", &blanket_header, true, ndim*nint1*ndeg1*nint2*ndeg2, ncols, 0, matFileName.c_str());
  }
}

void KNDataFile::openReadOnly(const std::string& fileName)
{
  mmapFileRead(fileName);
  lockRead();
  initHeaders();
//   std::cout << "mat4data address RD " << address << "\n";
  unlock();
}

KNDataFile::~KNDataFile()
{
  // at the moment we try to truncate the file
  if( wperm )
  {
    const size_t mapsz = matFileMemory.get_size();
    if ( mapsz < size )
    {
      try
      {
        matFileMemory.shrink_by( mapsz - size );
        // internal Boost 
        ipcdetail::truncate_file ( matFileMapping.get_mapping_handle().handle, size );
      }
      catch (interprocess_exception& exc)
      {
        std::cerr << "Boost::interprocess_exception. MAT=" << matFileName << " Problem: " << exc.what() << " THIS WILL ABORT THE PROGRAM\n";
        abort();
      }  
    }
  }
  address = nullptr;
}

void KNDataFile::setPar(size_t n, const KNVector& par)
{
  if (wperm && n < ncols)
  {
    if (par.size() <= npar)
    {
      for (size_t i = 0; i < par.size(); ++i)
        elem(par_offset, i, n) = par(i);
      for (size_t i = par.size(); i < npar; ++i)
        elem(par_offset, i, n) = 0.0;
    }
    else
    {
      P_MESSAGE1("Error while writing the MAT file.");
    }
  }
  else
  {
    P_MESSAGE1("Error while writing the MAT file.");
  }
}

void KNDataFile::setMul(size_t n, const KNVector& re, const KNVector& im)
{
  if (wperm && n < ncols && re.size() == im.size())
  {
    if (re.size() <= nmul)
    {
      for (size_t i = 0; i < re.size(); ++i)
      {
        elem(mul_offset, i, n) = re(i);
        elem_im(mul_offset, i, n) = im(i);
      }
      for (size_t i = re.size(); i < nmul; ++i)
      {
        elem(mul_offset, i, n) = 0.0;
        elem_im(mul_offset, i, n) = 0.0;
      }
    }
    else
    {
      P_MESSAGE1("Error while writing the MAT file.");
    }
  }
  else
  {
    P_MESSAGE1("Error while writing the MAT file.");
  }
}

void KNDataFile::setElem(size_t n, const KNVector& el)
{
  if (wperm && n < ncols)
  {
    if (el.size() == ndeg + 1)
    {
      for (size_t i = 0; i < ndeg + 1; ++i)
        elem(elem_offset, i, n) = el(i);
    }
    else
    {
      P_MESSAGE1("Error while writing the MAT file.");
    }
  }
  else
  {
    P_MESSAGE1("Error while writing the MAT file.");
  }
}

void KNDataFile::setMesh(size_t n, const KNVector& mesh)
{
  if (wperm && n < ncols)
  {
    if (mesh.size() == nint + 1)
    {
      for (size_t i = 0; i < nint + 1; ++i)
        elem(mesh_offset, i, n) = mesh(i);
    }
    else
    {
      P_MESSAGE1("Error while writing the MAT file.");
    }
  }
  else
  {
    P_MESSAGE1("Error while writing the MAT file.");
  }
}

void KNDataFile::setProfile(size_t n, const KNVector& prof)
{
  if (wperm && n < ncols)
  {
    if (prof.size() == ndim*(ndeg*nint + 1))
    {
      const size_t curr_npoints = static_cast<size_t>(elem(npoints_offset, 0, 0));
      if (n+1 > curr_npoints) elem(npoints_offset, 0, 0) = n+1;
      elem(ndim_offset, 0, n) = ndim;
      for (size_t i = 0; i < ndim*(ndeg*nint + 1); ++i) elem(prof_offset, i, n) = prof(i);
    }
    else
    {
      P_MESSAGE4("Error while writing the MAT file: (prof:size vs. ndim*(ndeg*nint + 1))", prof.size(), " vs. ", ndim*(ndeg*nint + 1));
    }
  }
  else
  {
    P_MESSAGE1("Error while writing the MAT file.");
  }
}

void KNDataFile::getBlanket(size_t n, KNVector& blanket) const
{
  const size_t curr_npoints = static_cast<size_t>(elem(npoints_offset, 0, 0));
  if ((blanket.size() == ndim*(ndeg1*nint1*ndeg2*nint2))&&(n < curr_npoints))
  {
    for (size_t i = 0; i < ndim*(ndeg1*nint1*ndeg2*nint2); ++i) blanket(i) = elem(blanket_offset, i, n);
  }
  else
  {
    P_MESSAGE1("Error while reading the MAT file.");
  }
}

// write in order size() = [ ndeg1*nint1, ndeg2*nint2, ndim ]
void KNDataFile::setBlanket(size_t n, const KNVector& blanket)
{
  if (wperm && n < ncols)
  {
    if (blanket.size() == ndim*(ndeg1*nint1*ndeg2*nint2))
    {
      const size_t curr_npoints = static_cast<size_t>(elem(npoints_offset, 0, 0));
      if (n+1 > curr_npoints) elem(npoints_offset, 0, 0) = n+1;
      elem(ndim_offset, 0, n) = ndim;
      elem(nint1_offset, 0, n) = nint1;
      elem(nint2_offset, 0, n) = nint2;
      elem(ndeg1_offset, 0, n) = ndeg1;
      elem(ndeg2_offset, 0, n) = ndeg2;
//    for( int i1 = 0; i1 < nint1; ++i1 ) {
//     for( int i2 = 0; i2 < nint2; ++i2 ) {
//      for( int j1 = 0; j1 < ndeg1; ++j1 ) {
//       for( int j2 = 0; j2 < ndeg2; ++j2 ) {
//        for( int p = 0; p < ndim; ++p ) {
//         elem( blanket_offset, j1 + ndeg1*(i1 + nint1*(j2 + ndeg2*(i2 + nint2*p))), n )
//          = blanket(p + ndim*(j1 + ndeg1*(j2 + ndeg2*(i1 + nint1*i2))));
//        }
//       }
//      }
//     }
//    }
      for (size_t i = 0; i < ndim*(ndeg1*nint1*ndeg2*nint2); ++i) elem(blanket_offset, i, n) = blanket(i);
    }
    else
    {
      P_MESSAGE1("Error while writing the MAT file.");
    }
  }
  else
  {
    P_MESSAGE1("Cannot write such point in MAT file.");
  }
}

void KNDataFile::getPar(size_t n, KNVector& par) const
{
  if (n < ncols)
  {
    if (par.size() <= npar)
    {
      for (size_t i = 0; i < par.size(); ++i)
        par(i) = elem(par_offset, i, n);
    }
    else
    {
      P_MESSAGE1("Error while reading the MAT file.");
    }
  }
  else
  {
    P_MESSAGE1("No such point in MAT file.");
  }
}

void KNDataFile::getMul(size_t n, KNVector& re, KNVector& im) const
{
  if (n < ncols && re.size() == im.size())
  {
    const size_t sz = std::min<size_t>(re.size(), nmul);
    for (size_t i = 0; i < sz; ++i)
    {
      re(i) = elem(mul_offset, i, n);
      im(i) = elem_im(mul_offset, i, n);
    }
  }
  else
  {
    P_MESSAGE1("No such point in MAT file.");
  }
}

void KNDataFile::getElem(size_t n, KNVector& el) const
{
  if (n < ncols)
  {
    if (el.size() == ndeg + 1)
    {
      for (size_t i = 0; i < ndeg + 1; ++i)
        el(i) = elem(elem_offset, i, n);
    }
    else
    {
      P_MESSAGE1("Error while reading the MAT file.");
    }
  }
  else
  {
    P_MESSAGE1("No such point in MAT file.");
  }
}

void KNDataFile::getMesh(size_t n, KNVector& mesh) const
{
  if (n < ncols)
  {
    if (mesh.size() == nint + 1)
    {
      for (size_t i = 0; i < nint + 1; ++i)
        mesh(i) = elem(mesh_offset, i, n);
    }
    else
    {
      P_MESSAGE1("Error while reading the MAT file.");
    }
  }
  else
  {
    P_MESSAGE1("No such point in MAT file.");
  }
}

void KNDataFile::getProfile(size_t n, KNVector& prof) const
{
  if (n < ncols)
  {
    if (prof.size() == ndim*(ndeg*nint + 1))
    {
      for (size_t i = 0; i < ndim*(ndeg*nint + 1); ++i)
        prof(i) = elem(prof_offset, i, n);
    }
    else
    {
      P_MESSAGE1("Error while reading the MAT file.");
    }
  }
  else
  {
    P_MESSAGE1("No such point in MAT file.");
  }
}

size_t KNDataFile::getUnstableMultipliers(size_t n) const
{
  KNVector mulRe(false), mulIm(false);
  const_cast<KNDataFile*>(this)->getMulReRef(n, mulRe);
  const_cast<KNDataFile*>(this)->getMulImRef(n, mulIm);
  const size_t lp = getNTrivMul(0);
  const size_t pd = getNTrivMul(1);
  const size_t ns = getNTrivMul(2);
  return unstableMultipliers(mulRe, mulIm, lp, pd, ns, n);
}

size_t KNDataFile::getNextBifurcation(size_t n, bool* stab) const
{
  KNVector mulRe(false), mulIm(false);
  const_cast<KNDataFile*>(this)->getMulReRef(n, mulRe);
  const_cast<KNDataFile*>(this)->getMulImRef(n, mulIm);
  const size_t lp = getNTrivMul(0);
  const size_t pd = getNTrivMul(1);
  const size_t ns = getNTrivMul(2);
  size_t p_ustab = unstableMultipliers(mulRe, mulIm, lp, pd, ns, n);
  for (size_t i = n + 1; i < getNPoints(); ++i)
  {
    const_cast<KNDataFile*>(this)->getMulReRef(i, mulRe);
    const_cast<KNDataFile*>(this)->getMulImRef(i, mulIm);
    size_t ustab = unstableMultipliers(mulRe, mulIm, lp, pd, ns, i);
    if (ustab != p_ustab) 
    {
      if (stab != nullptr)
      {
        if (ustab == 0 || p_ustab == 0) *stab = true;
        else *stab = false;
      }
      return i;
    }
    p_ustab = ustab;
  }
  return getNPoints();
}

BifType  KNDataFile::getBifurcationType(size_t n) const
{
  if (n > 0)
  {
    KNVector mulReA(false), mulImA(false);
    KNVector mulReB(false), mulImB(false);
    const_cast<KNDataFile*>(this)->getMulReRef(n-1, mulReA);
    const_cast<KNDataFile*>(this)->getMulImRef(n-1, mulImA);
    const_cast<KNDataFile*>(this)->getMulReRef(n, mulReB);
    const_cast<KNDataFile*>(this)->getMulImRef(n, mulImB);

    const size_t lp = getNTrivMul(0);
    const size_t pd = getNTrivMul(1);
    const size_t ns = getNTrivMul(2);

    return bifurcationType(mulReA, mulImA, mulReB, mulImB, lp, pd, ns, n-1, n);
  }
  return BifType::BifNone;
}

size_t KNAbstractData::findType(BifType type, size_t n) const
{
  size_t found = 0;
  for (size_t i=0; i<getNCols(); ++i)
  {
    if (getMagic(i) == type) found++;
    if (found == n) return i;
  }
  return getNCols();
}
