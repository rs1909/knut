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

#ifndef WIN32
#  include <sys/types.h>
#  include <sys/stat.h>
#  include <fcntl.h>
#  include <sys/mman.h>
#  include <errno.h>
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

#ifndef WIN32

void KNDataFile::lock() const
{
//  std::cout<<"L " << matFileName << "\n";
  if (flock(file, LOCK_EX) != 0) std::cerr<<"Error locking\n";
}

void KNDataFile::unlock() const
{
//  std::cout<<"U " << matFileName << "\n";
  if (flock(file, LOCK_UN) != 0) std::cerr<<"Error locking\n";
}

#else

void KNDataFile::lock() const
{
  // it is opened synchronous, so LockFileEx will wait until the lock is acquired
//   std::cout << "Trying to Lock " << matFileName <<"\n"; std::cout.flush();
  const_cast<OVERLAPPED*>(&fileOverlapped)->Offset = 0;
  const_cast<OVERLAPPED*>(&fileOverlapped)->OffsetHigh = 0;
  const_cast<OVERLAPPED*>(&fileOverlapped)->hEvent = 0;
  BOOL res = LockFileEx(file, LOCKFILE_EXCLUSIVE_LOCK, 0, filesize & 0xffffffff, filesize >> 32, const_cast<OVERLAPPED*>(&fileOverlapped));
  if( res != TRUE ) 
    P_MESSAGE3("MAT file '", matFileName, "' cannot be locked.");
}

void KNDataFile::unlock() const
{
//   std::cout << "UNLock " << matFileName <<"\n"; std::cout.flush();
  if( UnlockFileEx(file, 0,  filesize & 0xffffffff, filesize >> 32, const_cast<OVERLAPPED*>(&fileOverlapped)) != TRUE ) 
    P_MESSAGE3("MAT file '", matFileName, "' cannot be unlocked.");
}

#endif
off_t KNDataFile::findMatrix(const char* name, KNDataFile::header* found, bool test, int32_t r, int32_t c, int32_t imag, const char* fileName, int32_t type)
{
  struct header hd;
  off_t cur_off = 0;
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
      if (r != -1) P_ERROR_X5(hd.mrows == r, "Wrong number of rows of '", name, "' in file '", fileName, "'.");
      if (c != -1) P_ERROR_X5(hd.ncols == c, "Wrong number of columns of '", name, "' in file '", fileName, "'.");
      if ((imag != -1)&&(imag == 0)) if (hd.imagf != 0) P_MESSAGE5("KNMatrix '", name, "' in file '", fileName, "' has complex elements, but real was expected.");
      if ((imag != -1)&&(imag != 0)) if (hd.imagf == 0) P_MESSAGE5("KNMatrix '", name, "' in file '", fileName, "' has real elements, but complex was expected.");
      // check type
      P_ERROR_X5(hd.type == (byte_order()+type), "'", name, " in file '", fileName, "' is not a matrix of required type.");
      return cur_off;
    }
    cur_off += cur_size;
  }
  while ((size_t)cur_off < size);
  if (test) P_MESSAGE5("Could not find '", name, "' in file '", fileName, "'.");
  return -1;
}

#ifndef WIN32

static inline void *mmapFileWrite(int& file, const std::string& fileName, size_t size)
{
  if ((file = open(fileName.c_str(), O_CREAT | O_RDWR | O_TRUNC, S_IRUSR | S_IWUSR)) == -1)
  {
    P_MESSAGE5("Unable to create the MAT file '", fileName, "'. ", (const char *)strerror(errno), ".");
  }

  if (ftruncate(file, size) != 0)
  {
    P_MESSAGE5("Unable to resize the MAT file '", fileName, "'. ", (const char *)strerror(errno), ".");
  }

//  std::cout << "MMAP - WRITE\n";
  void *address;
  if ((address = mmap(0, size, PROT_WRITE | PROT_READ, MAP_SHARED, file, 0)) == MAP_FAILED)
  {
    P_MESSAGE5("Unable to map the MAT file '", fileName, "' to a memory location. ", (const char *)strerror(errno), ".");
  }
  return address;
}

static inline void *mmapFileRead(int& file, const std::string& fileName, size_t& size)
{
  if ((file = open(fileName.c_str(), O_RDONLY)) == -1)
  {
    P_MESSAGE5("Unable to open the MAT file '", fileName, "' for reading. ", (const char *)strerror(errno), ".");
  }

  struct stat filestat;
  if (fstat(file, &filestat) != 0)
  {
    P_MESSAGE5("Unable to stat the MAT file '", fileName, "'. ", (const char *)strerror(errno), ".");
  }
  size = filestat.st_size;

//  std::cout << "MMAP - READ\n";
  void *address;
  if ((address = mmap(0, size, PROT_READ, MAP_SHARED, file, 0)) == MAP_FAILED)
  {
    P_MESSAGE5("Unable to map the MAT file '", fileName, "' to a memory location. ", (const char *)strerror(errno), ".");
  }
  return address;
}

#else

static inline void *mmapFileWrite(HANDLE& file, HANDLE& mapHandle, const std::string& fileName, size_t size)
{
  if ((file = CreateFile(fileName.c_str(),
                         FILE_WRITE_DATA | FILE_READ_DATA,
                         FILE_SHARE_READ,
                         NULL,
                         CREATE_ALWAYS,
                         FILE_ATTRIBUTE_NORMAL,
                         NULL)) == NULL)
  {
    P_MESSAGE5("Unable to create the MAT file '", fileName, "'. Error code ", static_cast<int>(GetLastError()), ".");
  }

  if (SetFilePointer(file, size, NULL, FILE_BEGIN) == 0)
  {
    P_MESSAGE5("Unable to seek in the MAT file '", fileName, "'. Error code ", static_cast<int>(GetLastError()), ".");
  }
  P_ERROR_X5(SetEndOfFile(file), "Unable to truncate the MAT file '", fileName, "'. Error code ", static_cast<int>(GetLastError()), ".");

  if ((mapHandle = CreateFileMapping(file, NULL, PAGE_READWRITE, 0, size, fileName.c_str())) == 0)
  {
    P_MESSAGE5("Unable to map the MAT file '", fileName, "' to a memory location. Error code ", static_cast<int>(GetLastError()), ".");
  }

  void *address = MapViewOfFile(mapHandle, FILE_MAP_WRITE, 0, 0, 0);
  if (address != NULL) return address;
  else P_MESSAGE5("Unable to view the file map of '", fileName, "'. Error code ", static_cast<int>(GetLastError()), ".");
  return 0;
}

static inline void *mmapFileRead(HANDLE& file, HANDLE& mapHandle, const std::string& fileName, size_t& size)
{
  if ((file = CreateFile(fileName.c_str(),
                         FILE_READ_DATA,
                         FILE_SHARE_READ,
                         NULL,
                         OPEN_EXISTING,
                         FILE_ATTRIBUTE_NORMAL,
                         NULL)) == NULL)
  {
    P_MESSAGE5("Unable to create the MAT file '", fileName, "'. Error code ", static_cast<int>(GetLastError()), ".");
  }

  size = GetFileSize(file, NULL);
  if ((mapHandle = CreateFileMapping(file, NULL, PAGE_READONLY, 0, 0, fileName.c_str())) == 0)
  {
    P_MESSAGE5("Unable to map the MAT file '", fileName, "' to a memory location. Error code ", static_cast<int>(GetLastError()), ".");
  }

  void *address = MapViewOfFile(mapHandle, FILE_MAP_READ, 0, 0, 0);
  if (address != NULL) return address;
  else P_MESSAGE5("Unable to view the file map of '", fileName, "'. Error code ", static_cast<int>(GetLastError()), ".");
  return 0;
}

#endif

// returns the size of the whole data
inline size_t KNDataFile::createMatrixHeader(KNDataFile::header* hd, const char* name, int32_t rows, int32_t cols, int32_t type)
{
  hd->type = byte_order() + type;
  hd->mrows = rows;
  hd->ncols = cols;
  hd->imagf = 0;
  hd->namelen = static_cast<int32_t>(((strlen(name) + sizeof(KNDataFile::header) + 1) / sizeof(double) + 1) * sizeof(double) - sizeof(KNDataFile::header));
  return sizeof(KNDataFile::header) + hd->namelen + rows * cols * sizeof(double);
}

// returns the size of the whole data
inline size_t KNDataFile::createComplexMatrixHeader(KNDataFile::header* hd, const char* name, int32_t rows, int32_t cols, int32_t type)
{
  hd->type = byte_order() + type;
  hd->mrows = rows;
  hd->ncols = cols;
  hd->imagf = 1;
  hd->namelen = static_cast<int32_t>(((strlen(name) + sizeof(KNDataFile::header) + 1) / sizeof(double) + 1) * sizeof(double) - sizeof(KNDataFile::header));
  return sizeof(KNDataFile::header) + hd->namelen + 2 * rows * cols * sizeof(double);
}

inline void KNDataFile::writeMatrixHeader(void* address, size_t offset, KNDataFile::header* hd, const char* name)
{
  strncpy((char*)address + offset + sizeof(KNDataFile::header), name, hd->namelen);
  *((KNDataFile::header*)((char*)address + offset)) = *hd;
}

KNDataFile::KNDataFile(const std::string& fileName, const std::vector<std::string>& parNames, int steps_, int ndim_, int npar_, int nint_, int ndeg_, int nmul_)
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
                           static_cast<int32_t>(parNames.size()), static_cast<int32_t>(max_namesize), MAT_TEXT);

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

  const size_t approxSize = 8 * (sizeof(header) + 20) + sizeof(double) * (1 + ncols * (npar + 2 * nmul + 1 + (ndeg + 1) + (ndim + 1) * (ndeg * nint + 1)));
#ifndef WIN32
  address = mmapFileWrite(file, fileName, size);
#else
  address = mmapFileWrite(file, mapHandle, fileName, size);
#endif

  lock();
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
  for (int i = 0; i < getRows(parnames_offset); ++i)
  {
    for (int j = 0; j < getCols(parnames_offset); ++j)
    {
      buf[j] = (char)elem(parnames_offset, i, j);
    }
    buf[getCols(parnames_offset)] = '\0';
    parNames[i] = buf;
  }
  delete[] buf;
}

KNDataFile::KNDataFile(const std::string& fileName, const std::vector<std::string>& parNames, int steps_, int ndim_, int npar_, int nint1_, int nint2_, int ndeg1_, int ndeg2_)
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
                           static_cast<int32_t>(parNames.size()), static_cast<int32_t>(max_namesize), MAT_TEXT);

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

#ifndef WIN32
  address = mmapFileWrite(file, fileName, size);
#else
  address = mmapFileWrite(file, mapHandle, fileName, size);
#endif
  lock();
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

void KNDataFile::resizeMatrix(const char* name, int newcol, int32_t type)
{
  header mathead;
  off_t matoffset;
  if ((matoffset = findMatrix(name, &mathead, false, -1, -1, -1, "", type)) == -1) P_MESSAGE3("KNMatrix '", name, "' cannot be found.");
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
  mathead.ncols = newcol;
  writeMatrixHeader(address, matoffset, &mathead, name);
}

void KNDataFile::condenseData()
{
//   std::cout << "Writable " << wperm << "\n";
  // condense the fields so that it is only the size of npoints, where it counts
  if (wperm)
  {
    lock();
    const int npoints = getNPoints();
    if (!torus)
    {
      resizeMatrix("knut_prof", npoints);
      resizeMatrix("knut_mesh", npoints);
      resizeMatrix("knut_elem", npoints);
      resizeMatrix("knut_ndim", npoints);
      resizeMatrix("knut_magic", npoints);
      resizeMatrix("knut_mul", npoints);
      resizeMatrix("knut_par", npoints);
    }
    unlock();
  }
}

void KNDataFile::initHeaders()
{
  npoints_offset = findMatrix("knut_npoints", &npoints_header, true, 1, 1, 0, matFileName.c_str());
  
  par_offset = findMatrix("knut_par", &par_header, true, -1, -1, 0, matFileName.c_str());
  npar = par_header.mrows;
  ncols = par_header.ncols;
  
  parnames_offset = findMatrix("knut_parnames", &parnames_header, true, -1, -1, 0, matFileName.c_str(), MAT_TEXT);

  ndim_offset = findMatrix("knut_ndim", &ndim_header, true, 1, ncols, 0, matFileName.c_str());
  ndim = static_cast<int>(*((double*)((char*)address + ndim_offset + ndim_header.col_off(0))));

  if ((mul_offset = findMatrix("knut_mul", &mul_header, false, -1, ncols, 1, matFileName.c_str())) != -1)
  {
    torus = false;
    nmul = mul_header.mrows;
    
    // periodic solutions
    magic_offset = findMatrix("knut_magic", &magic_header, true, 1, ncols, 0, matFileName.c_str());
    
    ntrivmul_offset = findMatrix("knut_ntrivmul", &ntrivmul_header, true, 3, 1, 0, matFileName.c_str());

    elem_offset = findMatrix("knut_elem", &elem_header, true, -1, ncols, 0, matFileName.c_str());
    ndeg = elem_header.mrows - 1;

    mesh_offset = findMatrix("knut_mesh", &mesh_header, true, -1, ncols, 0, matFileName.c_str());
    nint = mesh_header.mrows - 1;

    prof_offset = findMatrix("knut_prof", &prof_header, true, ndim*(ndeg*nint + 1), ncols, 0, matFileName.c_str());
  }
  else
  {
    torus = true;
    // quasiperiodic solutions
    nint1_offset = findMatrix("knut_nint1", &nint1_header, true, 1, ncols, 0, matFileName.c_str());
    nint1 = static_cast<int>(*((double*)((char*)address + nint1_offset + nint1_header.col_off(0))));

    nint2_offset = findMatrix("knut_nint2", &nint2_header, true, 1, ncols, 0, matFileName.c_str());
    nint2 = static_cast<int>(*((double*)((char*)address + nint2_offset + nint2_header.col_off(0))));

    ndeg1_offset = findMatrix("knut_ndeg1", &ndeg1_header, true, 1, ncols, 0, matFileName.c_str());
    ndeg1 = static_cast<int>(*((double*)((char*)address + ndeg1_offset + ndeg1_header.col_off(0))));

    ndeg2_offset = findMatrix("knut_ndeg2", &ndeg2_header, true, 1, ncols, 0, matFileName.c_str());
    ndeg2 = static_cast<int>(*((double*)((char*)address + ndeg2_offset + ndeg2_header.col_off(0))));

    mesh1_offset = findMatrix("knut_mesh1", &mesh1_header, true, nint1*ndeg1, ncols, 0, matFileName.c_str());

    mesh2_offset = findMatrix("knut_mesh2", &mesh2_header, true, nint2*ndeg2, ncols, 0, matFileName.c_str());

    blanket_offset = findMatrix("knut_blanket", &blanket_header, true, ndim*nint1*ndeg1*nint2*ndeg2, ncols, 0, matFileName.c_str());
  }
}

void KNDataFile::openReadOnly(const std::string& fileName)
{
#ifndef WIN32
  address = mmapFileRead(file, fileName, size);
#else
  address = mmapFileRead(file, mapHandle, fileName, size);
#endif
  lock();
  initHeaders();
//   std::cout << "mat4data address RD " << address << "\n";
  unlock();
}

// WARNING This creates a large file and then shrinks it
// Use mmremap instead so that only the smaller file is written.
KNDataFile::~KNDataFile()
{
#ifndef WIN32
  size_t oldsize = size;
  condenseData();
//  std::cout << "MUNMAP\n";
  if (munmap(address, oldsize) != 0)
  {
    P_MESSAGE3("Unable to munmap the MAT file. ", strerror(errno), ".");
  }
  // at the moment we try to truncate the file
  if( wperm )
  {
    if( ftruncate( file, size) != 0 )
    {
      P_MESSAGE3("Unable to truncate the MAT file. ", strerror(errno), ".");
    }
  }
  if (close(file) != 0)
  {
    P_MESSAGE3("Unable to close the MAT file. ", strerror(errno), ".");
  }
#else
  if (address != 0)
  {
    UnmapViewOfFile(address);
    CloseHandle(mapHandle);
    // at the moment it won't be truncated
    //   if( wperm )
    //   {
    //    if( SetFilePointer( file, size, NULL, FILE_BEGIN ) == 0 )
    //    { P_ERROR_X2( false, "Unable to seek in the MAT file.", static_cast<int>(GetLastError()) ); }
    //    P_ERROR_X2( SetEndOfFile( file ), "Unable to truncate the MAT file.", static_cast<int>(GetLastError()) );
    //   }
    CloseHandle(file);
  }
#endif
//   std::cout << "mat4data address CLOSE " << address << "\n";
  address = 0;
}

void KNDataFile::setPar(int n, const KNVector& par)
{
  if (wperm && n < ncols)
  {
    if (par.size() <= npar)
    {
      for (int i = 0; i < par.size(); ++i)
        elem(par_offset, i, n) = par(i);
      for (int i = par.size(); i < npar; ++i)
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

void KNDataFile::setMul(int n, const KNVector& re, const KNVector& im)
{
  if (wperm && n < ncols && re.size() == im.size())
  {
    if (re.size() <= nmul)
    {
      for (int i = 0; i < re.size(); ++i)
      {
        elem(mul_offset, i, n) = re(i);
        elem_im(mul_offset, i, n) = im(i);
      }
      for (int i = re.size(); i < nmul; ++i)
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

void KNDataFile::setElem(int n, const KNVector& el)
{
  if (wperm && n < ncols)
  {
    if (el.size() == ndeg + 1)
    {
      for (int i = 0; i < ndeg + 1; ++i)
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

void KNDataFile::setMesh(int n, const KNVector& mesh)
{
  if (wperm && n < ncols)
  {
    if (mesh.size() == nint + 1)
    {
      for (int i = 0; i < nint + 1; ++i)
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

void KNDataFile::setProfile(int n, const KNVector& prof)
{
  if (wperm && n < ncols)
  {
    if (prof.size() == ndim*(ndeg*nint + 1))
    {
      const int curr_npoints = static_cast<int>(elem(npoints_offset, 0, 0));
      if (n+1 > curr_npoints) elem(npoints_offset, 0, 0) = n+1;
      elem(ndim_offset, 0, n) = ndim;
      for (int i = 0; i < ndim*(ndeg*nint + 1); ++i) elem(prof_offset, i, n) = prof(i);
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

void KNDataFile::getBlanket(int n, KNVector& blanket)
{
  const int curr_npoints = static_cast<int>(elem(npoints_offset, 0, 0));
  if ((blanket.size() == ndim*(ndeg1*nint1*ndeg2*nint2))&&(n < curr_npoints))
  {
    for (int i = 0; i < ndim*(ndeg1*nint1*ndeg2*nint2); ++i) blanket(i) = elem(blanket_offset, i, n);
  }
  else
  {
    P_MESSAGE1("Error while reading the MAT file.");
  }
}

// write in order size() = [ ndeg1*nint1, ndeg2*nint2, ndim ]
void KNDataFile::setBlanket(int n, const KNVector& blanket)
{
  if (wperm && n < ncols)
  {
    if (blanket.size() == ndim*(ndeg1*nint1*ndeg2*nint2))
    {
      const int curr_npoints = static_cast<int>(elem(npoints_offset, 0, 0));
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
      for (int i = 0; i < ndim*(ndeg1*nint1*ndeg2*nint2); ++i) elem(blanket_offset, i, n) = blanket(i);
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

void KNDataFile::getPar(int n, KNVector& par) const
{
  if (n < ncols)
  {
    if (par.size() <= npar)
    {
      for (int i = 0; i < par.size(); ++i)
        par(i) = elem(par_offset, i, n);
    }
    else
    {
      P_MESSAGE1("Error while reading the MAT file.");
    }
  }
  else
  {
    P_MESSAGE1("Error while reading the MAT file.");
  }
}

void KNDataFile::getMul(int n, KNVector& re, KNVector& im) const
{
  if (n < ncols && re.size() == im.size())
  {
    const int sz = std::min<int>(re.size(), nmul);
    for (int i = 0; i < sz; ++i)
    {
      re(i) = elem(mul_offset, i, n);
      im(i) = elem_im(mul_offset, i, n);
    }
  }
  else
  {
    P_MESSAGE1("Error while reading the MAT file.");
  }
}

void KNDataFile::getElem(int n, KNVector& el) const
{
  if (n < ncols)
  {
    if (el.size() == ndeg + 1)
    {
      for (int i = 0; i < ndeg + 1; ++i)
        el(i) = elem(elem_offset, i, n);
    }
    else
    {
      P_MESSAGE1("Error while reading the MAT file.");
    }
  }
  else
  {
    P_MESSAGE1("Error while reading the MAT file.");
  }
}

void KNDataFile::getMesh(int n, KNVector& mesh) const
{
  if (n < ncols)
  {
    if (mesh.size() == nint + 1)
    {
      for (int i = 0; i < nint + 1; ++i)
        mesh(i) = elem(mesh_offset, i, n);
    }
    else
    {
      P_MESSAGE1("Error while reading the MAT file.");
    }
  }
  else
  {
    P_MESSAGE1("Error while reading the MAT file.");
  }
}

void KNDataFile::getProfile(int n, KNVector& prof) const
{
  if (n < ncols)
  {
    if (prof.size() == ndim*(ndeg*nint + 1))
    {
      for (int i = 0; i < ndim*(ndeg*nint + 1); ++i)
        prof(i) = elem(prof_offset, i, n);
    }
    else
    {
      P_MESSAGE1("Error while reading the MAT file.");
    }
  }
  else
  {
    P_MESSAGE1("Error while reading the MAT file.");
  }
}

int KNDataFile::getUnstableMultipliers(int n) const
{
  KNVector mulRe(false), mulIm(false);
  const_cast<KNDataFile*>(this)->getMulReRef(n, mulRe);
  const_cast<KNDataFile*>(this)->getMulImRef(n, mulIm);
  const int lp = getNTrivMul(0);
  const int pd = getNTrivMul(1);
  const int ns = getNTrivMul(2);
  return unstableMultipliers(mulRe, mulIm, lp, pd, ns, n);
}

int KNDataFile::getNextBifurcation(int n, bool* stab) const
{
  KNVector mulRe(false), mulIm(false);
  const_cast<KNDataFile*>(this)->getMulReRef(n, mulRe);
  const_cast<KNDataFile*>(this)->getMulImRef(n, mulIm);
  const int lp = getNTrivMul(0);
  const int pd = getNTrivMul(1);
  const int ns = getNTrivMul(2);
  int p_ustab = unstableMultipliers(mulRe, mulIm, lp, pd, ns, n);
  for (int i = n + 1; i < getNPoints(); ++i)
  {
    const_cast<KNDataFile*>(this)->getMulReRef(i, mulRe);
    const_cast<KNDataFile*>(this)->getMulImRef(i, mulIm);
    int ustab = unstableMultipliers(mulRe, mulIm, lp, pd, ns, i);
    if (ustab != p_ustab) 
    {
      if (stab != 0)
      {
        if (ustab == 0 || p_ustab == 0) *stab = true;
        else *stab = false;
      }
      return i;
    }
    p_ustab = ustab;
  }
  return -1;
}

BifType  KNDataFile::getBifurcationType(int n) const
{
  if (n > 0)
  {
    KNVector mulReA(false), mulImA(false);
    KNVector mulReB(false), mulImB(false);
    const_cast<KNDataFile*>(this)->getMulReRef(n-1, mulReA);
    const_cast<KNDataFile*>(this)->getMulImRef(n-1, mulImA);
    const_cast<KNDataFile*>(this)->getMulReRef(n, mulReB);
    const_cast<KNDataFile*>(this)->getMulImRef(n, mulImB);

    const int lp = getNTrivMul(0);
    const int pd = getNTrivMul(1);
    const int ns = getNTrivMul(2);

    return bifurcationType(mulReA, mulImA, mulReB, mulImB, lp, pd, ns, n-1, n);
  }
  return BifNone;
}

int KNDataFile::findType(int32_t type, int n) const
{
  int found = 0;
  for (int i=0; i<getNCols(); ++i)
  {
    if (getMagic(i) == type) found++;
    if (found == n) return i;
  }
  return -1;
}
