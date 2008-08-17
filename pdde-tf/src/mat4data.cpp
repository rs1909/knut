// ------------------------------------------------------------------------- //
//
// This is part of PDDE-CONT
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
#include "pderror.h"
#include "mat4data.h"
#include "multipliers.h"

#define MAT_BIG_ENDIAN 1000
#define MAT_LITTLE_ENDIAN 0

static int32_t byte_order()
{
  int32_t tmp = 0x01020304;
  if ( ((char*)(&tmp))[0] == 1 ) return MAT_BIG_ENDIAN;
  else if ( ((char*)(&tmp))[0] == 4 ) return MAT_LITTLE_ENDIAN;
  else P_MESSAGE( "Fatal error. Unrecognized byte order." );
}

int mat4Data::findMatrix(const char* name, mat4Data::header* found)
{
  struct header hd;
  int cur_off = 0;
  int cur_size;
  do
  {
    memcpy(&hd, (char*)address + cur_off, sizeof(struct header));
    P_ERROR_X1(hd.type == byte_order(), "Not a matrix of double precision elements.");
    if (hd.imagf == 0)
      cur_size = sizeof(struct header) + hd.namelen * sizeof(char) + hd.mrows * hd.ncols * sizeof(double);
    else
      cur_size = sizeof(struct header) + hd.namelen * sizeof(char) + 2 * hd.mrows * hd.ncols * sizeof(double);
    if (strncmp(name, (char*)address + cur_off + sizeof(struct header), 20) == 0)
    {
      memcpy(found, &hd, sizeof(struct header));
      return cur_off;
    }
    cur_off += cur_size;
  }
  while (cur_off < size);
  return -1;
}

#ifndef WIN32

static inline void *mmapFileWrite(int& file, const std::string& fileName, int size)
{
  if ((file = open(fileName.c_str(), O_CREAT | O_RDWR | O_TRUNC, S_IRUSR | S_IWUSR)) == -1)
  {
    P_ERROR_X2(false, "Unable to create the MAT file.", strerror(errno));
  }

  if (ftruncate(file, size) != 0)
  {
    P_ERROR_X2(false, "Unable to resize the MAT file.", strerror(errno));
  }

  void *address;
  if ((address = mmap(0, size, PROT_WRITE | PROT_READ, MAP_SHARED, file, 0)) == MAP_FAILED)
  {
    P_ERROR_X2(false, "Unable to map the MAT file to a memory location.", strerror(errno));
  }
  return address;
}

static inline void *mmapFileRead(int& file, const std::string& fileName, int& size)
{
  if ((file = open(fileName.c_str(), O_RDONLY)) == -1)
  {
    P_ERROR_X2(false, "Unable to open the MAT file for reading.", strerror(errno));
  }

  struct stat filestat;
  if (fstat(file, &filestat) != 0)
  {
    P_ERROR_X2(false, "Unable to stat the MAT file.", strerror(errno));
  }
  int filesize = filestat.st_size;
  size = filesize;

  void *address;
  if ((address = mmap(0, filesize, PROT_READ, MAP_PRIVATE, file, 0)) == MAP_FAILED)
  {
    P_ERROR_X2(false, "Unable to map the MAT file to a memory location.", strerror(errno));
  }
  return address;
}

#else

static inline void *mmapFileWrite(HANDLE& file, HANDLE& mapHandle, const std::string& fileName, int size)
{
  if ((file = CreateFile(fileName.c_str(),
                         FILE_WRITE_DATA | FILE_READ_DATA,
                         FILE_SHARE_READ,
                         NULL,
                         CREATE_ALWAYS,
                         FILE_ATTRIBUTE_NORMAL,
                         NULL)) == NULL)
  {
    P_ERROR_X2(false, "Unable to create the MAT file.", static_cast<int>(GetLastError()));
  }

  if (SetFilePointer(file, size, NULL, FILE_BEGIN) == 0)
  {
    P_ERROR_X2(false, "Unable to seek in the MAT file.", static_cast<int>(GetLastError()));
  }
  P_ERROR_X2(SetEndOfFile(file), "Unable to truncate the MAT.", static_cast<int>(GetLastError()));

  if ((mapHandle = CreateFileMapping(file, NULL, PAGE_READWRITE, 0, size, fileName.c_str())) == 0)
  {
    P_ERROR_X2(false, "Unable to map the MAT file to a memory location.", static_cast<int>(GetLastError()));
  }

  void *address = MapViewOfFile(mapHandle, FILE_MAP_WRITE, 0, 0, 0);
  if (address != NULL) return address;
  else P_ERROR_X2(false, "Unable to view the file map.", static_cast<int>(GetLastError()));
  return 0;
}

static inline void *mmapFileRead(HANDLE& file, HANDLE& mapHandle, const std::string& fileName, int& size)
{
  if ((file = CreateFile(fileName.c_str(),
                         FILE_READ_DATA,
                         FILE_SHARE_READ,
                         NULL,
                         OPEN_EXISTING,
                         FILE_ATTRIBUTE_NORMAL,
                         NULL)) == NULL)
  {
    P_ERROR_X2(false, "Unable to create the MAT file.", static_cast<int>(GetLastError()));
  }

  size = GetFileSize(file, NULL);
  if ((mapHandle = CreateFileMapping(file, NULL, PAGE_READONLY, 0, 0, fileName.c_str())) == 0)
  {
    P_ERROR_X2(false, "Unable to map the MAT file to a memory location.", static_cast<int>(GetLastError()));
  }

  void *address = MapViewOfFile(mapHandle, FILE_MAP_READ, 0, 0, 0);
  if (address != NULL) return address;
  else P_ERROR_X2(false, "Unable to view the file map.", static_cast<int>(GetLastError()));
  return 0;
}

#endif

// returns the size of the whole data
inline int mat4Data::createMatrixHeader(mat4Data::header* hd, const char* name, int rows, int cols)
{
  hd->type = byte_order();
  hd->mrows = rows;
  hd->ncols = cols;
  hd->imagf = 0;
  hd->namelen = ((strlen(name) + sizeof(mat4Data::header) + 1) / sizeof(double) + 1) * sizeof(double) - sizeof(mat4Data::header);
  return sizeof(mat4Data::header) + hd->namelen + rows * cols * sizeof(double);
}

// returns the size of the whole data
inline int mat4Data::createComplexMatrixHeader(mat4Data::header* hd, const char* name, int rows, int cols)
{
  hd->type = byte_order();
  hd->mrows = rows;
  hd->ncols = cols;
  hd->imagf = 1;
  hd->namelen = ((strlen(name) + sizeof(mat4Data::header) + 1) / sizeof(double) + 1) * sizeof(double) - sizeof(mat4Data::header);
  return sizeof(mat4Data::header) + hd->namelen + 2 * rows * cols * sizeof(double);
}

inline int mat4Data::writeMatrixHeader(void* address, int offset, mat4Data::header* hd, const char* name)
{
  strncpy((char*)address + offset + sizeof(mat4Data::header), name, hd->namelen);
  *((mat4Data::header*)((char*)address + offset)) = *hd;
}

mat4Data::mat4Data(const std::string& fileName, int steps_, int ndim_, int npar_, int nint_, int ndeg_, int nmul_)
{
  wperm = true;
  torus = false;

  ncols = steps_;
  ndim = ndim_;
  npar = npar_;
  nint = nint_;
  ndeg = ndeg_;
  nmul = nmul_;

  // creating the matrices
  char npoints_string[] = "pdde_npoints";
  npoints_offset = 0;
  int npoints_size = createMatrixHeader(&npoints_header, npoints_string, 1, 1);

  char par_string[] = "pdde_par";
  par_offset = npoints_offset + npoints_size;
  int par_size = createMatrixHeader(&par_header, par_string, npar, ncols);

  char mul_string[] = "pdde_mul";
  mul_offset = par_offset + par_size;
  int mul_size = createComplexMatrixHeader(&mul_header, mul_string, nmul, ncols);

  char ntrivmul_string[] = "pdde_ntrivmul";
  ntrivmul_offset = mul_offset + mul_size;
  int ntrivmul_size = createMatrixHeader(&ntrivmul_header, ntrivmul_string, 3, 1);

  char ndim_string[] = "pdde_ndim";
  ndim_offset = ntrivmul_offset + ntrivmul_size;
  int ndim_size = createMatrixHeader(&ndim_header, ndim_string, 1, ncols);

  char elem_string[] = "pdde_elem";
  elem_offset = ndim_offset + ndim_size;
  int elem_size = createMatrixHeader(&elem_header, elem_string, ndeg + 1, ncols);

  char mesh_string[] = "pdde_mesh";
  mesh_offset = elem_offset + elem_size;
  int mesh_size = createMatrixHeader(&mesh_header, mesh_string, nint + 1, ncols);

  char prof_string[] = "pdde_prof";
  prof_offset = mesh_offset + mesh_size;
  int prof_size = createMatrixHeader(&prof_header, prof_string, ndim * (ndeg * nint + 1), ncols);
  size = prof_offset + prof_size;

  const int approxSize = 8 * (sizeof(header) + 20) + sizeof(double) * (1 + ncols * (npar + 2 * nmul + 1 + (ndeg + 1) + (ndim + 1) * (ndeg * nint + 1)));
#ifndef WIN32
  address = mmapFileWrite(file, fileName, size);
#else
  address = mmapFileWrite(file, mapHandle, fileName, size);
#endif

  writeMatrixHeader(address, npoints_offset, &npoints_header, npoints_string);
  writeMatrixHeader(address, par_offset, &par_header, par_string);
  writeMatrixHeader(address, mul_offset, &mul_header, mul_string);
  writeMatrixHeader(address, ntrivmul_offset, &ntrivmul_header, ntrivmul_string);
  writeMatrixHeader(address, ndim_offset, &ndim_header, ndim_string);
  writeMatrixHeader(address, elem_offset, &elem_header, elem_string);
  writeMatrixHeader(address, mesh_offset, &mesh_header, mesh_string);
  writeMatrixHeader(address, prof_offset, &prof_header, prof_string);
}

mat4Data::mat4Data(const std::string& fileName, int steps_, int ndim_, int npar_, int nint1_, int nint2_, int ndeg1_, int ndeg2_)
{
  wperm = true;
  torus = true;

  ncols = steps_;
  ndim = ndim_;
  npar = npar_;
  nint1 = nint1_;
  nint2 = nint2_;
  ndeg1 = ndeg1_;
  ndeg2 = ndeg2_;

  // creating the matrices
  char npoints_string[] = "pdde_npoints";
  npoints_offset = 0;
  int npoints_size = createMatrixHeader(&npoints_header, npoints_string, 1, 1);

  char par_string[] = "pdde_par";
  par_offset = npoints_offset + npoints_size;
  int par_size = createMatrixHeader(&par_header, par_string, npar, ncols);

  char ndim_string[] = "pdde_ndim";
  ndim_offset = par_offset + par_size;
  int ndim_size = createMatrixHeader(&ndim_header, ndim_string, 1, ncols);

  char nint1_string[] = "pdde_nint1";
  nint1_offset = ndim_offset + ndim_size;
  int nint1_size = createMatrixHeader(&nint1_header, nint1_string, 1, ncols);

  char nint2_string[] = "pdde_nint2";
  nint2_offset = nint1_offset + nint1_size;
  int nint2_size = createMatrixHeader(&nint2_header, nint2_string, 1, ncols);

  char ndeg1_string[] = "pdde_ndeg1";
  ndeg1_offset = nint2_offset + nint2_size;
  int ndeg1_size = createMatrixHeader(&ndeg1_header, ndeg1_string, 1, ncols);

  char ndeg2_string[] = "pdde_ndeg2";
  ndeg2_offset = ndeg1_offset + ndeg1_size;
  int ndeg2_size = createMatrixHeader(&ndeg2_header, ndeg2_string, 1, ncols);

  char mesh1_string[] = "pdde_mesh1";
  mesh1_offset = ndeg2_offset + ndeg2_size;
  int mesh1_size = createMatrixHeader(&mesh1_header, mesh1_string, nint1 * ndeg1, ncols);

  char mesh2_string[] = "pdde_mesh2";
  mesh2_offset = mesh1_offset + mesh1_size;
  int mesh2_size = createMatrixHeader(&mesh2_header, mesh2_string, nint2 * ndeg2, ncols);

  char blanket_string[] = "pdde_blanket";
  blanket_offset = mesh2_offset + mesh2_size;
  int blanket_size = createMatrixHeader(&blanket_header, blanket_string, ndim * nint1 * ndeg1 * nint2 * ndeg2, ncols);
  size = blanket_offset + blanket_size;

#ifndef WIN32
  address = mmapFileWrite(file, fileName, size);
#else
  address = mmapFileWrite(file, mapHandle, fileName, size);
#endif
  writeMatrixHeader(address, npoints_offset, &npoints_header, npoints_string);
  writeMatrixHeader(address, par_offset, &par_header, par_string);
  writeMatrixHeader(address, ndim_offset, &ndim_header, ndim_string);
  writeMatrixHeader(address, nint1_offset, &nint1_header, nint1_string);
  writeMatrixHeader(address, nint2_offset, &nint2_header, nint2_string);
  writeMatrixHeader(address, ndeg1_offset, &ndeg1_header, ndeg1_string);
  writeMatrixHeader(address, ndeg2_offset, &ndeg2_header, ndeg2_string);
  writeMatrixHeader(address, mesh1_offset, &mesh1_header, mesh1_string);
  writeMatrixHeader(address, mesh2_offset, &mesh2_header, mesh2_string);
  writeMatrixHeader(address, blanket_offset, &blanket_header, blanket_string);
}

mat4Data::mat4Data(const std::string& fileName)
{
  wperm = false;
  this->openReadOnly(fileName);
}

void mat4Data::resizeMatrix(const char* name, int newcol)
{
  header mathead;
  int matoffset;
  if ((matoffset = findMatrix(name, &mathead)) == -1) P_MESSAGE("No such matrix");
  P_ASSERT_X1(mathead.ncols >= newcol, "Cannot increase the size of a matrix");
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

void mat4Data::condenseData()
{
  // condense the fields so that it is only the size of npoints, where it counts
  if (wperm)
  {
    const int npoints = getNPoints();
    if (!torus)
    {
      resizeMatrix("pdde_prof", npoints);
      resizeMatrix("pdde_mesh", npoints);
      resizeMatrix("pdde_elem", npoints);
      resizeMatrix("pdde_ndim", npoints);
      resizeMatrix("pdde_mul", npoints);
      resizeMatrix("pdde_par", npoints);
    }
  }
}

void mat4Data::openReadOnly(const std::string& fileName)
{
#ifndef WIN32
  address = mmapFileRead(file, fileName, size);
#else
  address = mmapFileRead(file, mapHandle, fileName, size);
#endif

  if ((npoints_offset = findMatrix("pdde_npoints", &npoints_header)) == -1) P_MESSAGE("err-4");
  P_ERROR_X1(npoints_header.mrows == 1, "err-3");
  P_ERROR_X1(npoints_header.ncols == 1, "err-2");
  P_ERROR_X1(npoints_header.imagf == 0, "err-1");

  if ((par_offset = findMatrix("pdde_par", &par_header)) == -1) P_MESSAGE("err1");
  npar = par_header.mrows;
  ncols = par_header.ncols;
  if (par_header.imagf != 0) P_MESSAGE("err2");

  if ((ndim_offset = findMatrix("pdde_ndim", &ndim_header)) == -1) P_MESSAGE("err6");
  if (ndim_header.mrows != 1) P_MESSAGE("err7 ");
  if (ndim_header.ncols != ncols) P_MESSAGE("err9");
  if (ndim_header.imagf != 0) P_MESSAGE("err9");
  ndim = static_cast<int>(*((double*)((char*)address + ndim_offset + ndim_header.col_off(0))));

  if ((mul_offset = findMatrix("pdde_mul", &mul_header)) != -1)
  {
    torus = false;
    // periodic solutions
    nmul = mul_header.mrows;
    if (mul_header.ncols != ncols) P_MESSAGE("err4");
    if (mul_header.imagf != 1) P_MESSAGE("err5");

    if ((ntrivmul_offset = findMatrix("pdde_ntrivmul", &ntrivmul_header)) == -1) P_MESSAGE("err6");
    P_ERROR_X1(ntrivmul_header.mrows == 3, "err7");
    P_ERROR_X1(ntrivmul_header.ncols == 1, "err8");
    P_ERROR_X1(ntrivmul_header.imagf == 0, "err9");

    if ((elem_offset = findMatrix("pdde_elem", &elem_header)) == -1) P_MESSAGE("err18.1");
    if (elem_header.ncols != ncols) P_MESSAGE("err20.1");
    if (elem_header.imagf != 0) P_MESSAGE("err21.1");
    ndeg = elem_header.mrows - 1;

    if ((mesh_offset = findMatrix("pdde_mesh", &mesh_header)) == -1) P_MESSAGE("err18");
    if (mesh_header.ncols != ncols) P_MESSAGE("err20");
    if (mesh_header.imagf != 0) P_MESSAGE("err21");
    nint = mesh_header.mrows - 1;

    if ((prof_offset = findMatrix("pdde_prof", &prof_header)) == -1) P_MESSAGE("err22");
    if (prof_header.mrows != ndim*(ndeg*nint + 1)) P_MESSAGE("err23");
    if (prof_header.ncols != ncols) P_MESSAGE("err24");
    if (prof_header.imagf != 0) P_MESSAGE("err25");
  }
  else
  {
    torus = true;
    // quasiperiodic solutions
    if ((nint1_offset = findMatrix("pdde_nint1", &nint1_header)) == -1) P_MESSAGE("err26");
    if (nint1_header.mrows != 1) P_MESSAGE("err27");
    if (nint1_header.ncols != ncols) P_MESSAGE("err28");
    if (nint1_header.imagf != 0) P_MESSAGE("err29");
    nint1 = static_cast<int>(*((double*)((char*)address + nint1_offset + nint1_header.col_off(0))));

    if ((nint2_offset = findMatrix("pdde_nint2", &nint2_header)) == -1) P_MESSAGE("err30");
    if (nint2_header.mrows != 1) P_MESSAGE("err31");
    if (nint2_header.ncols != ncols) P_MESSAGE("err32");
    if (nint2_header.imagf != 0) P_MESSAGE("err33");
    nint2 = static_cast<int>(*((double*)((char*)address + nint2_offset + nint2_header.col_off(0))));

    if ((ndeg1_offset = findMatrix("pdde_ndeg2", &ndeg1_header)) == -1) P_MESSAGE("err34");
    if (ndeg1_header.mrows != 1) P_MESSAGE("err35");
    if (ndeg1_header.ncols != ncols) P_MESSAGE("err36");
    if (ndeg1_header.imagf != 0) P_MESSAGE("err37");
    ndeg1 = static_cast<int>(*((double*)((char*)address + ndeg1_offset + ndeg1_header.col_off(0))));

    if ((ndeg2_offset = findMatrix("pdde_ndeg2", &ndeg2_header)) == -1) P_MESSAGE("err38");
    if (ndeg2_header.mrows != 1) P_MESSAGE("err39");
    if (ndeg2_header.ncols != ncols) P_MESSAGE("err40");
    if (ndeg2_header.imagf != 0) P_MESSAGE("err41");
    ndeg2 = static_cast<int>(*((double*)((char*)address + ndeg2_offset + ndeg2_header.col_off(0))));

    if ((mesh1_offset = findMatrix("pdde_mesh1", &mesh1_header)) == -1) P_MESSAGE("err42");
    if (mesh1_header.mrows != nint1*ndeg1) P_MESSAGE("err43");
    if (mesh1_header.ncols != ncols) P_MESSAGE("err44");
    if (mesh1_header.imagf != 0) P_MESSAGE("err45");

    if ((mesh2_offset = findMatrix("pdde_mesh2", &mesh2_header)) == -1) P_MESSAGE("err46");
    if (mesh2_header.mrows != nint2*ndeg2) P_MESSAGE("err47");
    if (mesh2_header.ncols != ncols) P_MESSAGE("err48");
    if (mesh2_header.imagf != 0) P_MESSAGE("err49");

    if ((blanket_offset = findMatrix("pdde_blanket", &blanket_header)) == -1) P_MESSAGE("err50");
    if (blanket_header.mrows != ndim*nint1*ndeg1*nint2*ndeg2) P_MESSAGE("err51");
    if (blanket_header.ncols != ncols) P_MESSAGE("err52");
    if (blanket_header.imagf != 0) P_MESSAGE("err53");
  }
}

// WARNING This creates a large file and then shrinks it
// Use mmremap instead so that only the smaller file is written.
mat4Data::~mat4Data()
{
#ifndef WIN32
  int oldsize = size;
  condenseData();
  if (munmap(address, oldsize) != 0)
  {
    P_ERROR_X2(false, "Unable to munmap the MAT file.", strerror(errno));
  }
  // at the moment we try to truncate the file
  if( wperm )
  {
    if( ftruncate( file, size) != 0 )
    {
      P_ERROR_X2( false, "Unable to truncate the MAT file.", strerror( errno ) );
    }
  }
  if (close(file) != 0)
  {
    P_ERROR_X2(false, "Unable to close the MAT file.", strerror(errno));
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
}

void mat4Data::setPar(int n, const Vector& par)
{
  if (wperm && n < ncols)
  {
    if (par.Size() <= npar)
    {
      for (int i = 0; i < par.Size(); ++i)
        elem(par_offset, i, n) = par(i);
      for (int i = par.Size(); i < npar; ++i)
        elem(par_offset, i, n) = 0.0;
    }
    else
    {
      P_MESSAGE("Error while writing the MAT file.");
    }
  }
  else
  {
    P_MESSAGE("Error while writing the MAT file.");
  }
}

void mat4Data::setMul(int n, const Vector& re, const Vector& im)
{
  if (wperm && n < ncols && re.Size() == im.Size())
  {
    if (re.Size() <= nmul)
    {
      for (int i = 0; i < re.Size(); ++i)
      {
        elem(mul_offset, i, n) = re(i);
        elem_im(mul_offset, i, n) = im(i);
      }
      for (int i = re.Size(); i < nmul; ++i)
      {
        elem(mul_offset, i, n) = 0.0;
        elem_im(mul_offset, i, n) = 0.0;
      }
    }
    else
    {
      P_MESSAGE("Error while writing the MAT file.");
    }
  }
  else
  {
    P_MESSAGE("Error while writing the MAT file.");
  }
}

void mat4Data::setElem(int n, const Vector& el)
{
  if (wperm && n < ncols)
  {
    if (el.Size() == ndeg + 1)
    {
      for (int i = 0; i < ndeg + 1; ++i)
        elem(elem_offset, i, n) = el(i);
    }
    else
    {
      P_MESSAGE("Error while writing the MAT file.");
    }
  }
  else
  {
    P_MESSAGE("Error while writing the MAT file.");
  }
}

void mat4Data::setMesh(int n, const Vector& mesh)
{
  if (wperm && n < ncols)
  {
    if (mesh.Size() == nint + 1)
    {
      for (int i = 0; i < nint + 1; ++i)
        elem(mesh_offset, i, n) = mesh(i);
    }
    else
    {
      P_MESSAGE("Error while writing the MAT file.");
    }
  }
  else
  {
    P_MESSAGE("Error while writing the MAT file.");
  }
}

void mat4Data::setProfile(int n, const Vector& prof)
{
  if (wperm && n < ncols)
  {
    if (prof.Size() == ndim*(ndeg*nint + 1))
    {
      const int curr_npoints = static_cast<int>(elem(npoints_offset, 0, 0));
      if (n+1 > curr_npoints) elem(npoints_offset, 0, 0) = n+1;
      elem(ndim_offset, 0, n) = ndim;
      for (int i = 0; i < ndim*(ndeg*nint + 1); ++i) elem(prof_offset, i, n) = prof(i);
    }
    else
    {
      P_MESSAGE("Error while writing the MAT file.");
    }
  }
  else
  {
    P_MESSAGE("Error while writing the MAT file.");
  }
}

void mat4Data::getBlanket(int n, Vector& blanket)
{
  const int curr_npoints = static_cast<int>(elem(npoints_offset, 0, 0));
  if ((blanket.Size() == ndim*(ndeg1*nint1*ndeg2*nint2))&&(n < curr_npoints))
  {
    for (int i = 0; i < ndim*(ndeg1*nint1*ndeg2*nint2); ++i) blanket(i) = elem(blanket_offset, i, n);
  }
  else
  {
    P_MESSAGE("Error while reading the MAT file.");
  }
}

// write in order size() = [ ndeg1*nint1, ndeg2*nint2, ndim ]
void mat4Data::setBlanket(int n, const Vector& blanket)
{
  if (wperm && n < ncols)
  {
    if (blanket.Size() == ndim*(ndeg1*nint1*ndeg2*nint2))
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
      P_MESSAGE("Error while writing the MAT file.");
    }
  }
  else
  {
    P_MESSAGE("Error while writing the MAT file.");
  }
}

void mat4Data::getPar(int n, Vector& par) const
{
  if (n < ncols)
  {
    if (par.Size() <= npar)
    {
      for (int i = 0; i < par.Size(); ++i)
        par(i) = elem(par_offset, i, n);
    }
    else
    {
      P_MESSAGE("Error while reading the MAT file.");
    }
  }
  else
  {
    P_MESSAGE("Error while reading the MAT file.");
  }
}

void mat4Data::getMul(int n, Vector& re, Vector& im) const
{
  if (n < ncols && re.Size() == im.Size())
  {
    const int sz = std::min<int>(re.Size(), nmul);
    for (int i = 0; i < sz; ++i)
    {
      re(i) = elem(mul_offset, i, n);
      im(i) = elem_im(mul_offset, i, n);
    }
  }
  else
  {
    P_MESSAGE("Error while reading the MAT file.");
  }
}

void mat4Data::getElem(int n, Vector& el) const
{
  if (n < ncols)
  {
    if (el.Size() == ndeg + 1)
    {
      for (int i = 0; i < ndeg + 1; ++i)
        el(i) = elem(elem_offset, i, n);
    }
    else
    {
      P_MESSAGE("Error while reading the MAT file.");
    }
  }
  else
  {
    P_MESSAGE("Error while reading the MAT file.");
  }
}

void mat4Data::getMesh(int n, Vector& mesh) const
{
  if (n < ncols)
  {
    if (mesh.Size() == nint + 1)
    {
      for (int i = 0; i < nint + 1; ++i)
        mesh(i) = elem(mesh_offset, i, n);
    }
    else
    {
      P_MESSAGE("Error while reading the MAT file.");
    }
  }
  else
  {
    P_MESSAGE("Error while reading the MAT file.");
  }
}

void mat4Data::getProfile(int n, Vector& prof) const
{
  if (n < ncols)
  {
    if (prof.Size() == ndim*(ndeg*nint + 1))
    {
      for (int i = 0; i < ndim*(ndeg*nint + 1); ++i)
        prof(i) = elem(prof_offset, i, n);
    }
    else
    {
      P_MESSAGE("Error while reading the MAT file.");
    }
  }
  else
  {
    P_MESSAGE("Error while reading the MAT file.");
  }
}

int mat4Data::getNextBifurcation(int n) const
{
  Vector mulRe(false), mulIm(false);
  const_cast<mat4Data*>(this)->getMulReRef(n, mulRe);
  const_cast<mat4Data*>(this)->getMulImRef(n, mulIm);
  const int lp = getNTrivMul(0);
  const int pd = getNTrivMul(1);
  const int ns = getNTrivMul(2);
  int p_ustab = unstableMultipliers(mulRe, mulIm, lp, pd, ns);
  for (int i = n + 1; i < getNPoints(); ++i)
  {
    const_cast<mat4Data*>(this)->getMulReRef(i, mulRe);
    const_cast<mat4Data*>(this)->getMulImRef(i, mulIm);
    int ustab = unstableMultipliers(mulRe, mulIm, lp, pd, ns);
    if (ustab != p_ustab) return i;
    p_ustab = ustab;
  }
  return -1;
}

PtType  mat4Data::getBifurcationType(int n) const
{
  Vector mulRe(false), mulIm(false);
  const_cast<mat4Data*>(this)->getMulReRef(n, mulRe);
  const_cast<mat4Data*>(this)->getMulImRef(n, mulIm);
  const int lp = getNTrivMul(0);
  const int pd = getNTrivMul(1);
  const int ns = getNTrivMul(2);
  return bifurcationType(mulRe, mulIm, lp, pd, ns);
}
