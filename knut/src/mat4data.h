// ------------------------------------------------------------------------- //
//
// This is part of KNUT
// Copyright (c) 2006 by Robert Szalai
//
// For license, see the file COPYING in the root directory of the package
//
// ------------------------------------------------------------------------- //

#ifndef MAT4DATA_H
#define MAT4DATA_H

#include "matrix.h"
#include "pointtype.h"
#include <string>
#include <vector>
// HANDLE is defined here
#ifdef WIN32
#  include <windows.h>
#endif

#ifndef WIN32
#  include <sys/file.h>
#endif

// for the correct integer size
#ifndef _MSC_VER
# include <stdint.h>
#else
# define int32_t INT32
#endif

class KNDataFile
{
  public:

    // opens the file, determines its size, maps the memory from file, sets up variables
    // Constructor for periodic orbits
    KNDataFile(const std::string& fileName, const std::vector<std::string>& parNames, int steps_, int ndim_, int npar_, int nint_, int ndeg_, int nmul_);
    // Constructor for quasi-periodic orbits
    KNDataFile(const std::string& fileName, const std::vector<std::string>& parNames, int steps_, int ndim_, int npar_, int nint1_, int nint2_, int ndeg1_, int ndeg2_);
    // Constructor for opening an existing file
    KNDataFile(const std::string& fileName);
    // unmaps the memory, truncates the file if necessary, closes the file
    ~KNDataFile();
    
#ifndef WIN32
    void   lock() const { if (flock(file, LOCK_EX) != 0) std::cerr<<"Error locking\n"; }
    void   unlock() const { if (flock(file, LOCK_UN) != 0) std::cerr<<"Error locking\n"; }
#else
    // No locking on windows. Use a proper operating system instead.
    void   lock() const { }
    void   unlock() const { }
#endif
    const std::string& getFileName() const { return matFileName; }
    // resets the tables (no locking)
    void   initHeaders();
    void   setPar(int n, const KNVector& par);
    void   setParNames(const std::vector<std::string>& parNames);
    void   getParNames(std::vector<std::string>& parNames) const;
    void   setMul(int n, const KNVector& real, const KNVector& imag);
    void   setElem(int n, const KNVector& el);
    void   setMesh(int n, const KNVector& mesh);
    void   setProfile(int n, const KNVector& profile);
    void   setMesh1(int n, int j, double d)
    {
      elem(mesh1_offset, j, n) = d;
    }
    void   setMesh2(int n, int j, double d)
    {
      elem(mesh2_offset, j, n) = d;
    }
    void   getBlanket(int n, KNVector& blanket);
    void   setBlanket(int n, const KNVector& blanket);

    void   getPar(int n, KNVector& par) const;
    double getPar(int n, int j) const
    {
      return elem(par_offset, j, n);
    }
    void   getMul(int n, KNVector& real, KNVector& imag) const;
    double getMulRe(int n, int j) const
    {
      return elem(mul_offset, j, n);
    }
    double getMulIm(int n, int j) const
    {
      return elem_im(mul_offset, j, n);
    }
    void getMulReRef(int n, KNVector& el)
    {
      el.init(&elem(mul_offset, 0, n), nmul);
    }
    void getMulImRef(int n, KNVector& el)
    {
      el.init(&elem_im(mul_offset, 0, n), nmul);
    }
    // j == 0 : LP, j == 1 : PD, j == 2, NS
    int    getNTrivMul(const int j) const
    {
      return static_cast<int>(((double*)((char*)address + ntrivmul_offset + ntrivmul_header.col_off(0)))[j]);
    }
    void   setNTrivMul(const int j, int i)
    {
      ((double*)((char*)address + ntrivmul_offset + ntrivmul_header.col_off(0)))[j] = i;
    }
    void   getElem(int n, KNVector& el) const;
    double getElem(int n, int j) const
    {
      return elem(elem_offset, j, n);
    }
    void   getElemRef(int n, KNVector& el)
    {
      el.init(&elem(elem_offset, 0, n), ndeg + 1);
    }
    void   getMesh(int n, KNVector& mesh) const;
    double getMesh(int n, int j) const
    {
      return elem(mesh_offset, j, n);
    }
    void   getMeshRef(int n, KNVector& mesh)
    {
      mesh.init(&elem(mesh_offset, 0, n), nint + 1);
    }
    void   getProfile(int n, KNVector& profile) const;
    double getProfile(int n, int d, int j) const
    {
      return elem(prof_offset, d + ndim*j, n);
    }
    void   getProfileRef(int n, KNVector& profile)
    {
      profile.init(&elem(prof_offset, 0, n), ndim*(ndeg*nint + 1));
    }

    int  getNDim() const
    {
      return ndim;
    }
    int  getNInt() const
    {
      return nint;
    }
    int  getNDeg() const
    {
      return ndeg;
    }
    int  getNPar() const
    {
      return npar;
    }
    int  getNMul() const
    {
      return nmul;
    }
    int  getMeshLength() const
    {
      return nint + 1;
    }
    int  getElemLength() const
    {
      return ndeg + 1;
    }
    int  getNCols() const
    {
      return ncols;
    }
    int  getNPoints() const
    {
      return static_cast<int>(((double*)((char*)address + npoints_offset + npoints_header.col_off(0)))[0]);
    }
    int     getUnstableMultipliers(int n) const;
    int     getNextBifurcation(int n, bool* stab = 0) const;
    BifType getBifurcationType(int n) const;
    bool isTorus() const
    {
      return torus;
    }

  private:
    struct header
    {
      int32_t type;
      int32_t mrows;
      int32_t ncols;
      int32_t imagf;
      int32_t namelen;

      const char* name(const void* offset) const
      {
        return (char*)offset + namelen;
      }
      size_t size() const
      {
        return sizeof(struct header) + namelen;
      }
      size_t col_off(size_t i) const
      {
        return size() + i*mrows*sizeof(double);
      }
      size_t col_off_im(size_t i) const
      {
        return size() + (mrows*ncols + i*mrows)*sizeof(double);
      }
      size_t enddata(size_t i)
      {
        if (imagf == 0) return size() + i*mrows*sizeof(double);
        else return size() + 2*i*mrows*sizeof(double);
      }
      size_t getRows()
      {
        return (size_t)mrows;
      }
      size_t getCols()
      {
        return (size_t)ncols;
      }
    };
    static inline size_t createMatrixHeader(KNDataFile::header* hd, const char* name, int32_t rows, int32_t cols);
    static inline size_t createComplexMatrixHeader(KNDataFile::header* hd, const char* name, int32_t rows, int32_t cols);
    static inline void writeMatrixHeader(void* address, size_t offset, KNDataFile::header* hd, const char* name);

    off_t findMatrix(const char* name, KNDataFile::header* found, bool test=false, int32_t r=-1, int32_t c=-1, int32_t imag=-1, const char* fileName="");
    void openReadOnly(const std::string& fileName);

    void resizeMatrix(const char* name, int newcol);
    void condenseData();

    struct header *getHeader(size_t offset)
    {
      return (struct header*)((char*)address + offset);
    }
    const struct header *getHeader(size_t offset) const
      {
        return (struct header*)((char*)address + offset);
      }
    double& elem(size_t offset, size_t row, size_t col)
    {
      return ((double*)((char*)address + offset + getHeader(offset)->col_off(col)))[row];
    }
    const double& elem(size_t offset, size_t row, size_t col) const
    {
      return ((double*)((char*)address + offset + getHeader(offset)->col_off(col)))[row];
    }
    double& elem_im(size_t offset, size_t row, size_t col)
    {
      return ((double*)((char*)address + offset + getHeader(offset)->col_off_im(col)))[row];
    }
    const double& elem_im(size_t offset, size_t row, size_t col) const
    {
      return ((double*)((char*)address + offset + getHeader(offset)->col_off_im(col)))[row];
    }
    int32_t getRows(size_t offset) const
    {
      return getHeader(offset)->mrows;
    }
    int32_t getCols(size_t offset) const
    {
      return getHeader(offset)->ncols;
    }

#ifdef WIN32
    HANDLE file;
    HANDLE mapHandle;
#else
    int    file;
#endif
    const std::string matFileName;
    size_t filesize;
    void  *address;
    size_t size;
    int    ncols;
    bool   wperm;
    int    ndim;
    int    npar;
    int    nint;
    int    ndeg;
    int    nmul;

    int    nint1;
    int    nint2;
    int    ndeg1;
    int    ndeg2;

    bool   torus;

    off_t  npoints_offset;  //T
    header npoints_header;

    off_t  par_offset;      //T
    header par_header;
    
    off_t  parnames_offset;
    header parnames_header;

    off_t  mul_offset;
    header mul_header;

    off_t  ntrivmul_offset;
    header ntrivmul_header;

    off_t  ndim_offset;     //T
    header ndim_header;

    off_t  elem_offset;
    header elem_header;

    off_t  mesh_offset;
    header mesh_header;

    off_t  prof_offset;
    header prof_header;

    // for the torus
    off_t  nint1_offset;     //T
    header nint1_header;

    off_t  nint2_offset;     //T
    header nint2_header;

    off_t  ndeg1_offset;     //T
    header ndeg1_header;

    off_t  ndeg2_offset;     //T
    header ndeg2_header;

    off_t  mesh1_offset;     //T
    header mesh1_header;

    off_t  mesh2_offset;     //T
    header mesh2_header;

    off_t  blanket_offset;   //T
    header blanket_header;
};

#endif
