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
#include <cstddef>

// Implement a global lock: all file operations will block execution even if they access different files
#include <mutex>

// File mapping API
#include <boost/interprocess/file_mapping.hpp>
#include <boost/interprocess/mapped_region.hpp>

// for the correct integer size
// #ifndef _WIN32
# include <cstdint>
// #else
// # define int32_t INT32
// #endif

// TODO: make it possible to release data once it is not needed, e.g., after reading
// TODO: or just delete the data in 'basecomp.cpp'

using namespace boost::interprocess;

class KNAbstractData
{
  public:
    virtual void lockRead() const = 0;
    virtual void lockWrite() const = 0;
    virtual void unlock() const = 0;
    virtual void   setNTrivMul(const size_t j, size_t i) = 0;
    virtual void   setMagic(size_t n, BifType magic) = 0;
    virtual void   setPar(size_t n, const KNVector& par) = 0;
    virtual void   setMul(size_t n, const KNVector& real, const KNVector& imag) = 0;
    virtual void   setElem(size_t n, const KNVector& el) = 0;
    virtual void   setMesh(size_t n, const KNVector& mesh) = 0;
    virtual void   setProfile(size_t n, const KNVector& profile) = 0;
    virtual void   getPar(size_t n, KNVector& par) const = 0;
    virtual void   getMul(size_t n, KNVector& real, KNVector& imag) const = 0;
    virtual void   getMesh(size_t n, KNVector& mesh) const = 0;
    virtual size_t  getNPar() const = 0;
    virtual size_t  getNDim() const = 0;
    virtual size_t  getNInt() const = 0;
    virtual size_t  getNDeg() const = 0;
    virtual size_t  getNCols() const = 0;
    virtual BifType getMagic(size_t n) const = 0;
    virtual void   getProfile(size_t n, KNVector& profile) const = 0;
    size_t findType(BifType type, size_t n) const;
    // the extra bit for tori
    virtual void   setMesh1(size_t n, size_t j, double d) = 0;
    virtual void   setMesh2(size_t n, size_t j, double d) = 0;
    virtual void   getBlanket(size_t n, KNVector& blanket) const = 0;
    virtual void   setBlanket(size_t n, const KNVector& blanket) = 0;
};

class KNDataFile : public KNAbstractData
{
  public:

    // opens the file, determines its size, maps the memory from file, sets up variables
    // Constructor for periodic orbits
    KNDataFile(const std::string& fileName, const std::vector<std::string>& parNames, size_t steps_, size_t ndim_, size_t npar_, size_t nint_, size_t ndeg_, size_t nmul_);
    // Constructor for quasi-periodic orbits
    KNDataFile(const std::string& fileName, const std::vector<std::string>& parNames, size_t steps_, size_t ndim_, size_t npar_, size_t nint1_, size_t nint2_, size_t ndeg1_, size_t ndeg2_);
    // Constructor for opening an existing file
    KNDataFile(const std::string& fileName);
    // unmaps the memory, truncates the file if necessary, closes the file
    virtual ~KNDataFile();
    
    void lockRead() const override;
    void lockWrite() const override;
    void unlock() const override;

    const std::string& getFileName() const { return matFileName; }
    // resets the tables (no locking)
    void   initHeaders();
    void   setPar(size_t n, const KNVector& par) override;
    void   setParNames(const std::vector<std::string>& parNames);
    void   getParNames(std::vector<std::string>& parNames) const;
    void   setMul(size_t n, const KNVector& real, const KNVector& imag) override;
    void   setElem(size_t n, const KNVector& el) override;
    void   setMesh(size_t n, const KNVector& mesh) override;
    void   setProfile(size_t n, const KNVector& profile) override;
    void   setMesh1(size_t n, size_t j, double d) override
    {
      elem(mesh1_offset, j, n) = d;
    }
    void   setMesh2(size_t n, size_t j, double d) override
    {
      elem(mesh2_offset, j, n) = d;
    }
    void   getBlanket(size_t n, KNVector& blanket) const override;
    void   setBlanket(size_t n, const KNVector& blanket) override;
    void   setMagic(size_t n, BifType magic) override
    {
      elem(magic_offset, 0, n) = static_cast<double>(magic);
    }
    BifType getMagic(size_t n) const override
    {
      return static_cast<BifType>(elem(magic_offset, 0, n));
    }
    void   getPar(size_t n, KNVector& par) const override;
    double getPar(size_t n, size_t j) const
    {
      return elem(par_offset, j, n);
    }
    void   getMul(size_t n, KNVector& real, KNVector& imag) const override;
    double getMulRe(size_t n, size_t j) const
    {
      return elem(mul_offset, j, n);
    }
    double getMulIm(size_t n, size_t j) const
    {
      return elem_im(mul_offset, j, n);
    }
    void getMulReRef(size_t n, KNVector& el)
    {
      el.init(&elem(mul_offset, 0, n), nmul);
    }
    void getMulImRef(size_t n, KNVector& el)
    {
      el.init(&elem_im(mul_offset, 0, n), nmul);
    }
    // j == 0 : LP, j == 1 : PD, j == 2, NS
    size_t getNTrivMul(const size_t j) const
    {
      return static_cast<size_t>(((double*)((char*)address + ntrivmul_offset + ntrivmul_header.col_off(0)))[j]);
    }
    void   setNTrivMul(const size_t j, size_t i) override
    {
      ((double*)((char*)address + ntrivmul_offset + ntrivmul_header.col_off(0)))[j] = i;
    }
    void   getElem(size_t n, KNVector& el) const;
    double getElem(size_t n, size_t j) const
    {
      return elem(elem_offset, j, n);
    }
    void   getElemRef(size_t n, KNVector& el)
    {
      el.init(&elem(elem_offset, 0, n), ndeg + 1);
    }
    void   getMesh(size_t n, KNVector& mesh) const override;
    double getMesh(size_t n, size_t j) const
    {
      return elem(mesh_offset, j, n);
    }
    void   getMeshRef(size_t n, KNVector& mesh)
    {
      mesh.init(&elem(mesh_offset, 0, n), nint + 1);
    }
    void   getProfile(size_t n, KNVector& profile) const override;
    double getProfile(size_t n, size_t d, size_t j) const
    {
      return elem(prof_offset, d + ndim*j, n);
    }
    void   getProfileRef(size_t n, KNVector& profile)
    {
      profile.init(&elem(prof_offset, 0, n), ndim*(ndeg*nint + 1));
    }

    size_t  getNDim() const override
    {
      return ndim;
    }
    size_t  getNInt() const override
    {
      return nint;
    }
    size_t  getNDeg() const override
    {
      return ndeg;
    }
    size_t  getNPar() const override
    {
      return npar;
    }
    size_t  getNMul() const
    {
      return nmul;
    }
    size_t  getMeshLength() const
    {
      return nint + 1;
    }
    size_t  getElemLength() const
    {
      return ndeg + 1;
    }
    size_t  getNCols() const override
    {
      return ncols;
    }
    size_t  getNPoints() const
    {
      return static_cast<size_t>(((double*)((char*)address + npoints_offset + npoints_header.col_off(0)))[0]);
    }
    size_t  getUnstableMultipliers(size_t n) const;
    size_t  getNextBifurcation(size_t n, bool* stab = nullptr) const;
    BifType getBifurcationType(size_t n) const;
    // moved to Abstract subclass
//     size_t  findType(int32_t type, size_t n) const; // to find the n-th point of type; returns index
    bool isTorus() const
    {
      return torus;
    }

  private:
    struct header
    {
      int32_t type;
      uint32_t mrows;
      uint32_t ncols;
      uint32_t imagf;
      uint32_t namelen;

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
    static inline size_t createMatrixHeader(KNDataFile::header* hd, const char* name, size_t rows, size_t cols, int32_t type = 0);
    static inline size_t createComplexMatrixHeader(KNDataFile::header* hd, const char* name, size_t rows, size_t cols, int32_t type = 0);
    static inline void writeMatrixHeader(void* address, size_t offset, KNDataFile::header* hd, const char* name);

    size_t findMatrix(const char* name, KNDataFile::header* found, bool test=false, size_t r=0, size_t c=0, uint32_t imag=0, const char* fileName="", int32_t type=0);
    void openReadOnly(const std::string& fileName);

    void resizeMatrix(const char* name, size_t newcol, uint32_t imag, int32_t type=0);
    void condenseData();

    struct header *getHeader(size_t offset)
    {
      P_ERROR_X1(offset + sizeof(struct header) < size, "Matrix header is out of bound.");
      return (struct header*)((char*)address + offset);
    }
    const struct header *getHeader(size_t offset) const
    {
      P_ERROR_X1(offset + sizeof(struct header) < size, "Matrix header is out of bound.");
      return (struct header*)((char*)address + offset);
    }
    double& elem(size_t offset, size_t row, size_t col)
    {
      P_ERROR_X1(offset + getHeader(offset)->col_off(col) + sizeof(double)*row < size, "Matrix element is out of bound.");
      return ((double*)((char*)address + offset + getHeader(offset)->col_off(col)))[row];
    }
    double elem(size_t offset, size_t row, size_t col) const
    {
      P_ERROR_X1(offset + getHeader(offset)->col_off(col) + sizeof(double)*row < size, "Matrix element is out of bound.");
      return ((double*)((char*)address + offset + getHeader(offset)->col_off(col)))[row];
    }
    double& elem_im(size_t offset, size_t row, size_t col)
    {
      P_ERROR_X1(offset + getHeader(offset)->col_off(col) + sizeof(double)*row < size, "Matrix element is out of bound.");
      return ((double*)((char*)address + offset + getHeader(offset)->col_off_im(col)))[row];
    }
    const double& elem_im(size_t offset, size_t row, size_t col) const
    {
      P_ERROR_X1(offset + getHeader(offset)->col_off(col) + sizeof(double)*row < size, "Matrix element is out of bound.");
      return ((double*)((char*)address + offset + getHeader(offset)->col_off_im(col)))[row];
    }
    uint32_t getRows(size_t offset) const
    {
      return getHeader(offset)->mrows;
    }
    uint32_t getCols(size_t offset) const
    {
      return getHeader(offset)->ncols;
    }

    void mmapFileWrite(const std::string& fileName, size_t size);
    void mmapFileRead(const std::string& fileName);
//#ifdef _WIN32
//    HANDLE file;
//    OVERLAPPED fileOverlapped;
//    HANDLE mapHandle;
//    size_t filesize;
//#else
//    int    file;
//    struct flock fileLock;
//#endif

    file_mapping matFileMapping;
    mapped_region matFileMemory;
    
    static std::mutex fileLock;
    
    const std::string matFileName;
    void  *address;
    size_t size;
    size_t ncols;
    const bool wperm;
    size_t ndim;
    size_t npar;
    size_t nint;
    size_t ndeg;
    size_t nmul;

    size_t nint1;
    size_t nint2;
    size_t ndeg1;
    size_t ndeg2;

    bool   torus;

    size_t npoints_offset;  //T
    header npoints_header;

    size_t par_offset;      //T
    header par_header;
    
    size_t parnames_offset;
    header parnames_header;

    size_t mul_offset;
    header mul_header;

    size_t ntrivmul_offset;
    header ntrivmul_header;
    
    size_t magic_offset;
    header magic_header;

    size_t ndim_offset;     //T
    header ndim_header;

    size_t elem_offset;
    header elem_header;

    size_t mesh_offset;
    header mesh_header;

    size_t prof_offset;
    header prof_header;

    // for the torus
    size_t nint1_offset;     //T
    header nint1_header;

    size_t nint2_offset;     //T
    header nint2_header;

    size_t ndeg1_offset;     //T
    header ndeg1_header;

    size_t ndeg2_offset;     //T
    header ndeg2_header;

    size_t mesh1_offset;     //T
    header mesh1_header;

    size_t mesh2_offset;     //T
    header mesh2_header;

    size_t blanket_offset;   //T
    header blanket_header;
};

#endif
