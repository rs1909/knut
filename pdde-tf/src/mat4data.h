#ifndef MAT4DATA_H
#define MAT4DATA_H

#include <string>
#include "matrix.h"
#ifdef WIN32
#include <windows.h>
#endif

class mat4Data
{
	public:
		
		struct header
		{
			long type;
			long mrows;
			long ncols;
			long imagf;
			long namelen;
		
			const char* name( const void* offset ) const { return (char*)offset + namelen; }
			int size() const { return sizeof(struct header) + namelen; }
			int col_off( int i ) const { return size() + i*mrows*sizeof(double); }
			int col_off_im( int i ) const { return size() + (mrows*ncols + i*mrows)*sizeof(double); }
		};

		// opens the file, determines its size, maps the memory from file, sets up variables
		mat4Data( const std::string& fileName, int steps_, int ndim_, int npar_, int nint_, int ndeg_, int nmul_ );
		mat4Data( const std::string& fileName );
		// unmaps the memory, truncates the file if necessary, closes the file
		~mat4Data( );

		int  findMatrix( const char* name, mat4Data::header* found );
		void openReadOnly( const std::string& fileName );
		
		void   setPar( int n, const Vector& par );
		void   setMul( int n, const Vector& real, const Vector& imag );
		void   setMesh( int n, const Vector& mesh );
		void   setProfile( int n, const Vector& profile );
		void   getPar( int n, Vector& par ) const;
		double getPar( int n, int j ) const;
		void   getMul( int n, Vector& real,Vector& imag ) const;
		double getMulRe( int n, int j ) const;
		double getMulIm( int n, int j ) const;
		void   getMesh( int n, Vector& mesh ) const;
		double getMesh( int n, int j ) const;
		void   getProfile( int n, Vector& profile ) const;
// 		void   getProfile( int n, int i, Vector& profile ) const;
		double getProfile( int n, int d, int j ) const;
		
		int  getNDim() const { return ndim; }
		int  getNInt() const { return nint; }
		int  getNDeg() const { return ndeg; }
		int  getNPar() const { return npar; }
		int  getNMul() const { return nmul; }
		int  getMeshLength() const { return ndeg*nint+1; }
		int  getNCols() const { return ncols; }
		int  getNPoints() const { return static_cast<int>(((double*)((char*)address + npoints_offset + npoints_header.col_off(0)))[0]); }
		int  getNextBifurcation( int n, int aut ) const;
		int  getBifurcationType( int n, int aut ) const;
		
		struct header *getHeader( int offset ) { return (struct header*)((char*)address + offset); }
		const struct header *getHeader( int offset ) const { return (struct header*)((char*)address + offset); }
		double& elem( int offset, int row, int col )
		{
			return ((double*)((char*)address + offset + getHeader(offset)->col_off(col)))[row];
		}
		const double& elem( int offset, int row, int col ) const
		{
			return ((double*)((char*)address + offset + getHeader(offset)->col_off(col)))[row];
		}
		double& elem_im( int offset, int row, int col )
		{
			return ((double*)((char*)address + offset + getHeader(offset)->col_off_im(col)))[row];
		}
		const double& elem_im( int offset, int row, int col ) const
		{
			return ((double*)((char*)address + offset + getHeader(offset)->col_off_im(col)))[row];
		}

	private:
		inline int countUnstable( int n, int aut ) const;
		inline void findTrivialIndices( int n, int aut, int *imin, double* dmin ) const;

#ifdef WIN32
    HANDLE file;
    HANDLE mapHandle;
#else
	int    file;
#endif
	int    filesize;
	void  *address;
	int    size;
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
	
	int    npoints_offset;
	header npoints_header;
	
	int    par_offset;
	header par_header;

	int    mul_offset;
	header mul_header;

	int    ndim_offset;
	header ndim_header;

	int    nint_offset;
	header nint_header;

	int    ndeg_offset;
	header ndeg_header;
	
	int    mesh_offset;
	header mesh_header;

	int    prof_offset;
	header prof_header;
};

#endif
