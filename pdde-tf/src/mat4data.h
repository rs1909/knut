#ifndef MAT4DATA_H
#define MAT4DATA_H

#include <string>
#include "matrix.h"

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

		int  findMatrix( const char* name, mat4Data::header* found, int* sz );
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

	private:

	int    file;
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
	
	int    par_offset;
	char   par_name[20];
	header par_header;
	int    par_size;

	int    mul_offset;
	char   mul_name[20];
	header mul_header;
	int    mul_size;

	int    ndim_offset;
	char   ndim_name[20];
	header ndim_header;
	int    ndim_size;

	int    nint_offset;
	char   nint_name[20];
	header nint_header;
	int    nint_size;

	int    ndeg_offset;
	char   ndeg_name[20];
	header ndeg_header;
	int    ndeg_size;
	
	int    mesh_offset;
	char   mesh_name[20];
	header mesh_header;
	int    mesh_size;

	int    prof_offset;
	char   prof_name[20];
	header prof_header;
	int    prof_size;
};

#endif
