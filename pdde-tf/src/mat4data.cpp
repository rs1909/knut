#include <iostream>
#include <string>
#include <cmath>
#include <cfloat>

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

int mat4Data::findMatrix( const char* name, mat4Data::header* found )
{
	struct header hd;
	int cur_off = 0;
	int cur_size;
	do{
		memcpy( &hd, (char*)address + cur_off, sizeof(struct header) );
		P_ERROR_X( hd.type == 0, "not a double matrix");
		if( hd.imagf == 0 )
		  cur_size = sizeof(struct header) + hd.namelen*sizeof(char) + hd.mrows*hd.ncols*sizeof(double);
		else
		  cur_size = sizeof(struct header) + hd.namelen*sizeof(char) + 2*hd.mrows*hd.ncols*sizeof(double);
		if( strncmp( name, (char*)address + cur_off + sizeof(struct header), 20 ) == 0 )
		{
			memcpy( found, &hd, sizeof(struct header) );
			return cur_off;
		}
		cur_off += cur_size;
	}while( cur_off < size );
	return -1;
}

#ifndef WIN32

static inline void *mmapFileWrite( int& file, const std::string& fileName, int size )
{
	if( ( file = open( fileName.c_str(), O_CREAT | O_RDWR | O_TRUNC, S_IRUSR | S_IWUSR ) ) == -1 )
	{ P_ERROR_X2( false, "mmappedPointData::mmappedPointData: unable to open file\n", strerror( errno ) ); }
	
	if( ftruncate( file, size ) != 0 )
	{ P_ERROR_X2( false, "mmappedPointData::mmappedPointData: unable to truncate file\n", strerror( errno ) ); }
	
	void *address;
	if( ( address = mmap( 0, size, PROT_WRITE, MAP_SHARED, file, 0 ) ) == MAP_FAILED )
	{ P_ERROR_X2( false, "mmappedPointData::mmappedPointData: unable to mmap file\n", strerror( errno ) ); }
	return address;
}

// static inline void shrinkMap( void* address, int old_size, int new_size )
// {
// 	P_ERROR_X( old_size <= new_size, "mmappedPointData::mmappedPointData: larger new_size" );
// 	void* new_address;
// 	if( ( new_address = mremap( address, old_size, new_size, 0 ) ) != address )
// 	{ P_ERROR_X2( false, "mmappedPointData::mmappedPointData: unable remap file\n", strerror( errno ) ); }
// }

static inline void *mmapFileRead( int& file, const std::string& fileName, int& size )
{
	if( ( file = open( fileName.c_str(), O_RDONLY ) ) == -1 )
	{ P_ERROR_X2( false, "mmappedPointData::mmappedPointData: unable to open file\n", strerror( errno ) ); }
	
	struct stat filestat;
	if( fstat( file, &filestat ) != 0 )
	{ P_ERROR_X2( false, "mmappedPointData::mmappedPointData: unable to stat file\n", strerror( errno ) ); }
	int filesize = filestat.st_size;
	size = filesize;
	
	void *address;
	if( ( address = mmap( 0, filesize, PROT_READ, MAP_PRIVATE, file, 0 ) ) == MAP_FAILED )
	{ P_ERROR_X2( false, "mmappedPointData::mmappedPointData: unable to mmap file\n", strerror( errno ) ); }
	return address;
}

#else

static inline void *mmapFileWrite( HANDLE& file, HANDLE& mapHandle, const std::string& fileName, int size )
{
	if( ( file = CreateFile( fileName.c_str(),
		FILE_WRITE_DATA | FILE_READ_DATA,
		FILE_SHARE_READ,
		NULL,
		CREATE_ALWAYS,
		FILE_ATTRIBUTE_NORMAL,
		NULL ) ) == NULL )
	{ P_ERROR_X2( false, "mmappedPointData::mmappedPointData: unable open file\n", static_cast<int>(GetLastError()) ); }

	if( SetFilePointer( file, size, NULL, FILE_BEGIN ) == 0 )
	{ P_ERROR_X2( false, "mmappedPointData::mmappedPointData: unable to seek\n", static_cast<int>(GetLastError()) ); }
	P_ERROR_X2( SetEndOfFile( file ), "mmappedPointData::mmappedPointData: unable to SetEndOfFile\n", static_cast<int>(GetLastError()) );

	if( (mapHandle = CreateFileMapping( file, NULL, PAGE_READWRITE, 0, size, fileName.c_str() ) ) == 0 )
	{ P_ERROR_X2( false, "mmappedPointData::mmappedPointData: unable to map file\n", static_cast<int>(GetLastError()) ); }
	
	void *address = MapViewOfFile( mapHandle, FILE_MAP_WRITE, 0, 0, 0 );
	if( address != NULL ) return address;
	else P_ERROR_X2( false, "mmappedPointData::mmappedPointData: unable to map file\n", static_cast<int>(GetLastError()) );
	return 0;
}

static inline void *mmapFileRead( HANDLE& file, HANDLE& mapHandle, const std::string& fileName, int& size )
{
	if( ( file = CreateFile( fileName.c_str(), 
		FILE_READ_DATA,
		FILE_SHARE_READ, 
		NULL, 
		OPEN_EXISTING, 
		FILE_ATTRIBUTE_NORMAL, 
		NULL ) ) == NULL )
	{ P_ERROR_X2( false, "mmappedPointData::mmappedPointData: unable open file\n", static_cast<int>(GetLastError()) ); }

	size = GetFileSize( file, NULL );
	if( (mapHandle = CreateFileMapping( file, NULL, PAGE_READONLY, 0, 0, fileName.c_str() ) ) == 0 )
	{ P_ERROR_X2( false, "mmappedPointData::mmappedPointData: unable to map file\n", static_cast<int>(GetLastError()) ); }
	
	void *address = MapViewOfFile( mapHandle, FILE_MAP_READ, 0, 0, 0 );
	if( address != NULL ) return address;
	else P_ERROR_X2( false, "mmappedPointData::mmappedPointData: unable to map file\n", static_cast<int>(GetLastError()) );
	return 0;
}

#endif

// returns the size of the whole data
static inline int createMatrixHeader( void* address, int offset, mat4Data::header* hd, const char* name, int rows, int cols )
{
	strncpy( (char*)address + offset + sizeof(mat4Data::header), name, 20 );
	hd->type = 0;
	hd->mrows = rows;
	hd->ncols = cols;
	hd->imagf = 0;
	hd->namelen = ((strlen(name)+sizeof(mat4Data::header)+1)/sizeof(double) + 1)*sizeof(double) - sizeof(mat4Data::header);
	*((mat4Data::header*)( (char*)address + offset) ) = *hd;
	return sizeof(mat4Data::header) + hd->namelen + rows * cols * sizeof(double);
}

// returns the size of the whole data
static inline int createComplexMatrixHeader( void* address, int offset, mat4Data::header* hd, const char* name, int rows, int cols )
{
	strncpy( (char*)address + offset + sizeof(mat4Data::header), name, 20 );
	hd->type = 0;
	hd->mrows = rows;
	hd->ncols = cols;
	hd->imagf = 1;
	hd->namelen = ((strlen(name)+sizeof(mat4Data::header)+1)/sizeof(double) + 1)*sizeof(double) - sizeof(mat4Data::header);
	*((mat4Data::header*)( (char*)address + offset) ) = *hd;
	return sizeof(mat4Data::header) + hd->namelen + 2 * rows * cols * sizeof(double);
}

mat4Data::mat4Data( const std::string& fileName, int steps_, int ndim_, int npar_, int nint_, int ndeg_, int nmul_ )
{
	wperm = true;
	torus = false;
	
	ncols = steps_;
	ndim = ndim_;
	npar = npar_;
	nint = nint_;
	ndeg = ndeg_;
	nmul = nmul_;
	// for tori
// 	nint1 = params.getNInt1();
// 	nint2 = params.getNInt2();
// 	ndeg1 = params.getNDeg1();
// 	ndeg2 = params.getNDeg2();

	const int approxSize = 8*(sizeof(header) + 20) + sizeof(double)*(1+ncols*(npar+2*nmul+1+1+1+(ndim+1)*(ndeg*nint+1)));
 #ifndef WIN32
	address = mmapFileWrite( file, fileName, approxSize );
 #else
    address = mmapFileWrite( file, mapHandle, fileName, approxSize );
 #endif	
	// creating the matrices
	npoints_offset = 0;
	int npoints_size = createMatrixHeader( address, npoints_offset, &npoints_header, "pdde_npoints", 1, 1 );

	par_offset = npoints_offset + npoints_size;
	int par_size = createMatrixHeader( address, par_offset, &par_header, "pdde_par", npar, ncols );
	
	mul_offset = par_offset + par_size;
	int mul_size = createComplexMatrixHeader( address, mul_offset, &mul_header, "pdde_mul", nmul, ncols );
	
	ndim_offset = mul_offset + mul_size;
	int ndim_size = createMatrixHeader( address, ndim_offset, &ndim_header, "pdde_ndim", 1, ncols );
	
	nint_offset = ndim_offset + ndim_size;
	int nint_size = createMatrixHeader( address, nint_offset, &nint_header, "pdde_nint", 1, ncols );
	
	ndeg_offset = nint_offset + nint_size;
	int ndeg_size = createMatrixHeader( address, ndeg_offset, &ndeg_header, "pdde_ndeg", 1, ncols );

	mesh_offset = ndeg_offset + ndeg_size;
	int mesh_size = createMatrixHeader( address, mesh_offset, &mesh_header, "pdde_mesh", ndeg*nint+1, ncols );
	
	prof_offset = mesh_offset + mesh_size;
	int prof_size = createMatrixHeader( address, prof_offset, &prof_header, "pdde_prof", ndim*(ndeg*nint+1), ncols );
	size = prof_offset + prof_size;
	if( size > approxSize ) P_ERROR_X2(false,"Bad size approxiamion",size-approxSize);
}

mat4Data::mat4Data( const std::string& fileName, int steps_, int ndim_, int npar_, int nint1_, int nint2_, int ndeg1_, int ndeg2_ )
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
	
	const int approxSize = 10*(sizeof(header) + 20) + sizeof(double)*(1+ncols*(npar+1+1+1+1+1+ndeg1*nint1+ndeg2*nint2+ndim*ndeg1*nint1*ndeg2*nint2));
 #ifndef WIN32
	address = mmapFileWrite( file, fileName, approxSize );
 #else
    address = mmapFileWrite( file, mapHandle, fileName, approxSize );
 #endif	
	// creating the matrices
	npoints_offset = 0;
	int npoints_size = createMatrixHeader( address, npoints_offset, &npoints_header, "pdde_npoints", 1, 1 );
	
	par_offset = npoints_offset + npoints_size;
	int par_size = createMatrixHeader( address, par_offset, &par_header, "pdde_par", npar, ncols );
	
	ndim_offset = par_offset + par_size;
	int ndim_size = createMatrixHeader( address, ndim_offset, &ndim_header, "pdde_ndim", 1, ncols );
	
	nint1_offset = ndim_offset + ndim_size;
	int nint1_size = createMatrixHeader( address, nint1_offset, &nint1_header, "pdde_nint1", 1, ncols );
	
	nint2_offset = nint1_offset + nint1_size;
	int nint2_size = createMatrixHeader( address, nint2_offset, &nint2_header, "pdde_nint2", 1, ncols );
	
	ndeg1_offset = nint2_offset + nint2_size;
	int ndeg1_size = createMatrixHeader( address, ndeg1_offset, &ndeg1_header, "pdde_ndeg1", 1, ncols );

	ndeg2_offset = ndeg1_offset + ndeg1_size;
	int ndeg2_size = createMatrixHeader( address, ndeg2_offset, &ndeg2_header, "pdde_ndeg2", 1, ncols );

	mesh1_offset = ndeg2_offset + ndeg2_size;
	int mesh1_size = createMatrixHeader( address, mesh1_offset, &mesh1_header, "pdde_mesh1", nint1*ndeg1, ncols );

	mesh2_offset = mesh1_offset + mesh1_size;
	int mesh2_size = createMatrixHeader( address, mesh2_offset, &mesh2_header, "pdde_mesh2", nint2*ndeg2, ncols );

	blanket_offset = mesh2_offset + mesh2_size;
	int blanket_size = createMatrixHeader( address, blanket_offset, &blanket_header, "pdde_blanket", ndim*nint1*ndeg1*nint2*ndeg2, ncols );
	size = blanket_offset + blanket_size;
	P_ERROR_X2( size <= approxSize, "Bad size approxiamion",size-approxSize);
}

mat4Data::mat4Data( const std::string& fileName )
{
	wperm = false;
	this->openReadOnly( fileName );
}

void mat4Data::openReadOnly( const std::string& fileName )
{
 #ifndef WIN32
    address = mmapFileRead( file, fileName, size );
 #else
    address = mmapFileRead( file, mapHandle, fileName, size );
 #endif
	
	if( (npoints_offset = findMatrix( "pdde_npoints", &npoints_header )) == -1 ) P_MESSAGE("err-4");
	P_ERROR_X( npoints_header.mrows == 1, "err-3" );
	P_ERROR_X( npoints_header.ncols == 1, "err-2" );
	P_ERROR_X( npoints_header.imagf == 0, "err-1" );

	if( (par_offset = findMatrix( "pdde_par", &par_header )) == -1 ) P_MESSAGE("err1");
	npar = par_header.mrows;
	ncols = par_header.ncols;
	if( par_header.imagf != 0 ) P_MESSAGE("err2");

	if( (ndim_offset = findMatrix( "pdde_ndim", &ndim_header )) == -1 ) P_MESSAGE("err6");
	if( ndim_header.mrows != 1 ) P_MESSAGE("err7 ");
	if( ndim_header.ncols != ncols ) P_MESSAGE("err9");
	if( ndim_header.imagf != 0 ) P_MESSAGE("err9");
	ndim = static_cast<int>(*((double*)((char*)address + ndim_offset + ndim_header.col_off(0))));

	if( (mul_offset = findMatrix( "pdde_mul", &mul_header )) != -1 )
	{
		torus = false;
		// periodic solutions
		nmul = mul_header.mrows;
		if( mul_header.ncols != ncols ) P_MESSAGE("err4");
		if( mul_header.imagf == 0 ) P_MESSAGE("err5");
		
		if( (nint_offset = findMatrix( "pdde_nint", &nint_header )) == -1 ) P_MESSAGE("err10");
		if( nint_header.mrows != 1 ) P_MESSAGE("err11");
		if( nint_header.ncols != ncols ) P_MESSAGE("err12");
		if( nint_header.imagf != 0 ) P_MESSAGE("err13");
		nint = static_cast<int>(*((double*)((char*)address + nint_offset + nint_header.col_off(0))));
		
		if( (ndeg_offset = findMatrix( "pdde_ndeg", &ndeg_header )) == -1 ) P_MESSAGE("err14");
		if( ndeg_header.mrows != 1 ) P_MESSAGE("err15");
		if( ndeg_header.ncols != ncols ) P_MESSAGE("err16");
		if( ndeg_header.imagf != 0 ) P_MESSAGE("err17");
		ndeg = static_cast<int>(*((double*)((char*)address + ndeg_offset + ndeg_header.col_off(0))));
	
		if( (mesh_offset = findMatrix( "pdde_mesh", &mesh_header )) == -1 ) P_MESSAGE("err18");
		if( mesh_header.mrows != ndeg*nint+1 ) P_MESSAGE("err19");
		if( mesh_header.ncols != ncols ) P_MESSAGE("err20");
		if( mesh_header.imagf != 0 ) P_MESSAGE("err21");
		
		if( (prof_offset = findMatrix( "pdde_prof", &prof_header )) == -1 ) P_MESSAGE("err22");
		if( prof_header.mrows != ndim*(ndeg*nint+1) ) P_MESSAGE("err23");
		if( prof_header.ncols != ncols ) P_MESSAGE("err24");
		if( prof_header.imagf != 0 ) P_MESSAGE("err25");
	}else
	{
		torus = true;
		// quasiperiodic solutions
		if( (nint1_offset = findMatrix( "pdde_nint1", &nint1_header )) == -1 ) P_MESSAGE("err26");
		if( nint1_header.mrows != 1 ) P_MESSAGE("err27");
		if( nint1_header.ncols != ncols ) P_MESSAGE("err28");
		if( nint1_header.imagf != 0 ) P_MESSAGE("err29");
		nint1 = static_cast<int>(*((double*)((char*)address + nint1_offset + nint1_header.col_off(0))));
		
		if( (nint2_offset = findMatrix( "pdde_nint2", &nint2_header )) == -1 ) P_MESSAGE("err30");
		if( nint2_header.mrows != 1 ) P_MESSAGE("err31");
		if( nint2_header.ncols != ncols ) P_MESSAGE("err32");
		if( nint2_header.imagf != 0 ) P_MESSAGE("err33");
		nint2 = static_cast<int>(*((double*)((char*)address + nint2_offset + nint2_header.col_off(0))));
		
		if( (ndeg1_offset = findMatrix( "pdde_ndeg2", &ndeg1_header )) == -1 ) P_MESSAGE("err34");
		if( ndeg1_header.mrows != 1 ) P_MESSAGE("err35");
		if( ndeg1_header.ncols != ncols ) P_MESSAGE("err36");
		if( ndeg1_header.imagf != 0 ) P_MESSAGE("err37");
		ndeg1 = static_cast<int>(*((double*)((char*)address + ndeg1_offset + ndeg1_header.col_off(0))));
		
		if( (ndeg2_offset = findMatrix( "pdde_ndeg2", &ndeg2_header )) == -1 ) P_MESSAGE("err38");
		if( ndeg2_header.mrows != 1 ) P_MESSAGE("err39");
		if( ndeg2_header.ncols != ncols ) P_MESSAGE("err40");
		if( ndeg2_header.imagf != 0 ) P_MESSAGE("err41");
		ndeg2 = static_cast<int>(*((double*)((char*)address + ndeg2_offset + ndeg2_header.col_off(0))));

		if( (mesh1_offset = findMatrix( "pdde_mesh1", &mesh1_header )) == -1 ) P_MESSAGE("err42");
		if( mesh1_header.mrows != nint1*ndeg1 ) P_MESSAGE("err43");
		if( mesh1_header.ncols != ncols ) P_MESSAGE("err44");
		if( mesh1_header.imagf != 0 ) P_MESSAGE("err45");
		
		if( (mesh2_offset = findMatrix( "pdde_mesh2", &mesh2_header )) == -1 ) P_MESSAGE("err46");
		if( mesh2_header.mrows != nint2*ndeg2 ) P_MESSAGE("err47");
		if( mesh2_header.ncols != ncols ) P_MESSAGE("err48");
		if( mesh2_header.imagf != 0 ) P_MESSAGE("err49");
		
		if( (blanket_offset = findMatrix( "pdde_blanket", &blanket_header )) == -1 ) P_MESSAGE("err50");
		if( blanket_header.mrows != ndim*nint1*ndeg1*nint2*ndeg2 ) P_MESSAGE("err51");
		if( blanket_header.ncols != ncols ) P_MESSAGE("err52");
		if( blanket_header.imagf != 0 ) P_MESSAGE("err53");
	}
}

mat4Data::~mat4Data()
{
#ifndef WIN32
	if( munmap( address, size ) != 0 )
	{ P_ERROR_X2( false, "mmappedPointData::~mmappedPointData: unable to munmap file\n", strerror( errno ) ); }
	if( wperm )
	{
		if( ftruncate( file, size) != 0 )
		{ P_ERROR_X2( false, "mmappedPointData::mmappedPointData: unable to truncate file\n", strerror( errno ) ); }
	}
	if( close( file ) != 0 )
	{ P_ERROR_X2( false, "mmappedPointData::~mmappedPointData: unable to close file\n", strerror( errno ) ); }
#else
	if( address != 0 )
	{
		UnmapViewOfFile(address);
		CloseHandle(mapHandle);
		if( wperm )
		{
			if( SetFilePointer( file, size, NULL, FILE_BEGIN ) == 0 )
			{ P_ERROR_X2( false, "mmappedPointData::mmappedPointData: unable to seek\n", static_cast<int>(GetLastError()) ); }
			P_ERROR_X2( SetEndOfFile( file ), "mmappedPointData::mmappedPointData: unable to SetEndOfFile\n", static_cast<int>(GetLastError()) );
		}
		CloseHandle(file);
	}
#endif
}

void mat4Data::setPar( int n, const Vector& par )
{
	if( wperm && n < ncols )
	{
		if( par.Size() <= npar )
		{
			for( int i = 0; i < par.Size(); ++i)
				elem( par_offset, i, n ) = par(i);
			for( int i = par.Size(); i < npar; ++i)
				elem( par_offset, i, n ) = 0.0;
		}else
		{
			P_MESSAGE("setPar 1");
		}
	}else
	{
		P_MESSAGE("setPar 2");
	}
}

void mat4Data::setMul( int n, const Vector& re, const Vector& im )
{
	if( wperm && n < ncols && re.Size() == im.Size() )
	{
		if( re.Size() <= nmul )
		{
			for( int i = 0; i < re.Size(); ++i )
			{
				elem( mul_offset, i, n ) = re(i);
				elem_im( mul_offset, i, n ) = im(i);
			}
			for( int i = re.Size(); i < nmul; ++i )
			{
				elem( mul_offset, i, n ) = 0.0;
				elem_im( mul_offset, i, n ) = 0.0;
			}
		}else
		{
			P_MESSAGE("setMul 1");
		}
	}else
	{
		P_MESSAGE("setMul 2");
	}
}

void mat4Data::setMesh( int n, const Vector& mesh )
{
	if( wperm && n < ncols )
	{
		if( mesh.Size() == ndeg*nint+1 )
		{
			for( int i = 0; i < ndeg*nint+1; ++i)
				elem( mesh_offset, i, n ) = mesh(i);
		}else
		{
			P_MESSAGE("setMesh 1");
		}
	}else
	{
		P_MESSAGE("setMesh 2");
	}
}

void mat4Data::setProfile( int n, const Vector& prof )
{
	if( wperm && n < ncols )
	{
		if( prof.Size() == ndim*(ndeg*nint+1) )
		{
			const int curr_npoints = static_cast<int>( elem( npoints_offset, 0, 0 ) );
			if( n > curr_npoints ) elem(npoints_offset,0,0) = n;
			elem( ndim_offset, 0, n ) = ndim;
			elem( nint_offset, 0, n ) = nint;
			elem( ndeg_offset, 0, n ) = ndeg;
			for( int i = 0; i < ndim*(ndeg*nint+1); ++i ) elem( prof_offset, i, n ) = prof(i);
		}else
		{
			P_MESSAGE("setProf 1");
		}
	}else
	{
		P_MESSAGE("setProf 2");
	}
}

// write in order size() = [ ndeg1*nint1, ndeg2*nint2, ndim ]
void mat4Data::setBlanket( int n, const Vector& blanket )
{
	if( wperm && n < ncols )
	{
		if( blanket.Size() == ndim*(ndeg1*nint1*ndeg2*nint2) )
		{
			const int curr_npoints = static_cast<int>( elem( npoints_offset, 0, 0 ) );
			if( n > curr_npoints ) elem(npoints_offset,0,0) = n;
			elem( ndim_offset, 0, n ) = ndim;
			elem( nint1_offset, 0, n ) = nint1;
			elem( nint2_offset, 0, n ) = nint2;
			elem( ndeg1_offset, 0, n ) = ndeg1;
			elem( ndeg2_offset, 0, n ) = ndeg2;
// 			for( int i1 = 0; i1 < nint1; ++i1 ) {
// 				for( int i2 = 0; i2 < nint2; ++i2 ) {
// 					for( int j1 = 0; j1 < ndeg1; ++j1 ) {
// 						for( int j2 = 0; j2 < ndeg2; ++j2 ) {
// 							for( int p = 0; p < ndim; ++p ) {
// 								elem( blanket_offset, j1 + ndeg1*(i1 + nint1*(j2 + ndeg2*(i2 + nint2*p))), n )
// 									= blanket(p + ndim*(j1 + ndeg1*(j2 + ndeg2*(i1 + nint1*i2))));
// 							}
// 						}
// 					}
// 				}
// 			}
			for( int i = 0; i < ndim*(ndeg1*nint1*ndeg2*nint2); ++i ) elem( blanket_offset, i, n ) = blanket(i);
		}else
		{
			P_MESSAGE("setBlanket 1");
		}
	}else
	{
		P_MESSAGE("setBlanket 2");
	}
}

void mat4Data::getPar( int n, Vector& par ) const
{
	if( n < ncols )
	{
		if( par.Size() <= npar )
		{
			for( int i = 0; i < par.Size(); ++i)
				par(i) = elem( par_offset, i, n );
		}else
		{
			P_MESSAGE("getPar 1");
		}
	}else
	{
		P_MESSAGE("getPar 2");
	}
}

double mat4Data::getPar( int n, int j ) const
{
	return elem( par_offset, j, n );
}

void mat4Data::getMul( int n, Vector& re, Vector& im ) const
{
	if( n < ncols && re.Size() == im.Size() )
	{
		if( re.Size() <= nmul )
		{
			for( int i = 0; i < re.Size(); ++i )
			{
				re(i) = elem( mul_offset, i, n );
				im(i) = elem_im( mul_offset, i, n );
			}
		}else
		{
			P_MESSAGE("getMul 1");
		}
	}else
	{
		P_MESSAGE("getMul 2");
	}
}

double mat4Data::getMulRe( int n, int j ) const
{
	return elem( mul_offset, j, n );;
}

double mat4Data::getMulIm( int n, int j ) const
{
	return elem_im( mul_offset, j, n );
}

void mat4Data::getMesh( int n, Vector& mesh ) const
{
	if( n < ncols )
	{
		if( mesh.Size() == ndeg*nint+1 )
		{
			for( int i = 0; i < ndeg*nint+1; ++i)
				mesh(i) = elem( mesh_offset, i, n );
		}else
		{
			P_MESSAGE("getMesh 1");
		}
	}else
	{
		P_MESSAGE("getMesh 2");
	}
}

double mat4Data::getMesh( int n, int j ) const
{
	return elem( mesh_offset, j, n );
}

void mat4Data::getProfile( int n, Vector& prof ) const
{
	if( n < ncols )
	{
		if( prof.Size() == ndim*(ndeg*nint+1) )
		{
			for( int i = 0; i < ndim*(ndeg*nint+1); ++i)
				prof(i) = elem( prof_offset, i, n );
		}else
		{
			P_MESSAGE("getProf 1");
		}
	}else
	{
		P_MESSAGE("getProf 2");
	}
}

double mat4Data::getProfile( int n, int d, int j ) const
{
	return elem( prof_offset, j, n );
}

inline void mat4Data::findTrivialIndices( int n, int aut, int *imin, double* dmin ) const
{
	for( int i = 0; i < aut; ++i )
	{
		imin[i] = -1;
		dmin[i] = DBL_MAX;
	}
	for( int i = 0; i < getNMul(); ++i )
	{
		const double mabs = fabs(sqrt((getMulRe(n,i)-1.0)*(getMulRe(n,i)-1.0) + getMulIm(n,i)*getMulIm(n,i)));
		for( int j = 0; j < aut; ++j )
		{
			if( dmin[j] > mabs )
			{
				dmin[j] = mabs;
				imin[j] = i;
				break;
			}
		}
	}
}

inline int mat4Data::countUnstable( int n, int aut ) const
{
	int imin[aut];
	double dmin[aut];
	findTrivialIndices( n, aut, imin, dmin );
	int ustab = 0;
	for( int i = 0; i < getNMul(); ++i )
	{
		const double mabs = (getMulRe(n,i)*getMulRe(n,i) + getMulIm(n,i)*getMulIm(n,i));
		bool ok = true;
		for( int j = 0; j < aut; ++j ) if( i == imin[j] ) ok = false;
		if( ok && mabs >= 1.0 )  ++ustab;
	}
	return ustab;
}

int mat4Data::getNextBifurcation( int n, int aut ) const
{
	int p_ustab = countUnstable( n, aut );
	for( int i = n+1; i < getNPoints(); ++i )
	{
		int ustab = countUnstable( i, aut );
		if( ustab != p_ustab ) return i;
		p_ustab = ustab;
	}
	return -1;
}

int mat4Data::getBifurcationType( int n, int aut ) const
{
	int imin[aut];
	double dmin[aut];
	findTrivialIndices( n, aut, imin, dmin );
	double dminLP = DBL_MAX, dminPD = DBL_MAX, dminNS = DBL_MAX;
	int iminLP = -1, iminPD = -1, iminNS = -1;
	for( int i = 0; i < getNMul(); ++i )
	{
		const double mre = getMulRe(n,i);
		const double mim = getMulIm(n,i);
		bool ok = true;
		for( int j = 0; j < aut; ++j ) if( i == imin[j] ) ok = false;
		if( ok )
		{
			const double LPabs = fabs(sqrt((mre-1.0)*(mre-1.0)+mim*mim));
			const double PDabs = fabs(sqrt((mre+1.0)*(mre+1.0)+mim*mim));
			const double NSabs = fabs(sqrt(mre*mre+mim*mim)-1.0);
			if( (dminLP > LPabs) && (mim == 0.0) ) { dminLP = LPabs; iminLP = i; }
			if( (dminPD > PDabs) && (mim == 0.0) ) { dminPD = PDabs; iminPD = i; }
			if( (dminNS > NSabs) && (mim != 0.0) ) { dminNS = NSabs; iminNS = i; }
		}
	}
	if( (dminLP < dminPD) && (dminLP < dminNS) ) return 1; // BifTFLP;
	else if( (dminPD < dminLP) && (dminPD < dminNS) ) return 2; // BifTFPD;
	else if( (dminNS < dminPD) && (dminNS < dminLP) ) return 3; // BifTFNS;
	else return 0; //SolTF;
}
