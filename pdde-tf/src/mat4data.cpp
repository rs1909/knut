#include <iostream>
#include <string>

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
	
	if( lseek( file, size-1, SEEK_SET ) == -1 )
	{ P_ERROR_X2( false, "mmappedPointData::mmappedPointData: unable to seek file\n", strerror( errno ) ); }
	
	if( write( file, "\0", 1 ) == -1 )
	{ P_ERROR_X2( false, "mmappedPointData::mmappedPointData: unable to write file\n", strerror( errno ) ); }

	void *address;
	if( ( address = mmap( 0, size, PROT_WRITE, MAP_SHARED, file, 0 ) ) == MAP_FAILED )
	{ P_ERROR_X2( false, "mmappedPointData::mmappedPointData: unable to mmap file\n", strerror( errno ) ); }
	return address;
}

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
    
    if( SetFilePointer( file, size-1, NULL, FILE_BEGIN ) == 0 )
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
	if( size > approxSize ) P_ERROR_X2(false,"BAd size approxiamion",size-approxSize);
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

	if( (mul_offset = findMatrix( "pdde_mul", &mul_header )) == -1 ) P_MESSAGE("err3");
	nmul = mul_header.mrows;
	if( mul_header.ncols != ncols ) P_MESSAGE("err4");
	if( mul_header.imagf == 0 ) P_MESSAGE("err5");
	
	if( (ndim_offset = findMatrix( "pdde_ndim", &ndim_header )) == -1 ) P_MESSAGE("err6");
	if( ndim_header.mrows != 1 ) P_MESSAGE("err7 ");
	if( ndim_header.ncols != ncols ) P_MESSAGE("err9");
	if( ndim_header.imagf != 0 ) P_MESSAGE("err9");
	ndim = static_cast<int>(*((double*)((char*)address + ndim_offset + ndim_header.col_off(0))));
	
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
}

mat4Data::~mat4Data()
{
#ifndef WIN32
	if( munmap( address, size ) != 0 )
	{ P_ERROR_X2( false, "mmappedPointData::~mmappedPointData: unable to munmap file\n", strerror( errno ) ); }
	if( close( file ) != 0 )
	{ P_ERROR_X2( false, "mmappedPointData::~mmappedPointData: unable to close file\n", strerror( errno ) ); }
#else
     if( address != 0 )
     {
         UnmapViewOfFile(address);
         CloseHandle(file);
         CloseHandle(mapHandle);
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
