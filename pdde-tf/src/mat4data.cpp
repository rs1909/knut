#include <iostream>
#include <string>

#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <sys/mman.h>

#include "matrix.h"
#include "pderror.h"
#include "mat4data.h"

int mat4Data::findMatrix( const char* name, mat4Data::header* found, int *sz )
{
	struct header hd;
	int cur_off = 0;
	int cur_size;
	do{
		memcpy( &hd, (char*)address + cur_off, sizeof(struct header) );
		if( hd.type != 0 )
		  { std::cout<<"not a double matrix\n"; PDError(-1); }
		if( hd.imagf == 0 )
		  cur_size = sizeof(struct header) + hd.namelen*sizeof(char) + hd.mrows*hd.ncols*sizeof(double);
		else
		  cur_size = sizeof(struct header) + hd.namelen*sizeof(char) + 2*hd.mrows*hd.ncols*sizeof(double);
		if( strncmp( name, (char*)address + cur_off + sizeof(struct header), 20 ) == 0 )
		{
			memcpy( found, &hd, sizeof(struct header) );
			*sz = cur_size;
			return cur_off;
		}
		cur_off += cur_size;
	}while( cur_off < size );
	return -1;
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
	
	par_offset = 0;
	strncpy( par_name, "pdde_par", 20 );
	par_header.type = 0;
	par_header.mrows = npar;
	par_header.ncols = ncols;
	par_header.imagf = 0;
	par_header.namelen = ((strlen(par_name)+sizeof(struct header)+1)/sizeof(double) + 1)*sizeof(double) - sizeof(struct header);
	par_size = sizeof(struct header) + par_header.namelen + ncols * npar * sizeof(double);
	
	mul_offset = par_offset + par_size;
	strncpy( mul_name, "pdde_mul", 20 );
	mul_header.type = 0;
	mul_header.mrows = nmul;
	mul_header.ncols = ncols;
	mul_header.imagf = 1;
	mul_header.namelen = ((strlen(mul_name)+sizeof(struct header)+1)/sizeof(double) + 1)*sizeof(double) - sizeof(struct header);
	mul_size = sizeof(struct header) + mul_header.namelen + ncols * 2 * nmul * sizeof(double);
	
	ndim_offset = mul_offset + mul_size;
	strncpy( ndim_name, "pdde_ndim", 20 );
	ndim_header.type = 0;
	ndim_header.mrows = 1;
	ndim_header.ncols = ncols;
	ndim_header.imagf = 0;
	ndim_header.namelen = ((strlen(ndim_name)+sizeof(struct header)+1)/sizeof(double) + 1)*sizeof(double) - sizeof(struct header);
	ndim_size = sizeof(struct header) + ndim_header.namelen + ncols * 1 * sizeof(double);

	nint_offset = ndim_offset + ndim_size;
	strncpy( nint_name, "pdde_nint", 20 );
	nint_header.type = 0;
	nint_header.mrows = 1;
	nint_header.ncols = ncols;
	nint_header.imagf = 0;
	nint_header.namelen = ((strlen(nint_name)+sizeof(struct header)+1)/sizeof(double) + 1)*sizeof(double) - sizeof(struct header);
	nint_size = sizeof(struct header) + nint_header.namelen + ncols * 1 * sizeof(double);

	ndeg_offset = nint_offset + nint_size;
	strncpy( ndeg_name, "pdde_ndeg", 20 );
	ndeg_header.type = 0;
	ndeg_header.mrows = 1;
	ndeg_header.ncols = ncols;
	ndeg_header.imagf = 0;
	ndeg_header.namelen = ((strlen(ndeg_name)+sizeof(struct header)+1)/sizeof(double) + 1)*sizeof(double) - sizeof(struct header);
	ndeg_size = sizeof(struct header) + ndeg_header.namelen + ncols * 1 * sizeof(double);

	mesh_offset = ndeg_offset + ndeg_size;
	strncpy( mesh_name, "pdde_mesh", 20 );
	mesh_header.type = 0;
	mesh_header.mrows = ndeg*nint+1;
	mesh_header.ncols = ncols;
	mesh_header.imagf = 0;
	mesh_header.namelen = ((strlen(mesh_name)+sizeof(struct header)+1)/sizeof(double) + 1)*sizeof(double) - sizeof(struct header);
	mesh_size = sizeof(struct header) + mesh_header.namelen + ncols * (ndeg*nint+1) * sizeof(double);

	prof_offset = mesh_offset + mesh_size;
	strncpy( prof_name, "pdde_prof", 20 );
	prof_header.type = 0;
	prof_header.mrows = ndim*(ndeg*nint+1);
	prof_header.ncols = ncols;
	prof_header.imagf = 0;
	prof_header.namelen = ((strlen(prof_name)+sizeof(struct header)+1)/sizeof(double) + 1)*sizeof(double) - sizeof(struct header);
	prof_size = sizeof(struct header) + prof_header.namelen + ncols * ndim*(ndeg*nint+1) * sizeof(double);
	size = prof_offset + prof_size;
	
	// creating the mmapped file
	if( ( file = open( fileName.c_str(), O_CREAT | O_RDWR | O_TRUNC, S_IRUSR | S_IWUSR ) ) == -1 )
	{ perror("mmappedPointData::mmappedPointData: unable to open file\n"); throw(-1); }
	
	if( lseek( file, size-1, SEEK_SET ) == -1 )
	{ perror("mmappedPointData::mmappedPointData: unable to seek file\n"); throw(-1); }
	
	if( write( file, "\0", 1 ) == -1 )
	{ perror("mmappedPointData::mmappedPointData: unable to write file\n"); throw(-1); }
	
	if( ( address = mmap( 0, size, PROT_WRITE, MAP_SHARED, file, 0 ) ) == MAP_FAILED )
	{ perror("mmappedPointData::mmappedPointData: unable to mmap file\n"); throw(-1); }

	*((struct header*)( (char*)address + par_offset) ) = par_header;
	*((struct header*)( (char*)address + mul_offset) ) = mul_header;
	*((struct header*)( (char*)address + ndim_offset) ) = ndim_header;
	*((struct header*)( (char*)address + nint_offset) ) = nint_header;
	*((struct header*)( (char*)address + ndeg_offset) ) = ndeg_header;
	*((struct header*)( (char*)address + mesh_offset) ) = mesh_header;
	*((struct header*)( (char*)address + prof_offset) ) = prof_header;
	strncpy( (char*)address + par_offset + sizeof(struct header), par_name, 20 );
	strncpy( (char*)address + mul_offset + sizeof(struct header), mul_name, 20 );
	strncpy( (char*)address + ndim_offset + sizeof(struct header), ndim_name, 20 );
	strncpy( (char*)address + nint_offset + sizeof(struct header), nint_name, 20 );
	strncpy( (char*)address + ndeg_offset + sizeof(struct header), ndeg_name, 20 );
	strncpy( (char*)address + mesh_offset + sizeof(struct header), mesh_name, 20 );
	strncpy( (char*)address + prof_offset + sizeof(struct header), prof_name, 20 );

// 	std::cout<<"NDIM "<<ndim<<" NINT "<<nint<<" NDEG "<<ndeg<<" NPAR "<<npar<<" NMUL "<<nmul<<"\n";
}

mat4Data::mat4Data( const std::string& fileName )
{
	wperm = false;
	this->openReadOnly( fileName );
}

void mat4Data::openReadOnly( const std::string& fileName )
{
	if( ( file = open( fileName.c_str(), O_RDONLY ) ) == -1 )
	{ perror("mmappedPointData::mmappedPointData: unable to open file\n"); throw(-1); }
	
	struct stat filestat;
	if( fstat( file, &filestat ) != 0 )
	{ perror("mmappedPointData::mmappedPointData: unable to stat file\n"); throw(-1); }
	filesize = filestat.st_size;
	size = filesize;
	
	if( ( address = mmap( 0, filesize, PROT_READ, MAP_PRIVATE, file, 0 ) ) == MAP_FAILED )
	{ perror("mmappedPointData::mmappedPointData: unable to mmap file\n"); throw(-1); }
	
	if( (par_offset = findMatrix( "pdde_par", &par_header, &par_size )) == -1 ) { std::cout<<"err1"; PDError(-1); }
	npar = par_header.mrows; std::cout<<"NPAR "<<npar<<"\n";
	ncols = par_header.ncols;
	if( par_header.imagf != 0 ) { std::cout<<"err2"; PDError(-1); }

	if( (mul_offset = findMatrix( "pdde_mul", &mul_header, &mul_size )) == -1 ) { std::cout<<"err3"; PDError(-1); }
	nmul = mul_header.mrows; std::cout<<"NMUL "<<nmul<<"\n";
	if( mul_header.ncols != ncols ) { std::cout<<"err4"; PDError(-1); }
	if( mul_header.imagf == 0 ) { std::cout<<"err5"; PDError(-1); }
	
	if( (ndim_offset = findMatrix( "pdde_ndim", &ndim_header, &ndim_size )) == -1 ) { std::cout<<"err6"; PDError(-1); }
	if( ndim_header.mrows != 1 ) { std::cout<<"err7 "<<ndim_header.mrows<<"\n"; PDError(-1); }
	if( ndim_header.ncols != ncols ) { std::cout<<"err9"; PDError(-1); }
	if( ndim_header.imagf != 0 ) { std::cout<<"err9"; PDError(-1); }
	ndim = *((double*)((char*)address + ndim_offset + ndim_header.col_off(0)));
// 	std::cout<<"NDIM "<<ndim<<"\n";
	
	if( (nint_offset = findMatrix( "pdde_nint", &nint_header, &nint_size )) == -1 ) { std::cout<<"err10"; PDError(-1); }
	if( nint_header.mrows != 1 ) { std::cout<<"err11"; PDError(-1); }
	if( nint_header.ncols != ncols ) { std::cout<<"err12"; PDError(-1); }
	if( nint_header.imagf != 0 ) { std::cout<<"err13"; PDError(-1); }
	nint = *((double*)((char*)address + nint_offset + nint_header.col_off(0)));
// 	std::cout<<"NINT "<<nint<<"\n";
	
	if( (ndeg_offset = findMatrix( "pdde_ndeg", &ndeg_header, &ndeg_size )) == -1 ) { std::cout<<"err14"; PDError(-1); }
	if( ndeg_header.mrows != 1 ) { std::cout<<"err15"; PDError(-1); }
	if( ndeg_header.ncols != ncols ) { std::cout<<"err16"; PDError(-1); }
	if( ndeg_header.imagf != 0 ) { std::cout<<"err17"; PDError(-1); }
	ndeg = *((double*)((char*)address + ndeg_offset + ndeg_header.col_off(0)));
// 	std::cout<<"NDEG "<<ndeg<<"\n";

	if( (mesh_offset = findMatrix( "pdde_mesh", &mesh_header, &mesh_size )) == -1 ) { std::cout<<"err18"; PDError(-1); }
	if( mesh_header.mrows != ndeg*nint+1 ) { std::cout<<"err19 "<<mesh_header.mrows<<"\n"; PDError(-1); }
	if( mesh_header.ncols != ncols ) { std::cout<<"err20"; PDError(-1); }
	if( mesh_header.imagf != 0 ) { std::cout<<"err21"; PDError(-1); }
	
	if( (prof_offset = findMatrix( "pdde_prof", &prof_header, &prof_size )) == -1 ) { std::cout<<"err22"; PDError(-1); }
	if( prof_header.mrows != ndim*(ndeg*nint+1) ) { std::cout<<"err23"; PDError(-1); }
	if( prof_header.ncols != ncols ) { std::cout<<"err24"; PDError(-1); }
	if( prof_header.imagf != 0 ) { std::cout<<"err25"; PDError(-1); }
	
// 	std::cout<<"NDIM "<<ndim<<" NINT "<<nint<<" NDEG "<<ndeg<<" NPAR "<<npar<<" NMUL "<<nmul<<"\n";
}

mat4Data::~mat4Data()
{
	if( munmap( address, size ) != 0 )
	{ perror("mmappedPointData::~mmappedPointData: unable to munmap file\n"); throw(-1); }
	if( close( file ) != 0 )
	{ perror("mmappedPointData::~mmappedPointData: unable to close file\n"); throw(-1); }
}

void mat4Data::setPar( int n, const Vector& par )
{
	if( wperm && n < ncols )
	{
		if( par.Size() <= npar )
		{
			for( int i = 0; i < par.Size(); ++i)
				((double*)( (char*)address + par_offset + par_header.col_off(n) ))[i] = par(i);
			for( int i = par.Size(); i < npar; ++i)
				((double*)( (char*)address + par_offset + par_header.col_off(n) ))[i] = 0.0;
		}else
		{
			std::cout<<"setPar 1";
			PDError(-1);
		}
	}else
	{
		std::cout<<"setPar 2";
		PDError(-1);
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
				((double*)( (char*)address + mul_offset + mul_header.col_off(n) ))[i] = re(i);
				((double*)( (char*)address + mul_offset + mul_header.col_off_im(n) ))[i] = im(i);
			}
			for( int i = re.Size(); i < nmul; ++i )
			{
				((double*)( (char*)address + mul_offset + mul_header.col_off(n) ))[i] = 0.0;
				((double*)( (char*)address + mul_offset + mul_header.col_off_im(n) ))[i] = 0.0;
			}
		}else
		{
			std::cout<<"setMul 1";
			PDError(-1);
		}
	}else
	{
		std::cout<<"setMul 2";
		PDError(-1);
	}
}

void mat4Data::setMesh( int n, const Vector& mesh )
{
	if( wperm && n < ncols )
	{
		if( mesh.Size() == ndeg*nint+1 )
		{
			for( int i = 0; i < ndeg*nint+1; ++i)
				((double*)( (char*)address + mesh_offset + mesh_header.col_off(n) ))[i] = mesh(i);
		}else
		{
			std::cout<<"setMesh 1";
			PDError(-1);
		}
	}else
	{
		std::cout<<"setMesh 2";
		PDError(-1);
	}
}

void mat4Data::setProfile( int n, const Vector& prof )
{
	if( wperm && n < ncols )
	{
		if( prof.Size() == ndim*(ndeg*nint+1) )
		{
			((double*)((char*)address + ndim_offset + ndim_header.col_off(n)))[0] = ndim;
			((double*)((char*)address + nint_offset + nint_header.col_off(n)))[0] = nint;
			((double*)((char*)address + ndeg_offset + ndeg_header.col_off(n)))[0] = ndeg;
			for( int i = 0; i < ndim*(ndeg*nint+1); ++i)
				((double*)( (char*)address + prof_offset + prof_header.col_off(n) ))[i] = prof(i);
		}else
		{
			std::cout<<"setProf 1";
			PDError(-1);
		}
	}else
	{
		std::cout<<"setProf 2";
		PDError(-1);
	}
}

void mat4Data::getPar( int n, Vector& par )
{
	if( n < ncols )
	{
		if( par.Size() <= npar )
		{
			for( int i = 0; i < par.Size(); ++i)
				par(i) = ((double*)( (char*)address + par_offset + par_header.col_off(n) ))[i];
		}else
		{
			std::cout<<"getPar 1";
			PDError(-1);
		}
	}else
	{
		std::cout<<"getPar 2";
		PDError(-1);
	}
}

void mat4Data::getMul( int n, Vector& re, Vector& im )
{
	if( n < ncols && re.Size() == im.Size() )
	{
		if( re.Size() <= nmul )
		{
			for( int i = 0; i < re.Size(); ++i )
			{
				re(i) = ((double*)( (char*)address + mul_offset + mul_header.col_off(n) ))[i];
				im(i) = ((double*)( (char*)address + mul_offset + mul_header.col_off_im(n) ))[i];
			}
		}else
		{
			std::cout<<"getMul 1";
			PDError(-1);
		}
	}else
	{
		std::cout<<"getMul 2";
		PDError(-1);
	}
}

void mat4Data::getMesh( int n, Vector& mesh )
{
	if( n < ncols )
	{
		if( mesh.Size() == ndeg*nint+1 )
		{
			for( int i = 0; i < ndeg*nint+1; ++i)
				mesh(i) = ((double*)( (char*)address + mesh_offset + mesh_header.col_off(n) ))[i];
		}else
		{
			std::cout<<"getMesh 1";
			PDError(-1);
		}
	}else
	{
		std::cout<<"getMesh 2";
		PDError(-1);
	}
}

void mat4Data::getProfile( int n, Vector& prof )
{
	if( n < ncols )
	{
		if( prof.Size() == ndim*(ndeg*nint+1) )
		{
			for( int i = 0; i < ndim*(ndeg*nint+1); ++i)
				prof(i) = ((double*)( (char*)address + prof_offset + prof_header.col_off(n) ))[i];
		}else
		{
			std::cout<<"getProf 1";
			PDError(-1);
		}
	}else
	{
		std::cout<<"getProf 2";
		PDError(-1);
	}
}
