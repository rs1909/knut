#-------------------------------------------------------------------------------
# UMFPACK makefile (for GNU make or original make)
#-------------------------------------------------------------------------------

# UMFPACK requires the AMD package to be in ../AMD

# compile everything (except hb, fortran, and fortran64), including AMD and the
# MATLAB mexFunctions
all:
	( cd ../AMD ; make )
	( cd Source ; make )
	- cat Doc/License

# compile just the C-callable libraries and demo programs (not mexFunctions)
lib:
	( cd ../AMD ; make lib )
	( cd Source ; make )
	- cat Doc/License

# compile the FORTRAN interface and demo program
fortran:
	( cd Demo   ; make fortran )

# compile the 64-bit FORTRAN interface and demo program
fortran64:
	( cd Demo   ; make fortran64 )

# compile the Harwell/Boeing demo program
hb:
	( cd Demo   ; make hb )

# remove object files, but keep the compiled programs and library archives
clean:
	( cd ../AMD ; make clean )
	( cd Source ; make clean )

# clean, and then remove compiled programs and library archives
purge:
	( cd ../AMD ; make purge )
	( cd Source ; make purge )

# create PDF documents for the original distribution
doc:

# get ready for distribution
dist: purge
	( cd ../AMD ; make dist )

DSDIRS  = Doc Include Lib Make Source
DSFILES = Makefile README

distdir:
	list='$(DSDIRS)'; for subdir in $$list; do \
	mkdir $(distdir)/$$subdir; \
	( cd $$subdir; make distdir="$(distdir)/$$subdir" distdir ) \
	done;
	list='$(DSFILES)'; for ff in $$list; do \
	cp -pR $$ff $(distdir); \
	done;
	