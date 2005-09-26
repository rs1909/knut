#-------------------------------------------------------------------------------
# UMFPACK makefile (for GNU make or original make)
#-------------------------------------------------------------------------------

# UMFPACK requires the AMD package to be in ../AMD

# compile everything (except hb, fortran, and fortran64), including AMD and the
# MATLAB mexFunctions
all:
	( cd ../AMD ; $(MAKE) )
	( cd Source ; $(MAKE) )
	- cat Doc/License

# compile just the C-callable libraries and demo programs (not mexFunctions)
lib:
	( cd ../AMD ; $(MAKE) lib )
	( cd Source ; $(MAKE) )
	- cat Doc/License

# compile the FORTRAN interface and demo program
fortran:
	( cd Demo   ; $(MAKE) fortran )

# compile the 64-bit FORTRAN interface and demo program
fortran64:
	( cd Demo   ; $(MAKE) fortran64 )

# compile the Harwell/Boeing demo program
hb:
	( cd Demo   ; $(MAKE) hb )

# remove object files, but keep the compiled programs and library archives
clean:
	( cd ../AMD ; $(MAKE) clean )
	( cd Source ; $(MAKE) clean )

# clean, and then remove compiled programs and library archives
purge:
	( cd ../AMD ; $(MAKE) purge )
	( cd Source ; $(MAKE) purge )

# create PDF documents for the original distribution
doc:

# get ready for distribution
dist: purge
	( cd ../AMD ; $(MAKE) dist )

DSDIRS  = Doc Include Lib Make Source
DSFILES = Makefile README

distdir:
	list='$(DSDIRS)'; for subdir in $$list; do \
	mkdir $(distdir)/$$subdir; \
	( cd $$subdir; $(MAKE) distdir="$(distdir)/$$subdir" distdir ) \
	done;
	list='$(DSFILES)'; for ff in $$list; do \
	cp -pR $$ff $(distdir); \
	done;
	