#-------------------------------------------------------------------------------
# AMD Makefile (for GNU Make or original make)
#-------------------------------------------------------------------------------

# Compile everything, including the C-callable routine and the mexFunctions.
# Do not compile the FORTRAN versions.
all:
	( cd Source ; $(MAKE) )

# compile just the C-callable libraries and demo programs (not mexFunctions)
lib:
	( cd Source ; $(MAKE) )

# compile the FORTRAN libraries and demo programs (not compiled by "$(MAKE) all")
fortran:
	( cd Source ; $(MAKE) fortran )

# compile a FORTRAN demo program that calls the C version of AMD
# (not compiled by "$(MAKE) all")
cross:
	( cd Demo   ; $(MAKE) cross )

# remove object files, but keep the compiled programs and library archives
clean:
	( cd Source ; $(MAKE) clean )

# clean, and then remove compiled programs and library archives
purge:
	( cd Source ; $(MAKE) purge )

doc:

# get ready for distribution
dist: purge

DSDIRS = Include Lib Make Source
DSFILES = Makefile README

distdir:
	list='$(DSDIRS)'; for subdir in $$list; do \
	mkdir $(distdir)/$$subdir; \
	( cd $$subdir; $(MAKE) distdir="$(distdir)/$$subdir" distdir ) \
	done;
	cp -pR $(DSFILES) $(distdir)  
                                      