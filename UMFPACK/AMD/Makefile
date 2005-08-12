#-------------------------------------------------------------------------------
# AMD Makefile (for GNU Make or original make)
#-------------------------------------------------------------------------------

# Compile everything, including the C-callable routine and the mexFunctions.
# Do not compile the FORTRAN versions.
all:
	( cd Source ; make )

# compile just the C-callable libraries and demo programs (not mexFunctions)
lib:
	( cd Source ; make )

# compile the FORTRAN libraries and demo programs (not compiled by "make all")
fortran:
	( cd Source ; make fortran )

# compile a FORTRAN demo program that calls the C version of AMD
# (not compiled by "make all")
cross:
	( cd Demo   ; make cross )

# remove object files, but keep the compiled programs and library archives
clean:
	( cd Source ; make clean )

# clean, and then remove compiled programs and library archives
purge:
	( cd Source ; make purge )

doc:

# get ready for distribution
dist: purge

DSDIRS = Lib Source
DSFILES = Include Make  Makefile README

distdir:
	list='$(DSDIRS)'; for subdir in $$list; do \
	mkdir $(distdir)/$$subdir; \
	( cd $$subdir; make distdir="$(distdir)/$$subdir" distdir ) \
	done;
	cp -pR $(DSFILES) $(distdir)  
                                      