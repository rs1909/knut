all:
	( cd UMFPACK; make )

DSDIRS = AMD UMFPACK
DSFILES = Makefile

distdir:
	list='$(DSDIRS)'; for subdir in $$list; do \
	mkdir $(distdir)/$$subdir; \
	( cd $$subdir; make distdir="$(distdir)/$$subdir" distdir ) \
	done;
	cp -pR $(DSFILES) $(distdir)
	   
clean:
	( cd UMFPACK; make clean )

install:
	

