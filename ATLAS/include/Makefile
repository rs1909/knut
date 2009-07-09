DSDIRS = 
DSFILES = cblas.h clapack.h Makefile

distdir:
	list='$(DSDIRS)'; for subdir in $$list; do \
	mkdir $(distdir)/$$subdir; \
	( cd $$subdir; $(MAKE) distdir="../$(distdir)/$$subdir" distdir ) \
	done;
	list='$(DSFILES)'; for ff in $$list; do \
	cp -pR $$ff $(distdir); \
	done;

all:
	
clean:
	
install:
	

