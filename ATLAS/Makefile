DSDIRS = include lib
DSFILES = Makefile

distdir:
	list='$(DSDIRS)'; for subdir in $$list; do \
	mkdir $(distdir)/$$subdir; \
	( cd $$subdir; make distdir="$(distdir)/$$subdir" distdir ) \
	done;
	list='$(DSFILES)'; for ff in $$list; do \
	cp -pR $$ff $(distdir); \
	done;

all:
	
clean:
	
install:
	

