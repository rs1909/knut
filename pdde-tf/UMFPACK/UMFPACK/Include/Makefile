DSDIRS = 
DSFILES = Makefile \
	umfpack_col_to_triplet.h   umfpack_numeric.h          umfpack_save_numeric.h \
	umfpack_defaults.h         umfpack_qsymbolic.h        umfpack_save_symbolic.h \
	umfpack_free_numeric.h     umfpack_report_control.h   umfpack_scale.h \
	umfpack_free_symbolic.h    umfpack_report_info.h      umfpack_solve.h \
	umfpack_get_determinant.h  umfpack_report_matrix.h    umfpack_symbolic.h \
	umfpack_get_lunz.h         umfpack_report_numeric.h   umfpack_tictoc.h \
	umfpack_get_numeric.h      umfpack_report_perm.h      umfpack_timer.h \
	umfpack_get_symbolic.h     umfpack_report_status.h    umfpack_transpose.h \
	umfpack.h                  umfpack_report_symbolic.h  umfpack_triplet_to_col.h \
	umfpack_load_numeric.h     umfpack_report_triplet.h   umfpack_wsolve.h \
	umfpack_load_symbolic.h    umfpack_report_vector.h    umppack_global.h

distdir:
	list='$(DSDIRS)'; for subdir in $$list; do \
	mkdir $(distdir)/$$subdir; \
	( cd $$subdir; $(MAKE) distdir="$(distdir)/$$subdir" distdir ) \
	done;
	list='$(DSFILES)'; for ff in $$list; do \
	cp -pR $$ff $(distdir); \
	done;

all:
	
clean:
	
install:
	

