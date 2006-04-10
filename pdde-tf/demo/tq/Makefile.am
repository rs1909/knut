NAME = tq

exdatadir = $(pkgdatadir)/demo/$(NAME)
dist_exdata_DATA = sys-$(NAME).cpp sys-$(NAME).so
dist_exdata_DATA += cfile-lp-1-a cfile-lp-2-a  cfile-ns      cfile-pdsw-2  cfile-tr \
		    cfile-lp-1-b cfile-lp-2-b  cfile-pdsw-1  cfile-start

sys-$(NAME).so: sys-$(NAME).cpp
	$(SHELL) ../../bin/pcompile sys-$(NAME).cpp

test:
	pdde -c cfile-start -o out-start -b br-start
test-tr:
	pdde -c cfile-tr -i out-start -o out-tr -b br-tr
test-pdsw:
	pdde -c cfile-pdsw-1 -i out-start -o out-pdsw-1 -b br-pdsw-1
	pdde -c cfile-pdsw-2 -i out-start -o out-pdsw-2 -b br-pdsw-2
test-lp-1-a:
	pdde -c cfile-lp-1-a -i out-pdsw-1 -o out-lp-1-a -b br-lp-1-a
test-lp-1-b:
	pdde -c cfile-lp-1-b -i out-pdsw-1 -o out-lp-1-b -b br-lp-1-b
test-lp-2-a:
	pdde -c cfile-lp-2-a -i out-pdsw-2 -o out-lp-2-a -b br-lp-2-a
test-lp-2-b:
	pdde -c cfile-lp-2-b -i out-pdsw-2 -o out-lp-2-b -b br-lp-2-b
test-ns:
	pdde -c cfile-ns -i out-start -o out-ns -b br-ns
