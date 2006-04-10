NAME = glass

exdatadir = $(pkgdatadir)/demo/$(NAME)
dist_exdata_DATA = sys-$(NAME).cpp sys-$(NAME).so
dist_exdata_DATA += cfile-start cfile-pdsw cfile-pd

sys-$(NAME).so: sys-$(NAME).cpp
	$(SHELL) ../../bin/pcompile sys-$(NAME).cpp

clean:
	-test -z sys-$(NAME).so || rm -f sys-$(NAME).so

test: sys-$(NAME).so
	pdde -c cfile-start -o in-start -b br-start
	pdde -c cfile-pdsw -i in-start -o in-pdsw -b br-pdsw
test-pd:
	pdde -c cfile-pd -i in-start -o in-pd -b br-pd
