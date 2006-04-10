NAME = milloneup

exdatadir = $(pkgdatadir)/demo/$(NAME)
dist_exdata_DATA = sys-$(NAME).cpp sys-$(NAME).so
dist_exdata_DATA += cfile-start cfile-ns cfile-tr

sys-$(NAME).so: sys-$(NAME).cpp
	$(SHELL) ../../bin/pcompile sys-$(NAME).cpp

clean:
	-test -z sys-$(NAME).so || rm -f sys-$(NAME).so

test: $(NAME).so
	pdde -c cfile-start -o out-start -b br-start
	pdde -c cfile-ns -i out-start -o out-ns -b br-ns
test-tr:
	pdde -c cfile-tr -i out-start -o out-tr -b br-tr
