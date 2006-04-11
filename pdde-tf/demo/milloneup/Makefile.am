NAME = milloneup

exdatadir = $(pkgdatadir)/demo/$(NAME)
dist_exdata_DATA = sys-$(NAME).cpp
nodist_exdata_DATA = sys-$(NAME).so
dist_exdata_DATA += cfile-start cfile-ns cfile-tr

sys-$(NAME).so: sys-$(NAME).cpp
	$(CXX) @MOD_CXXFLAGS@ -I../../src sys-$(NAME).cpp -o sys-$(NAME).so

clean:
	-test -z sys-$(NAME).so || rm -f sys-$(NAME).so

test: $(NAME).so
	pdde -c cfile-start -o out-start -b br-start
	pdde -c cfile-ns -i out-start -o out-ns -b br-ns
test-tr:
	pdde -c cfile-tr -i out-start -o out-tr -b br-tr
