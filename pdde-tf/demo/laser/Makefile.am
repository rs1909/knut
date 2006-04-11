NAME = laser

exdatadir = $(pkgdatadir)/demo/$(NAME)
dist_exdata_DATA = sys-$(NAME).cpp
nodist_exdata_DATA = sys-$(NAME).so
dist_exdata_DATA += cfile-start cfile-ns cfile-pd cfile-lp in-start

sys-$(NAME).so: sys-$(NAME).cpp
	$(CXX) @MOD_CXXFLAGS@ -I../../src sys-$(NAME).cpp -o sys-$(NAME).so

clean:
	-test -z sys-$(NAME).so || rm -f sys-$(NAME).so

test: $(NAME).so
	pdde -c cfile-start -i in-start -o out-start -b br-start
test-ns:
	pdde -c cfile-ns -i out-start -o out-ns -b branch-ns
test-pd:
	pdde -c cfile-pd -i out-start -o out-pd -b branch-pd
test-lp:
	pdde -c cfile-lp -i out-start -o out-lp -b branch-lp
