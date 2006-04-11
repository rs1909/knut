NAME = logistic

exdatadir = $(pkgdatadir)/demo/$(NAME)
dist_exdata_DATA = sys-$(NAME).cpp
nodist_exdata_DATA = sys-$(NAME).so

clean:
	-test -z sys-$(NAME).so || rm -f sys-$(NAME).so

sys-$(NAME).so: sys-$(NAME).cpp
	$(CXX) @MOD_CXXFLAGS@ -I../../src sys-$(NAME).cpp -o sys-$(NAME).so
