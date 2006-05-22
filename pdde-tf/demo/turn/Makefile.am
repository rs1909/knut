NAME = turn

exdatadir = $(pkgdatadir)/demo/$(NAME)
dist_exdata_DATA = sys-$(NAME).cpp
nodist_exdata_DATA = sys-$(NAME).so
dist_exdata_DATA += cfile-start cfile-ns cfile-nssw

sys-$(NAME).so: sys-$(NAME).cpp
	$(CXX) @MOD_CXXFLAGS@ -I../../src sys-$(NAME).cpp -o sys-$(NAME).so

clean:
	-test -z sys-$(NAME).so || rm -f sys-$(NAME).so

test: sys-$(NAME).so
	pdde -c cfile-start -o out-start -b br-start
	pdde -c cfile-ns -i out-start -o out-ns -b br-ns
	pdde -c cfile-nssw -i out-start -o out-nssw -b br-nssw
                
