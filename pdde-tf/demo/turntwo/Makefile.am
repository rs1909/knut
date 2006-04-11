NAME = turntwo

exdatadir = $(pkgdatadir)/demo/$(NAME)
dist_exdata_DATA = sys-$(NAME).cpp
nodist_exdata_DATA = sys-$(NAME).so
dist_exdata_DATA += cfile-start cfile-nssw cfile-ns-1 cfile-ns-2 cfile-lp-1 cfile-lp-2

sys-$(NAME).so: sys-$(NAME).cpp
	$(CXX) @MOD_CXXFLAGS@ -I../../src sys-$(NAME).cpp -o sys-$(NAME).so

clean:
	-test -z sys-$(NAME).so || rm -f sys-$(NAME).so

test:
	pdde -c cfile-start -o out-start -b br-start
test-nssw:
	pdde -c cfile-nssw -i out-start -o out-nssw -b br-nssw
test-ns:
	pdde -c cfile-ns-1 -i out-start -o out-ns-1 -b br-ns-1
	pdde -c cfile-ns-2 -i out-start -o out-ns-2 -b br-ns-2
test-lp:
	pdde -c cfile-lp-1 -i out-nssw -o out-lp-1 -b br-lp-1
	pdde -c cfile-lp-2 -i out-nssw -o out-lp-2 -b br-lp-2
