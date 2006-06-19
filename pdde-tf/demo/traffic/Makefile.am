NAME = traffic3cars

exdatadir = $(pkgdatadir)/demo/$(NAME)
dist_exdata_DATA = sys-$(NAME).cpp
nodist_exdata_DATA = sys-$(NAME).so
dist_exdata_DATA += cfile-start-a cfile-start-b \
	            cfile-ns-a cfile-ns-b \
	            cfile-nssw-a cfile-nssw-b \
	            cfile-lp-a cfile-lp-b \
	            README sys-traffic50cars.cpp

clean:
	-test -z sys-$(NAME).so || rm -f sys-$(NAME).so

sys-$(NAME).so: sys-$(NAME).cpp
	$(CXX) @MOD_CXXFLAGS@ -I../../src sys-$(NAME).cpp -o sys-$(NAME).so

test1a:
	pdde -c cfile-start-a -o out-start-a -b br-start-a

test1b:
	pdde -c cfile-start-b -o out-start-b -b br-start-b

test2a: out-start-a
	pdde -c cfile-nssw-a -i out-start-a -o out-nssw-a -b br-nssw-a

test2b: out-start-b
	pdde -c cfile-nssw-b -i out-start-b -o out-nssw-b -b br-nssw-b

test3a: out-nssw-a
	pdde -c cfile-lp-a -i out-nssw-a -o out-lp-a -b br-lp-a

test3b: out-nssw-b
	pdde -c cfile-lp-b -i out-nssw-b -o out-lp-b -b br-lp-b
