NAME = traffic

exdatadir = $(pkgdatadir)/demo/$(NAME)
dist_exdata_DATA = sys-$(NAME).cpp sys-$(NAME).so
dist_exdata_DATA += cfile-start-a cfile-start-b \
	            cfile-ns-a cfile-ns-b \
	            cfile-nssw-a cfile-nssw-b \
	            cfile-lp-a cfile-lp-b

clean:
	-test -z sys-$(NAME).so || rm -f sys-$(NAME).so

sys-$(NAME).so: sys-$(NAME).cpp
	$(SHELL) ../../bin/pcompile sys-$(NAME).cpp

test2: out-start-a
	pdde -c cfile-nssw-a -i out-start-a -o out-nssw-a -b br-nssw-a

test3: out-start-b
	pdde -c cfile-nssw-b -i out-start-b -o out-nssw-b -b br-nssw-b

test4: out-nssw-a
	pdde -c cfile-lp-a -i out-nssw-a -o out-lp-a -b br-lp-a

test5: out-nssw-b
	pdde -c cfile-lp-b -i out-nssw-b -o out-lp-b -b br-lp-b
