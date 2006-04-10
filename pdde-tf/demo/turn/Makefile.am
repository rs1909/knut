NAME = turn

exdatadir = $(pkgdatadir)/demo/$(NAME)
dist_exdata_DATA = sys-$(NAME).cpp sys-$(NAME).so
dist_exdata_DATA += cfile-start cfile-ns cfile-nssw

sys-$(NAME).so: sys-$(NAME).cpp
	$(SHELL) ../../bin/pcompile sys-$(NAME).cpp

test: $(NAME).so
	pdde -c cfile-start -o out-start -b br-start
	pdde -c cfile-ns -i out-start -o out-ns -b br-ns
	pdde -c cfile-nssw -i out-start -o out-nssw -b br-nssw
                