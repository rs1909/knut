NAME = milltwofull

exdatadir = $(pkgdatadir)/demo/$(NAME)
dist_exdata_DATA = sys-$(NAME).cpp
nodist_exdata_DATA = sys-$(NAME).so

clean:
	-test -z sys-$(NAME).so || rm -f sys-$(NAME).so

sys-$(NAME).so: sys-$(NAME).cpp
	$(SHELL) ../../bin/pcompile sys-$(NAME).cpp
