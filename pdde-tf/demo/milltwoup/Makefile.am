NAME = milltwoup

exdatadir = $(pkgdatadir)/demo/$(NAME)
dist_exdata_DATA = sys-$(NAME).cpp sys-$(NAME).so

sys-$(NAME).so: sys-$(NAME).cpp
	$(SHELL) ../../bin/pcompile sys-$(NAME).cpp
