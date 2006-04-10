SUBDIRS = fig html
EXTRA_DIST = macros.tex ddedoc.cls manual-bibl.bib manual.tex
pkgdata_DATA = manual.pdf

doc: manual.dvi

manual.dvi: macros.tex ddedoc.cls manual-bibl.bib manual.tex
	latex manual
	bibtex manual
	latex manual
	latex manual
	
manual.pdf: manual.dvi
	dvipdfm manual.dvi
	