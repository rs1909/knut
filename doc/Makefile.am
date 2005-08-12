SUBDIRS = fig
EXTRA_DIST = macros.tex ddedoc.cls manual-bibl.bib manual.tex

doc: macros.tex ddedoc.cls manual-bibl.bib manual.tex
	latex manual
	bibtex manual
	latex manual
	latex manual
	