include Makefile.rules

default: pdf

FNAME=imprints

SUBDIRS=

STAT_PICS_SRC= figs/static/storage_rpp8 figs/static/storage_rpp16 figs/static/storage_rpp32 figs/static/storage_rpp64 figs/static/storagetime_rpp8 figs/static/storagetime_rpp16 figs/static/storagetime_rpp32 figs/static/storagetime_rpp64 figs/static/percentage figs/static/entropy figs/static/cumul figs/static/qtime figs/static/improv figs/static/stats figs/static/cum_queries
STAT_PICS_PDF=$(STAT_PICS_SRC:%=%.pdf)
STAT_PICS_PS=$(STAT_PICS_SRC:%=%.eps)
STAT_PICS_MPS=$(STAT_PICS_SRC:%=%.mps)

DYN_PICS_SRC= figs/example1 figs/cfpstruct figs/example2 figs/imprint figs/example3
DYN_PICS_PDF=$(DYN_PICS_SRC:%=%.pdf)
DYN_PICS_PS=$(DYN_PICS_SRC:%=%.eps)
DYN_PICS_MPS=$(DYN_PICS_SRC:%=%.mps)

PICS_PDF=$(DYN_PICS_PDF) $(STAT_PICS_PDF)
PICS_PS=$(DYN_PICS_PS) $(STAT_PICS_PS)
PICS_MPS=$(DYN_PICS_MPS) $(STAT_PICS_MPS)

TEX_SRCS = imprints.tex introduction.tex background.tex construction.tex operators.tex evaluation.tex conclusions.tex figs/static/snap1.tex figs/static/snap2.tex figs/static/snap3.tex figs/static/snap4.tex figs/static/snap5.tex updates.tex

vps: ps
	evince $(FNAME).ps &

vpdf: pdf
	evince $(FNAME).pdf &

all: ps pdf
pdf: $(FNAME).pdf
ps: $(FNAME).ps

$(FNAME).pdf : $(PICS_PS) $(TEX_SRCS) bibliography.bib
#	for i in $(SUBDIRS) ; do make -C $$i ; done
	for i in aux bbl blg log ; do rm -f $(FNAME).$$i ; done
	latex $(FNAME)
	bibtex $(FNAME) -terse
	latex $(FNAME)
#	perl -i -p -n0 -e 's/\n\n\\bib/\n\\vspace{-1.5mm}\n\\bib/g;' $(FNAME).bbl
#	perl -i -p -e 'if (/vspace/g && !$$i) {$$_="\n";$$i++}' $(FNAME).bbl 
	latex $(FNAME)
	dvips -t letter -Ppdf $(FNAME).dvi -o $(FNAME).ps
	ps2pdf $(FNAME).ps

####
####
#make ps to achieve the letter format and the ps2pdf!!!
####
$(FNAME).ps : $(PICS_PS) $(TEX_SRCS) bibliography.bib
#	for i in $(SUBDIRS) ; do make -C $$i ; done
	for i in aux bbl blg log ; do rm -f $(FNAME).$$i ; done
	latex $(FNAME)
	bibtex $(FNAME) -terse
	latex $(FNAME)
#	perl -i -p -n0 -e 's/\n\n\\bib/\n\\vspace{-1.5mm}\n\\bib/g;' $(FNAME).bbl
#	perl -i -p -e 'if (/vspace/g && !$$i) {$$_="\n";$$i++}' $(FNAME).bbl 
	latex $(FNAME)
	dvips -t letter -Ppdf $(FNAME).dvi -o $(FNAME).ps

spell : $(TEX_SRCS)
	for i in $(TEX_SRCS) ; do aspell -t -p $$PWD/aspell-ignore -c $$i ; done

clean:
#	for i in $(SUBDIRS) ; do make -C $$i clean ; done
	rm -f $(FNAME).[ablpd]* *.log \
	$(DYN_PICS_PDF) $(DYN_PICS_PS) $(DYN_PICS_MPS) *~ \
	$(ID).ps $(ID).pdf *.bak  \
	comment.cut \
	*/*.mpx */*.log */*.mps
