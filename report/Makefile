PROJECT = report.pdf abstract.pdf report_se499.pdf
TEMP    = tmp

all: $(PROJECT)

$(TEMP):
	-rm -rf $(TEMP)
	mkdir $(TEMP)

clean:
	rm -rf $(TEMP)/* $(PROJECT)

%.pdf: $(TEMP) %.tex %.bib
	pdflatex -output-directory $(TEMP) $*.tex
	cd $(TEMP) && makeindex $*.nlo -s nomencl.ist -o $*.nls
	pdflatex -output-directory $(TEMP) $*.tex
	biber --output_directory $(TEMP) $*.bcf
	pdflatex -output-directory $(TEMP) $*.tex
	cp $(TEMP)/$@ $@

%.pdf: $(TEMP) %.tex
	pdflatex -output-directory $(TEMP) $*.tex
	pdflatex -output-directory $(TEMP) $*.tex
	cp $(TEMP)/$@ $@

.PHONY: all clean
