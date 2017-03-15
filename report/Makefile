PROJECT = report.pdf
TEMP    = tmp

all: $(PROJECT)

$(TEMP):
	-rm -rf $(TEMP)
	mkdir $(TEMP)

clean:
	rm -rf $(TEMP)/* $(PROJECT)

%.pdf: $(TEMP) %.tex %.bib
	pdflatex -output-directory $(TEMP) $*.tex
	biber --output_directory $(TEMP) $*.bcf
	pdflatex -output-directory $(TEMP) $*.tex
	pdflatex -output-directory $(TEMP) $*.tex
	cp $(TEMP)/$@ $@

.PHONY: all clean
