# Build the report: make  (needs pdflatex and bibtex)
NHFOAM_report.pdf: *.tex 8-References.bib ../figures/*.pdf
	pdflatex -interaction=nonstopmode 0-main > /dev/null
	bibtex 0-main > /dev/null
	pdflatex -interaction=nonstopmode 0-main > /dev/null
	pdflatex -interaction=nonstopmode 0-main > /dev/null
	cp 0-main.pdf NHFOAM_report.pdf

clean:
	rm -f *.aux *.log *.bbl *.blg *.toc *.out 0-main.pdf
