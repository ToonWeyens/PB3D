pdflatex flowchart.tex
pdftocairo -svg flowchart.pdf flowchart.svg
pdftoppm -png flowchart.pdf > flowchart.png
rm flowchart.pdf
