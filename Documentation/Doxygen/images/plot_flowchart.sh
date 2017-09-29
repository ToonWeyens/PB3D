lualatex flowchart.tex
ln -sf ../../Fonts/Lato* .
pdftocairo -svg flowchart.pdf flowchart.svg
pdftoppm -png flowchart.pdf > flowchart.png
rm Lato*
