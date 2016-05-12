# Based on http://www.gnuplotting.org/interpolation-of-heat-maps/, http://gnuplot.sourceforge.net/demo/heatmaps.4.gnu and http://stackoverflow.com/questions/23505711/gnuplot-how-to-avoid-white-margin-in-pdf-through-adjusting-page-size-and-other
# general settings and variables
n = 10
unset key
set grid; set border 4095 front linetype -1 linewidth 1.000; set terminal pdf size 5,4; set output 'alpha_study.pdf';

# 2D plot
set lmargin at screen 0.07;
set rmargin at screen 0.80;
set tmargin at screen 0.97;
set bmargin at screen 0.07;
set xrange [ -0.000000 : n-1 ] nowriteback
set yrange [ -0.000000 : 10.00000 ] nowriteback
set cblabel "Eigenvalue" 
set pm3d map
set pm3d interpolate 0,0
splot "alpha_study.dat" matrix rowheaders columnheaders using 1:2:3 every :2

# 3D plot
set output 'alpha_study_3D.pdf'
unset view
unset lmargin
unset rmargin
unset tmargin
unset bmargin
set xrange [1:n]
set yrange [0:1]
set zrange [-0.2:0]
set ticslevel 0
set xlabel 'theta limit index'
set ylabel 'alpha [pi]'
set zlabel 'EV'
splot for [i=1:n] "alpha_study.dat" u (i):1:i+1 every ::1 with lines lw 2

# plot of variances
unset pm3d
set output 'alpha_study_sigma.pdf'
set xrange [1:n]
set yrange [-0.7:0]
set xlabel 'theta limit index'
set ylabel 'sigma / average'
plot 'alpha_study.dat' matrix using 1:3 every 1:n:1:14 with lines lw 2
