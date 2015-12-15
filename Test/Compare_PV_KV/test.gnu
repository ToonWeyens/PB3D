set grid; set border 4095 front linetype -1 linewidth 1.000; set terminal pdf size 4,2; set output 'PV_KV_comparison.pdf';
set logscale y
plot 'test.dat' u 1:2 with linespoints lw 2 title 'X*PX delta', 'test.dat' u 1:3 with linespoints lw 2 title 'X*KX delta', 'test.dat' u 1:($2/$3) with linespoints lw 3 title 'X*PX/X*KX'
