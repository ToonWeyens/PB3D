gprof PB3D gmon.out > gmon.txt
cat gmon.txt | ~/Scripts/gprof2dot/jrfonseca.gprof2dot/gprof2dot.py -n 0.01 -e 0.01 > profile.xdot
rm gmon.txt
xdot profile.xdot
rm profile.xdot
