set xrange [0:1]
set yrange [0:2]
set xlabel "x"
set ylabel "y"
set zlabel "V"
set ticslevel 0
set dgrid3d 30, 30
set hidden3d
set term x11 0
set title "RHO"
splot "outfileRHO.txt" u 1:2:3 with lines
set term x11 1
set title "SOR"
splot "outfileSOR.txt" u 1:2:3 with lines

#sp "outfileSOR.txt" u 1:2:3:4 w pm3d
pause -1
