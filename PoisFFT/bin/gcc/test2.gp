set xrange [0:1]
set yrange [0:2]
set zrange [0:3]
set xlabel "x"
set ylabel "y"
set zlabel "z"
splot 'test.txt' u 1:2:3:4 palette pt 9
pause -1
