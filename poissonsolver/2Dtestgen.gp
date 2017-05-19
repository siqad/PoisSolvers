#Set grid and label parameters
set xrange [0:1]
set yrange [0:2]
set xlabel "x"
set ylabel "y"
set zlabel "V"
set ticslevel 0
set dgrid3d 30, 30
set hidden3d

#Save as .pngs
set terminal png
set output "rho.png"
set title "RHO - sin(2*pi*x/Lx)"
splot "outfileRHO.txt" u 1:2:3 with lines
set output "V.png"
set title "V (slice at z = Lz/2), obtained by SOR"
splot "outfileSOR.txt" u 1:2:3 with lines

#show both in gnuplot terminal simultaneously
set term x11 0
set title "RHO - sin(2*pi*x/Lx)"
splot "outfileRHO.txt" u 1:2:3 with lines
#the last one to be plotted is interactive
set term x11 1
set title "V (slice at z = Lz/2), obtained by SOR"
splot "outfileSOR_GEN.txt" u 1:2:3 with lines

pause -1
