#Set grid and label parameters
set view map
set xrange [0:1]
set yrange [0:1]
set xlabel "x"
set ylabel "z"
set zlabel "V"
set ticslevel 0
set hidden3d
set palette model RGB
set palette defined
#Save as .pngs
set terminal png
set output "rho.png"
set title "RHO"
splot "outfileRHO.txt" u 1:2:3 with pm3d
set output "V.png"
set title "V (slice at y = Ly/2), obtained by SOR"
splot "outfileMG.txt" u 1:2:3 with pm3d
set output "eps.png"
set title "Permittivity"
splot "outfileEPS.txt" u 1:2:3 with pm3d

#Save as .svgs
set terminal svg
set output "rho.svg"
set title "RHO"
splot "outfileRHO.txt" u 1:2:3 with pm3d
set output "V.svg"
set title "V (slice at y = Ly/2), obtained by SOR"
splot "outfileMG.txt" u 1:2:3 with pm3d
set output "eps.svg"
set title "Permittivity"
splot "outfileEPS.txt" u 1:2:3 with pm3d

#show both in gnuplot terminal simultaneously
set term x11 0
set title "Permittivity"
splot "outfileEPS.txt" u 1:2:3 with pm3d
set term x11 1
set title "RHO"
splot "outfileRHO.txt" u 1:2:3 with pm3d
#the last one to be plotted is interactive
set term x11 2
set title "V (slice at y = Ly/2), obtained by SOR"
splot "outfileMG.txt" u 1:2:3 with pm3d

pause -1
