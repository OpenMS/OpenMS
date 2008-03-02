set ylabel "ion count"
set xlabel "concentration"
set key left Left reverse
set terminal postscript eps 
set output "AdditiveSeries_1_gnuplot_tmp.eps"
plot "AdditiveSeries_1_gnuplot_tmp.dat"  w points ps 2 pt 1 lt 8 title "data" ,  0.0782404+0.210347*x lt 2 lw 3 title "linear regression: 0.0782404 + 0.210347 * x" , "AdditiveSeries_1_gnuplot_tmp.dat"  w points ps 2 pt 1 lt 8 notitle , "AdditiveSeries_1_gnuplot_tmp.err"  using ($1):(0) w points pt 13 ps 2 lt 1 title "x-intercept: -0.371959" , "AdditiveSeries_1_gnuplot_tmp.err"  w xerrorbars lw 3 lt 1 title "95% interval: [ -0.403779, -0.341014 ]"
