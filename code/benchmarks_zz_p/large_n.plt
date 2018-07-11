set xlabel "n" offset 14,1.5
set ylabel "time in sec."
set terminal pdf
set key left top

set out outfile
set style line 1 lt 1 lw 8 lc 4
set style line 2 lt 1 lw 4 lc 7
set style line 3 lt 0 lw 8 lc 8
set xtics 0,1000,5000

set multiplot 
set size ratio 1.5
set origin -0.45,0
set yrange [0:4]
plot infile1 using 1:3 with lines ls 1 title "alpha=5",\
     infile2 using 1:3 with lines ls 2 title "alpha=50",\
     infile3 using 1:3 with lines ls 3 title "dense"
set size ratio 1.5
set origin 0.25,0
set yrange [0:40]
plot infile4 using 1:3 with lines ls 1 title "alpha=5",\
     infile5 using 1:3 with lines ls 2 title "alpha=50",\
     infile6 using 1:3 with lines ls 3 title "dense"
unset multiplot
