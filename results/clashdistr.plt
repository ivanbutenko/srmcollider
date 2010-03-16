set terminal png enhanced # size 300 300
set boxwidth 1
set samples 1001  # high quality
set border 31 linewidth .3 # thin border
set xlabel "# clashes"
set ylabel "occurence"
set xrange [-0.6:50]
set yrange [-0.6:30000]
#set xtics 1000, 500
set output "mzclashes.png"
plot "pairlendist.csv" notitle with boxes lw 2
