#set xrange [-0.7:20]
set xrange [-0.7:100.7]
set yrange [0:20]
set xtics 5, 10

set terminal png enhanced # size 300 300
set boxwidth 10
set samples 1001  # high quality
set border 31 linewidth .3 # thin border
set xlabel "unique transitions / %"
set ylabel "occurence / %"

set output "unique_ratios.png"
plot "unique_ratios.csv" notitle with boxes lt -1 lw 2
