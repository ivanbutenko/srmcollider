set terminal png enhanced # size 300 300
set samples 1001  # high quality
set border 31 linewidth .3 # thin border
set xlabel "unique transitions"
set ylabel "peptides"
set output "unique_transitions.png"
#set logscale y

set xrange [-0.7:50.7]
set style line 1 lw 5 lt 7
plot "unique_transitions.csv" notitle lt -1 lw 2
