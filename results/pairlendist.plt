#reset
#blue = "#babaff"
#set isosample 2, 100
#set table 'pairlendist.csv'
#splot 1-exp(-y/2.0)
#unset table
#
#unset colorbox

set xrange [-0.7:20]
set yrange [0:5]

#A = 4
#B = 88
#C = 0
#D = 5
#E = -0.6
#
#eps=0.05*E
#eps2=0.005*(D-C)
#
#f(x) = ( x < 5 ? x : x - B )

set terminal png enhanced # size 300 300
set boxwidth 1
set samples 1001  # high quality
set border 31 linewidth .3 # thin border
set xlabel "# clashes"
set ylabel "occurence / percent"

#set ytics 0, 1, 4
#set ytics add (gprintf("%.0f", 6+B) 6)
#set arrow 1 from A-eps2,     -E     to A+eps2, -E           nohead lc rgb "#ffffff" front
#set arrow 2 from A-eps2,     E      to A+eps2, E            nohead lc rgb "#ffffff" front
#set arrow 3 from A-eps-eps2, -E-eps to A+eps-eps2, -E+eps   nohead front
#set arrow 4 from A-eps+eps2, -E-eps to A+eps+eps2, -E+eps   nohead front
#set arrow 5 from A-eps-eps2, E-eps  to A+eps-eps2, E+eps    nohead front
#set arrow 6 from A-eps+eps2, E-eps  to A+eps+eps2, E+eps    nohead front
#set xtics 1000, 500
set output "pairlen.png"
#plot [3:50] "pairlendist.csv" notitle with boxes lw 2, [0:3] "pairlendist.csv" notitle with boxes lw 2
plot "pairlendist.csv" notitle with boxes lw 2
#plot "pairlendist.csv" u (f($2))



