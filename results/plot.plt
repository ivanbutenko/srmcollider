set terminal png enhanced # size 300 300
set samples 1001  # high quality
set border 31 linewidth .3 # thin border
set xlabel "m/z"
set ylabel "occurence"
set xtics 1000, 500
set output "mzdist.png"
plot [800:5000] "mzdist.csv" notitle
