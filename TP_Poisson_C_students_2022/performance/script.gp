set xlabel "La taille du matrice"
set ylabel "Temps d'éxécution (s)"
set style data histograms
set style fill solid
set boxwidth 1
set terminal png size 1000,600

set title "Comparaison entre les methodes LAPACK"

set output "perf.png"

plot "perf.dat" using 2:xtic(1) title "DGBTRF" lt rgb "blue",\
     "" using 3 title "DGBTRS" lt rgb "green",\
     "" using 4 title "DGBSV" lt rgb "red"


