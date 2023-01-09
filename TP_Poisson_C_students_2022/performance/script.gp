set xlabel "La taille du matrice"
set ylabel "Temps d'éxécution (s)"
set terminal png size 1000,600

set title "Comparaison entre les methodes LAPACK"

set output "perf_Lapack.png"

plot "perf.dat" using 2:xtic(1) title "DGBTRF + DGBTRS" with lines lt rgb "blue",\
     "" using 3 title "DGBSV" with lines lt rgb "green"

########################################################################################

set xlabel "La taille du matrice"
set ylabel "Temps d'éxécution (s)"
set terminal png size 1000,600

set title "Comparaison entre la factorisation LU avec et sans DGBTRF"

set output "perf_LU.png"

plot "perf_LU.dat" using 2:xtic(1) title "Sans DGBTRF" with lines lt rgb "blue",\
     "" using 3 title "Avec DGBTRF" with lines lt rgb "green"



