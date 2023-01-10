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

########################################################################################

set xlabel "La taille du matrice"
set ylabel "Temps d'éxécution (s)"
set terminal png size 1000,600

set title "Comparaison entre les methodes iteratives"

set output "perf_iter.png"

plot "perf_iter.dat" using 2:xtic(1) title "RICHARDSON ALPHA" with lines lt rgb "blue",\
     "" using 3 title "RICHARDSON JACOBI" with lines lt rgb "green",\
     "" using 4 title "RICHARDSON GAUSS" with lines lt rgb "red"

########################################################################################

set xlabel "La taille du matrice"
set ylabel "Temps d'éxécution (s)"
set terminal png size 1000,600

set title "l’erreur par rapport à la solution analytique"

set output "perf_err.png"

plot "perf_err.dat" using 2:xtic(1) title "RICHARDSON ALPHA" with lines lt rgb "blue",\
     "" using 3 title "RICHARDSON JACOBI" with lines lt rgb "green",\
     "" using 4 title "RICHARDSON GAUSS" with lines lt rgb "red"

########################################################################################

set xlabel "La taille du matrice"
set ylabel "Temps d'éxécution (s)"
set terminal png size 1000,600

set title "Comparaison entre le nombre d'iteration"

set output "perf_nbr.png"

plot "perf_nbr.dat" using 2:xtic(1) title "RICHARDSON JACOBI" with lines lt rgb "blue",\
     "" using 3 title "RICHARDSON GAUSS" with lines lt rgb "green"




