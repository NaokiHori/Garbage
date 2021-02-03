
set xlabel 'k'
set ylabel 'E_{ii}'

set logscale x
set logscale y

plot \
  'output/E11.dat' u 1:2 t 'E11' lc rgb '#FF0000' dt 1 w l, \
  'output/E22.dat' u 1:2 t 'E22' lc rgb '#0000FF' dt 1 w l

