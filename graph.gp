set style data lines
set logscale xy
set key left bottom Left title 'Legend' box
set title "U-238 SLBW Data From First 3 Resonances @ 0K"
set xlabel "Energy (eV)"
set ylabel "Cross Section (b)"
plot 'data.dat' using 1:2 title 'sigma_f', 'data.dat' using 1:3 title 'sigma_n', 'data.dat' using 1:4 title 'sigma_t'
