A SLBW model and doppler broadened C application.

Calculates the (n,gamma) (n,n) and total cross sections for U-238 based
on the first 3 s-wave SLBW resonances and stores the data to the "data.dat"
file. The include gnuplot script "graph.gp" can be used to generate the plot.

The code can be run using the included makefile as follows:

>$ make
>$ make run

If gnuplot is installed, a plot of the output can be generated using:

>$ make graph
