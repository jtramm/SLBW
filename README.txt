A SLBW model and doppler broadened C application.

Calculates the (n,gamma) (n,n) and total cross sections for U-238 based
on the first 3 s-wave SLBW resonances and stores the data to the "data.dat"
file. The include gnuplot script "graph.gp" can be used to generate the plot.

I'm using my own custom Faddeeva function, as the MIT implementation is way
too slow and is too heavy weight for this application, as the imaginary
component of Z is always positive and |z| is usually fairly large. 
See more details at https://github.com/jtramm/FNF

The code can be run using the included makefile as follows:

>$ make
>$ make run

If gnuplot is installed, a plot of the output can be generated using:

>$ make graph
