The problem consists in a square (or rectangular) cavity filled with a fluid, in which
the top lid is moving at constant velocity while the other walls remain still.
The grid is uniform.
We use a backward-staggered grid with pressure nodes and velocity nodes in between.
Using the SIMPLE algorithm we solve the discretized momentum and pressure
correction equations with an underrelaxation factor for optimizing the convergence.

It uses the libraries Eigen and gnuplot-iostream.h

Compile it with:
g++ -o -std=c++11 example example.cpp -lboost_iostreams -lboost_system -lboost_filesystem
