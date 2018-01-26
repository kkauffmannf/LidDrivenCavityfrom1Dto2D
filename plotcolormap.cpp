/* plotcolormap.cpp
 *
 * Plots a bounded heatmap (color map)
 * The first two columns of the array "values"
 *
 *
 *  Created on: May 9, 2017
 *      Author: Karla Kauffmann
 *      E-mail: karla.kauffmann@cchen.cl
 */

#include <iostream>
#include <string>
#include <vector>

//This library is used to plot with gnuplot in run time
#include "gnuplot-iostream.h"
#include "globals.h"

using namespace std;

void plotcolormap()
{

	// initializes the gnuplot figure gp
	Gnuplot gp;

	double Ny = 200;
	double pi = atan(1)*4; //definition of pi
    double delta_x = L/(N-1); // delta x
    double delta_y; // delta y
    vector<double> Area_y(N); // Cross area, that decreases linearly as a function of x
    for (int i=0;i<N;i++) {
 	   Area_y[i] = Area_in - ((Area_in-Area_out)/(N-1)*i);
    }

    // vector that stores the values of x y and z in 3 columns
    int index;
    vector< vector<double> > values((N*Ny), vector<double>(3));

    // intializing with zeroes
    for(int counter1=0;counter1<(N*Ny);counter1++) {
    		  values[counter1][0]=0.0;
    		  values[counter1][1]=0.0;
    		  values[counter1][2]=0.0;
       	 }

   // for each value of x, y loops from -r to r where r is the radius of the area
   // corresponding to that position x.
   for(int counter1=0;counter1<(N-1);counter1++) {
	   for(int counter2=0;counter2<Ny;counter2++){
		  index = Ny*counter1+counter2;
		  delta_y = 2*sqrt(Area_y[counter1]/pi)/(Ny-1);
		  values[index][0]=counter1*delta_x;
	      values[index][1]=delta_y*(-(Ny-1)/2.0+counter2);
	      values[index][2]= velocity[counter1];
	      if(counter1>=(N-2)){
	      values[index][2]= velocity[N-2];
	      }
       }
   	 }

	gp << "set view map\n";
	gp << "set dgrid3d 200,200,2\n";

	//quality of the grid. Ten times interpolated in x and y.
	//Default is too "pixelated"
//	gp << "set pm3d interpolate 10,10\n";

	//set color palette to the matlab default
	gp << "set palette defined (0 0 0 0.5, 1 0 0 1, 2 0 0.5 1, 3 0 1 1, 4 0.5 1 0.5, 5 1 1 0, 6 1 0.5 0, 7 1 0 0, 8 0.5 0 0)\n";

	//set color palette to greyscale
    //gp << "set palette grey\n";

	//set axis fonts and titles
	gp << "set xlabel 'position in x [m]' font 'Times-Roman,16'\n";
	gp << "set ylabel 'position in y [m]' font 'Times-Roman,16'\n";
	gp << "set xtics font 'Times-Roman,14'\n";
	gp << "set ytics font 'Times-Roman,14'\n";
	gp << "set title 'Velocity of the fluid in a nozzle' font 'Times-Roman,20'\n";
	gp << "unset key\n";

	//plot with contours
//	gp << "set parametric\n";
//	gp << "unset surface\n";
//	gp << "set contour base\n";
//	gp << "set cntrparam levels auto 40\n";
//	gp << "unset key\n";

    // store color plot of values[][] in a temporary file tmp.txt
	gp << "set table 'tmp.txt'\n";
	gp << "splot " << gp.file1d(values) << " using 1:2:3\n";
	gp << "unset table\n";

	// set left, right, upper and lower bounds of plot
	// by setting z to NaN if out of bounds.
	gp << "Areay(x) = " << Area_in << " - (" << Area_in << " - " << Area_out << ")/(" << L << " - 0.0)*x\n";
	gp << "bounds(x,y,z) = x > " << L << " || x < 0.0 || y > sqrt(Areay(x)/" << pi << ") || y < -sqrt(Areay(x)/" << pi << ") ? NaN : z\n";
	gp << "set multiplot\n";

	//plot the graphic
	gp << "plot 'tmp.txt' u 1:2:(bounds($1,$2,$3))  w image\n";
	gp << "unset multiplot\n";

    //Generates pdf figure
	gp << "set term pdf\nset output 'velocity_nozzle.pdf'\nreplot\nset term x11" << endl;

	// delete the temporary file
	gp << "!rm 'tmp.txt'\n";

}
