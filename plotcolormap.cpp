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

//	double pi = atan(1)*4; //definition of pi

    // vector that stores the values of x y and velocity in 3 columns
    vector< vector<double> > velocities((Nx*Ny), vector<double>(3));

    // vector that stores the values of x y and pressure in 3 columns
    vector< vector<double> > pressures((Nx*Ny), vector<double>(3));

    // intializing with zeroes
    for(int counter1=0;counter1<(Nx*Ny);counter1++) {
    	velocities[counter1][0]=0.0;
    	velocities[counter1][1]=0.0;
    	velocities[counter1][2]=0.0;
    }

   // for each value of x, y loops from -r to r where r is the radius of the area
   // corresponding to that position x.
   int index=0;
   for(int counter1=0;counter1<Ny;counter1++) {
	   for(int counter2=0;counter2<Nx;counter2++){
		   velocities[index][0]=position_v_velocity_node_y[counter1];
		   velocities[index][1]=position_u_velocity_node_x[counter2];
		   velocities[index][2]=sqrt(u_velocity[counter2][counter1]*u_velocity[counter2][counter1] + v_velocity[counter2][counter1]*v_velocity[counter2][counter1]);
//		   velocities[index][2]=u_velocity[counter2][counter1];
		   index++;
       }
   }

	gp << "set view map\n";
	gp << "set dgrid3d 200,200\n";

	//quality of the grid. Ten times interpolated in x and y.
	//Default is too "pixelated"
//	gp << "set pm3d interpolate 10,10\n";
	gp << "set pm3d at b map\n";

	//set color palette to the matlab default
	gp << "set palette defined (0 0 0 0.5, 1 0 0 1, 2 0 0.5 1, 3 0 1 1, 4 0.5 1 0.5, 5 1 1 0, 6 1 0.5 0, 7 1 0 0, 8 0.5 0 0)\n";

	//set color palette to greyscale
    //gp << "set palette grey\n";

	//set axis fonts and titles
	gp << "set xlabel 'position in x [m]' font 'Times-Roman,16'\n";
	gp << "set ylabel 'position in y [m]' font 'Times-Roman,16'\n";
	gp << "set xtics font 'Times-Roman,14'\n";
	gp << "set ytics font 'Times-Roman,14'\n";
	gp << "set title 'Velocity of the fluid in the cavity' font 'Times-Roman,20'\n";
	gp << "unset key\n";

	//set ranges
	gp << "set xrange[0.0:" << Lx << "]\n";
	gp << "set yrange[0.0:" << Ly << "]\n";

//	//plot with contours
////	gp << "set parametric\n";
////	gp << "unset surface\n";
////	gp << "set contour base\n";
////	gp << "set cntrparam levels auto 40\n";
////	gp << "unset key\n";
//

//	gp << "splot " << gp.file1d(velocities) << " using 2:1:3\n";
	gp << "plot " << gp.file1d(velocities) << " using 2:1:3 with image pixels\n" << endl;

    //Generates pdf figure
//	gp << "set term pdf\nset output 'velocity_cavity.pdf'\nreplot\nset term x11" << endl;

}
