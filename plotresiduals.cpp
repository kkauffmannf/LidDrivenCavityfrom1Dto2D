
/*
 * plotresiduals.cpp
 *
 * Plots the momentum residual vs the iteration number
 *
 *
 *  Created on: Jun , 2017
 *      Author: Karla Kauffmann
 *      E-mail: karla.kauffmann@cchen.cl
 *
 *
 *
 */

#include <iostream>
#include <string>
#include <vector>

//This library is used to plot with gnuplot in run time
#include "gnuplot-iostream.h"
#include "globals.h"

using namespace std;

void plotresiduals()
{

	// initializes the gnuplot figure gp
	Gnuplot gp;

	// vector that stores the values of the residual
	vector< vector<double> > values(i_iter,vector<double>(3));

	for(int counter1=0;counter1<i_iter;counter1++) {

	 	values[counter1][0]= counter1;
	 	values[counter1][1]= x_momentum_residual_sum[counter1];
	 	values[counter1][2]= pressure_residual_sum[counter1];

	}

	//set axis fonts and titles
	gp << "set xlabel 'Number of iterations' font 'Times-Roman,16'\n";
	gp << "set ylabel 'Value of the residuals' font 'Times-Roman,12'\n";
	gp << "set xtics font 'Times-Roman,14'\n";
	gp << "set ytics font 'Times-Roman,14'\n";
	gp << "set title 'Value of the residuals on each iteration' font 'Times-Roman,20'\n";
//	gp << "unset key\n";

	//plot the graphic

	/* setting the logarithmic scale in y-axis */
	gp << "set logscale y\n";

	/* This line is to ensure the correct viewing of superscripts in axis, since in
	 * some devices it fails to interpret 10^{x} as a superscript */
	gp << "set termopt enhanced\n";

	/* to better view the logarithmic scale as powers of 10 */
	gp << "set format y '10^{%L}'\n";

	/* plotting the values */
	gp << "plot " << gp.file1d(values) << "u 1:2 w lp title 'x-momentum', " << gp.file1d(values) << "u 1:3 w lp title 'pressure'\n";

    //Generates pdf figure
//	gp << "set term pdf\nset output 'momentumresidual.pdf'\nreplot\nset term x11" << endl;

}
