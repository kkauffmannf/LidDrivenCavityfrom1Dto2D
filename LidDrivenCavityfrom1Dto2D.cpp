/*
 * LidDrivenCavityfrom1Dto2D.cpp
 *
 *  The problem consists in a square (or rectangular) cavity filled with a fluid, in which
 * the top lid is moving at constant velocity while the other walls remain still.
 * The grid is uniform.
 * We use a backward-staggered grid with pressure nodes and velocity nodes in between.
 * Using the SIMPLE algorithm we solve the discretized momentum and pressure
 * correction equations with an underrelaxation factor for optimizing the convergence.
 * *
 * It uses the libraries Eigen and gnuplot-iostream.h
 *
 *
 * Compile it with:
 * g++ -o -std=c++11 example example.cpp -lboost_iostreams -lboost_system -lboost_filesystem
 *
 *  Created on: May 22, 2017
 *      Author: Karla Kauffmann
 *      E-mail: karla.kauffmann@cchen.cl
 */

#include <iostream>
#include <vector>

// These libraries are used to measure the execution time of the program
#include <chrono>
#include <ctime>

// This library is used to maximize the output precision
#include <limits>

//This library is used for linear algebra calculations (matrix solver)
#include <Eigen/Dense>

#include "globals.h" /* global variables are declared here to be used across the different source files */

using namespace std;
using namespace Eigen;
using namespace std::chrono;

/* These are the global variables declared in global.h */

/* Variables in input.txt to be defined at the start by the user. */
int N = 5; /* number of pressure nodes for the CV computation */
double L = 2.0; /* length of the rod */
double p_in = 10.0; /* stagnation pressure at inlet of the nozzle */
double p_out = 0.0; /* static pressure at the exit of the nozzle */
double density = 1.0; /* density of the fluid */
double Area_in = 0.5; /* cross-sectional area at the inlet */
double Area_out = 0.1; /* cross-sectional area at the exit */
double m_dot = 1.0; /* guessed mass flow rate to generate the an initial velocity field */
double urfu = 0.8; /* under-relaxation factor for velocity u. Cannot be 0.0 since it goes in the denominator
                      of a_p in the momentum equations */
double urfp = 0.8; /* under-relaxation factor for pressure */
double residual_threshold = 1.0E-5; /* the value of the residual to end the iterations */

/* Physical quantities */
vector<double> position_pressure_node; /* position of the pressure nodes */
vector<double> position_velocity_node; /* position of the velocity nodes*/
vector<double> Area_pressure_node;  /* cross section at pressure nodes */
vector<double> Area_velocity_node; /* cross section at velocity nodes */
vector<double> velocity;
vector<double> velocity_old; /* the value of the velocity in the previous iteration */
vector<double> pressure;
vector<double> pressure_old; /* the value of the pressure in the previous iteration */
vector<double> dx; /* parameter d for the pressure correction equation */

/* Variables used for the iterations */
int i_iter = 0; /* number of iterations */
int MAX_ITER = 1000000; /* set the maximum number of iterations to store in the residual vector */
vector<double> x_momentum_residual_sum; /* sum of the residuals of the x-momentum equation per iteration*/
vector<double> pressure_residual_sum; /* sum of the residuals of the pressure equation per iteration*/


/* Main execution of the program */
int main()
{

   /* Display numbers with maximum precision possible */
	cout.precision(numeric_limits<double>::digits10 + 1);


   /* Variable for tracking execution time */
	high_resolution_clock::time_point t1 = high_resolution_clock::now();

   /* reads the file input.txt and stores the values */
   readinput();


   /* Physical quantities vector size definition. They are resized here because they depend on the input value of N.
    * Pressure and velocity nodes vectors have different lengths because a staggered grid is used  */
   position_pressure_node.resize(N);
   position_velocity_node.resize(N-1);
   Area_pressure_node.resize(N);
   Area_velocity_node.resize(N-1);
   velocity.resize(N-1);
   velocity_old.resize(N-1);
   pressure.resize(N);
   pressure_old.resize(N);
   dx.resize(N-1);

   /* The residuals are stored for each iteration with a maximum iteration number of MAX_ITER */
   x_momentum_residual_sum.resize(MAX_ITER);
   pressure_residual_sum.resize(MAX_ITER);


   /* Initializes the variables */
   initialization();


   /* We obtain the guessed velocities (u_star*) by solving the system of
	* momentum equations. */
   VectorXd u_star;
   momentum_equation_solve(u_star,i_iter);


   /* Now that we solved the momentum equations we have to solve the
    * pressure correction equations
    */


   /* for the zeroth iteration we don't have pressure correction yet, so we use the
    * initially guessed pressure */
   VectorXd pressure_prime = VectorXd::Map(pressure.data(), pressure.size());
   pressure_correction_equation_solve(u_star,pressure_prime,i_iter);

   /* We have the velocities from the momentum equation and the pressures corrections,
    * so we now proceed to correct the pressure and velocities  */
   correct_pressure_and_velocities(u_star,pressure_prime);

   /* We apply the underrelaxation factors for velocity and pressure */
   underrelaxation();





      /* Start next iteration until the residuals are lower than the threshold */
      while( (x_momentum_residual_sum[i_iter] > residual_threshold) || (pressure_residual_sum[i_iter] > residual_threshold) ) {


    	      /* the new velocities and pressures are the old ones in the next iteration */
    	      for(int i=0;i<(N-1);i++){
    	         velocity_old[i] = velocity[i];
    	      }

    	      for(int i=0;i<N;i++){
    	         pressure_old[i] = pressure[i];
    	      }


    	      /* Solving the momentum equation */
    	      VectorXd u_star;
              momentum_equation_solve(u_star,(i_iter + 1));


              /* Solving the pressure equation */
              pressure_correction_equation_solve(u_star,pressure_prime,(i_iter + 1));


              /* Correcting the pressure and velocity */
              correct_pressure_and_velocities(u_star,pressure_prime);


              /* Applying the underrelaxation factor */
              underrelaxation();


         /* Prints out the residuals */
	     cout << i_iter << "\t" << x_momentum_residual_sum[i_iter] << "\t" << pressure_residual_sum[i_iter] << endl;

	     /* Advancing the iteration number */
	     i_iter++;

	 /*************************** */


      }

      /* print out the position, the velocity and pressure */
      cout << endl;
      cout << "# position" << "\t" << "velocity" << "\t" << "pressure" << endl;
      for(int i=0;i<(N-1);i++){
    	  cout << position_velocity_node[i] << "\t" << velocity[i] << "\t" << pressure[i] << endl;
   	      }

     /* color map of the velocities in 1 D */
     plotcolormap();


     /* Residuals vs iterations */
     plotresiduals();

     /* Ending the execution and printing the execution time */
     high_resolution_clock::time_point t2 = high_resolution_clock::now();

     duration<double> time_span = duration_cast<duration<double>>(t2 - t1);

     cout << endl;
     cout << "It took me " << time_span.count() << " seconds.";
     cout << endl;


}
