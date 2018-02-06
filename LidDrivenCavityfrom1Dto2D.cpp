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
int Nx = 130; /* number of pressure nodes for the CV computation in the y direction */
int Ny = 130; /* number of pressure nodes for the CV computation in the x direction */
double Lx = 1.0; /* length of the cavity in the x direction */
double Ly = 1.0; /* length of the cavity in the y direction */
double p_init = 0.0;
double density = 1.0; /* density of the fluid */
double lid_velocity = 1.0; /* velocity of the lid */
double Reynolds_num = 100.0; /* Reynolds number */
double urfu = 0.8; /* under-relaxation factor for velocity u. Cannot be 0.0 since it goes in the denominator
                      of a_p in the momentum equations */
double urfv = 0.8; /* under-relaxation factor for velocity v. Cannot be 0.0 since it goes in the denominator
                      of a_p in the momentum equations */
double urfp = 0.8; /* under-relaxation factor for pressure */

double residual_threshold = 1.0E-5; /* the value of the residuals to end the iterations */

/* Physical quantities */
vector<double> position_pressure_node_x; /* position of the pressure nodes in x */
vector<double> position_pressure_node_y; /* position of the pressure nodes in y */
vector<double> position_u_velocity_node_x; /* position of the u velocity nodes in x */
vector<double> position_u_velocity_node_y; /* position of the u velocity nodes in y */
vector<double> position_v_velocity_node_x; /* position of the v velocity nodes in x */
vector<double> position_v_velocity_node_y; /* position of the v velocity nodes in y */
vector<vector<double>> Area_pressure_node; /* Cross-section for the pressure nodes */
vector<vector<double>> Area_velocity_node_u; /* Cross-section for the velocity nodes for u velocity with (i,J) indexes */
vector<vector<double>> Area_velocity_node_v; /* Cross-section for the velocity nodes for v velocity with (I,j) indexes */
vector<vector<double>> u_velocity; /* velocity in the x direction */
vector<vector<double>> u_velocity_old; /* the value of the u_velocity in the previous iteration */
vector<vector<double>> v_velocity; /* velocity in the y direction */
vector<vector<double>> v_velocity_old; /* the value of the v_velocity in the previous iteration */
vector<vector<double>> pressure;
vector<vector<double>> pressure_old; /* the value of the pressure in the previous iteration */
vector<vector<double>> d_u; /* parameter d for the pressure correction equation for u velocity with (i,J) indexes */
vector<vector<double>> d_v; /* parameter d for the pressure correction equation for v velocity with (I,j) indexes */

/* Variables used for the iterations */
int i_iter = 0; /* number of iterations */
//int MAX_ITER = 1000000; /* set the maximum number of iterations to store in the residual vector */
int MAX_ITER = 1000; /* set the maximum number of iterations to store in the residual vector */
vector<double> x_momentum_residual_sum; /* sum of the residuals of the x-momentum equation per iteration*/
vector<double> y_momentum_residual_sum; /* sum of the residuals of the y-momentum equation per iteration*/
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
   position_pressure_node_x.resize(Nx);
   position_pressure_node_y.resize(Ny);
   position_u_velocity_node_x.resize(Nx);
   position_u_velocity_node_y.resize(Ny);
   position_v_velocity_node_x.resize(Nx);
   position_v_velocity_node_y.resize(Ny);
   Area_pressure_node.resize(Nx,vector<double>(Ny));
   Area_velocity_node_u.resize(Nx,vector<double>(Ny));
   Area_velocity_node_v.resize(Nx,vector<double>(Ny));
   u_velocity.resize(Nx,vector<double>(Ny));
   u_velocity_old.resize(Nx,vector<double>(Ny));
   v_velocity.resize(Nx,vector<double>(Ny));
   v_velocity_old.resize(Nx,vector<double>(Ny));
   pressure.resize(Nx,vector<double>(Ny));
   pressure_old.resize(Nx,vector<double>(Ny));
   d_u.resize(Nx,vector<double>(Ny));
   d_v.resize(Nx,vector<double>(Ny));

   /* The residuals are stored for each iteration with a maximum iteration number of MAX_ITER */
   x_momentum_residual_sum.resize(MAX_ITER);
   y_momentum_residual_sum.resize(MAX_ITER);
   pressure_residual_sum.resize(MAX_ITER);


   /* Initializes the variables */
   initialization();


   /* We obtain the guessed velocities (u_star*) by solving the system of
	* momentum equations. */
   MatrixXd u_star;
   MatrixXd v_star;
   u_star.resize(Nx,Ny);
   u_star=MatrixXd::Zero(Nx,Ny);
   v_star.resize(Nx,Ny);
   v_star=MatrixXd::Zero(Nx,Ny);
   momentum_equation_solve(u_star, v_star, i_iter);

//   cout << u_star << endl;
//   cout << v_star << endl;
//
//   for(int i=0;i<Nx;i++){
//	   for(int j=0;j<Ny;j++){
//		   u_velocity[i][j]=u_star(i,j);
//		   v_velocity[i][j]=v_star(i,j);
//	   }
//   }


   //////////////////////////////////////////////////////////////////////////////


   /* Now that we solved the momentum equations we have to solve the
    * pressure correction equations
    */


   /* for the zeroth iteration we don't have pressure correction yet, so we use the
    * initially guessed pressure. We map the array of the vector type to a MatrixXd.
    * MatrixXD pressure_prime = vector<vector<double>> pressure */
   MatrixXd pressure_prime (pressure.size(), pressure[0].size());

   /* cast unsigned int pressure.size to int to be able to compare */
   int signedIntsize = (int) pressure.size();
   for (int i = 0; i < signedIntsize; ++i){
       pressure_prime.row(i) = VectorXd::Map(&pressure[i][0], pressure[0].size());
   }
   pressure_correction_equation_solve(u_star,v_star, pressure_prime,i_iter);

   ///////////////////////////////////////////////////////////////////////////////////////


   /* We have the velocities from the momentum equation and the pressures corrections,
    * so we now proceed to correct the pressure and velocities  */
   correct_pressure_and_velocities(u_star,v_star,pressure_prime);

   /* We apply the underrelaxation factors for velocity and pressure */
   underrelaxation(pressure_prime);





      /* Start next iteration until the residuals are lower than the threshold */
      while((i_iter < MAX_ITER) && ((x_momentum_residual_sum[i_iter] > residual_threshold) || (y_momentum_residual_sum[i_iter] > residual_threshold) || (pressure_residual_sum[i_iter] > residual_threshold)) ) {

    	      /* the new velocities and pressures are the old ones in the next iteration */
    	      for(int i=0;i<Nx;i++){
    	    	  for(int j=0;j<Ny;j++){
    	    		  u_velocity_old[i][j] = u_velocity[i][j];
    	    		  v_velocity_old[i][j] = v_velocity[i][j];
    	    		  pressure_old[i][j] = pressure[i][j];
    	    	  }
    	      }

    	      /* Solving the momentum equation */
    	      u_star=MatrixXd::Zero(Nx,Ny);
    	      v_star=MatrixXd::Zero(Nx,Ny);
              momentum_equation_solve(u_star,v_star,(i_iter + 1));


              /* Solving the pressure equation */
              pressure_correction_equation_solve(u_star,v_star,pressure_prime,(i_iter + 1));


              /* Correcting the pressure and velocity */
              correct_pressure_and_velocities(u_star,v_star,pressure_prime);


              /* Applying the underrelaxation factor */
              underrelaxation(pressure_prime);


         /* Prints out the residuals */
	     cout << i_iter << "\t" << x_momentum_residual_sum[i_iter] <<  "\t" << y_momentum_residual_sum[i_iter] << "\t" << pressure_residual_sum[i_iter] << endl;

	     /* Advancing the iteration number */
	     i_iter++;
	     if (i_iter==MAX_ITER){
	    	 cout << "The program has reached the maximum of iterations allowed" << endl;
	     }

	 /*************************** */


      }

//      /* print out the position, the velocity and pressure */
//      cout << endl;
//      cout << "# position" << "\t" << "velocity" << "\t" << "pressure" << endl;
//      for(int i=0;i<(Nx-1);i++){
//    	  cout << position_velocity_node[i] << "\t" << velocity[i] << "\t" << pressure[i] << endl;
//   	      }

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
