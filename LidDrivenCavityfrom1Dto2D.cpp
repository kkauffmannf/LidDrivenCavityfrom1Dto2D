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
#include <fstream>
#include <vector>

// These libraries are used to measure the execution time of the program
#include <chrono>
#include <ctime>
#include <cmath>

// This library is used to maximize the output precision
#include <limits>

////This library is used for linear algebra calculations (matrix solver)
//#include <Eigen/Dense>

#include "globals.h" /* global variables are declared here to be used across the different source files */

using namespace std;
//using namespace Eigen;
using namespace std::chrono;

/* These are the global variables declared in global.h */

/* Variables in input.txt to be defined at the start by the user. */
int Nodesx = 130; /* number of pressure nodes for the CV computation in the x direction */
int Nodesy = 130; /* number of pressure nodes for the CV computation in the y direction */
int ngc = 1; /* number of guard cells */
int Npx = Nodesx + 2*ngc; /* number of pressure nodes plus guard cells in the x direction */
int Npy = Nodesy + 2*ngc; /* number of pressure nodes plus guard cells in the y direction */
int Nux = Nodesx + 2*ngc + 1; /* number of u-velocity nodes plus guard cells in the x direction */
int Nuy = Nodesy + 2*ngc; /* number of u-velocity plus guard cells in the y direction */
int Nvx = Nodesx + 2*ngc; /* number of v-velocity nodes plus guard cells in the x direction */
int Nvy = Nodesy + 2*ngc + 1; /* number of v-velocity nodes plus guard cells in the y direction */
double Lx = 1.0; /* length of the cavity in the x direction */
double Ly = 1.0; /* length of the cavity in the y direction */
double p_init = 0.0;
double density = 1.0; /* density of the fluid */
double lid_velocity = 1.0; /* velocity of the lid */
double Reynolds_num = 100.0; /* Reynolds number */
double urfu = 0.2; /* under-relaxation factor for velocity u. Cannot be 0.0 since it goes in the denominator
                      of a_p in the momentum equations */
double urfv = 0.2; /* under-relaxation factor for velocity v. Cannot be 0.0 since it goes in the denominator
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
//int MAX_ITER = 4; /* set the maximum number of iterations to store in the residual vector */
vector<double> x_momentum_residual_sum; /* sum of the residuals of the x-momentum equation per iteration*/
vector<double> y_momentum_residual_sum; /* sum of the residuals of the y-momentum equation per iteration*/
vector<double> pressure_residual_sum; /* sum of the residuals of the pressure equation per iteration*/
double x_momentum_residual_sum_norm; /* norm sum of the residuals of the x-momentum equation (1st iteration)*/
double y_momentum_residual_sum_norm; /* norm of sum of the residuals of the y-momentum equation (1st iteration)*/
double pressure_residual_sum_norm; /* norm of the sum of the residuals of the pressure equation (1st iteration)*/

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
   Npx = Nodesx + 2*ngc;
   Npy = Nodesy + 2*ngc;
   Nux = Nodesx + 1 + 2*ngc;
   Nuy = Nodesy + 2*ngc;
   Nvx = Nodesx + 2*ngc;
   Nvy = Nodesy + 1 + 2*ngc;
   position_pressure_node_x.resize(Npx);
   position_pressure_node_y.resize(Npy);
   position_u_velocity_node_x.resize(Nux);
   position_u_velocity_node_y.resize(Nuy);
   position_v_velocity_node_x.resize(Nvx);
   position_v_velocity_node_y.resize(Nvy);
   Area_velocity_node_u.resize(Nux,vector<double>(Nuy));
   Area_velocity_node_v.resize(Nvx,vector<double>(Nvy));
   u_velocity.resize(Nux,vector<double>(Nuy));
   u_velocity_old.resize(Nux,vector<double>(Nuy));
   v_velocity.resize(Nvx,vector<double>(Nvy));
   v_velocity_old.resize(Nvx,vector<double>(Nvy));
   pressure.resize(Npx,vector<double>(Npy));
   pressure_old.resize(Npx,vector<double>(Npy));
   d_u.resize(Nux,vector<double>(Nuy));
   d_v.resize(Nvx,vector<double>(Nvy));

   /* The residuals are stored for each iteration with a maximum iteration number of MAX_ITER */
   x_momentum_residual_sum.resize(MAX_ITER);
   y_momentum_residual_sum.resize(MAX_ITER);
   pressure_residual_sum.resize(MAX_ITER);

   if(ngc > 1){
	   cout << "ngc > 1 has not been implemented yet. Please set ngc = 1" << endl;
	   return 0;
   }


   /* Initializes the variables */
   initialization();

   /* Imposes the boundary conditions on u_velocity, v_velocity and pressure guard cells and borders */
   boundary_conditions();

//   /* We obtain the guessed velocities (u_star*) by solving the system of
//	* momentum equations. */
   vector<vector<double>> u_star(Nux, vector<double>(Nuy));
   vector<vector<double>> v_star(Nux, vector<double>(Nuy));
//   MatrixXd u_star;
//   MatrixXd v_star;
//   u_star.resize(Nx,Ny);
//   u_star=MatrixXd::Zero(Nx,Ny);
//   v_star.resize(Nx,Ny);
//   v_star=MatrixXd::Zero(Nx,Ny);

   momentum_equation_solve(u_star, v_star, i_iter);








   //////////////////////////////////////////////////////////////////////////////


   /* Now that we solved the momentum equations we have to solve the
    * pressure correction equations
    */

   vector<vector<double>> pressure_prime(Npx, vector<double>(Npy));

   pressure_correction_equation_solve(u_star,v_star, pressure_prime,i_iter);

   ///////////////////////////////////////////////////////////////////////////////////////

   /* We have the velocities from the momentum equation and the pressures corrections,
    * so we now proceed to correct the pressure and velocities  */
   correct_pressure_and_velocities(u_star,v_star,pressure_prime);

      /* Start next iteration until the residuals are lower than the threshold */
      while((i_iter < MAX_ITER) && ((x_momentum_residual_sum[i_iter] > residual_threshold) || (y_momentum_residual_sum[i_iter] > residual_threshold) || (pressure_residual_sum[i_iter] > residual_threshold)) ) {

			  /* Imposes the boundary conditions on u_velocity, v_velocity and pressure guard cells and borders */
			  boundary_conditions();

    	      /* the new velocities and pressures are the old ones in the next iteration */
    	      for(int i=0;i<Npx;i++){
    	    	  for(int j=0;j<Npy;j++){

    	    		  pressure_old[i][j] = pressure[i][j];
    	    	  }
    	      }

    	      for(int i=0;i<Nux;i++){
    	    	  for(int j=0;j<Nuy;j++){
    	    		  u_velocity_old[i][j] = u_velocity[i][j];
    	    		  u_star[i][j]=u_velocity_old[i][j];
    	    	  }
    	      }

    	      for(int i=0;i<Nvx;i++){
    	    	  for(int j=0;j<Nvy;j++){
    	    		  v_velocity_old[i][j] = v_velocity[i][j];
    	    		  v_star[i][j]=v_velocity_old[i][j];
    	    	  }
    	      }

    	      /* Solving the momentum equation */
              momentum_equation_solve(u_star,v_star,(i_iter + 1));

              vector<vector<double>> pressure_prime(Npx, vector<double>(Npy));

              /* Solving the pressure equation */
              pressure_correction_equation_solve(u_star,v_star,pressure_prime,(i_iter + 1));

              /* Correcting the pressure and velocity */
              correct_pressure_and_velocities(u_star,v_star,pressure_prime);


         /* Prints out the residuals */
	     cout << i_iter << "\t" << x_momentum_residual_sum[i_iter] <<  "\t" << y_momentum_residual_sum[i_iter] << "\t" << pressure_residual_sum[i_iter] << endl;

	     /* Advancing the iteration number */
	     i_iter++;
	     if (i_iter==MAX_ITER){
	    	 cout << "The program has reached the maximum of iterations allowed" << endl;
	     }

	 /*************************** */
      }

  	  /* Imposes the boundary conditions on u_velocity, v_velocity and pressure guard cells and borders */
  	  boundary_conditions();
//
//////      /* print out the position, the velocity and pressure */
//////      cout << endl;
//////      cout << "# position" << "\t" << "velocity" << "\t" << "pressure" << endl;
//////      for(int i=0;i<(Nx-1);i++){
//////    	  cout << position_velocity_node[i] << "\t" << velocity[i] << "\t" << pressure[i] << endl;
//////   	      }
////

//  	   for(int i=0;i<Nux;i++){
//  		   for(int j=0;j<Nuy;j++){
//  			   cout << u_velocity[i][j] <<" ";
//  		   }
//  		   cout << endl;

//  	   }
     /* color map of the velocities in 2 D */
     plotcolormap();


     /* Residuals vs iterations */
     plotresiduals();

     /* Ending the execution and printing the execution time */
     high_resolution_clock::time_point t2 = high_resolution_clock::now();

     duration<double> time_span = duration_cast<duration<double>>(t2 - t1);

     cout << endl;
     cout << "It took me " << time_span.count() << " seconds.";
     cout << endl;

     /* write to a file the u_velocity along vertical line through center
      * write to a file the v_velocity along horizontal line through center*/
     ofstream myfile[2];

     myfile[0].open("u_velocity_center.txt");
     myfile[1].open("v_velocity_center.txt");
     if (myfile[0].is_open() && myfile[1].is_open())
     {
         for(int j=ngc; j<(Nuy - ngc); j++) {
      	     myfile[0] << (j - ngc) * Ly/(Nodesy-0.5) << "\t" << u_velocity[floor(Nux/2)][j] << "\n";
         }
         for(int i=ngc; i<(Nvx - ngc); i++) {
        	 myfile[1] << (i - ngc) * Lx/(Nodesx-0.5) << "\t" << v_velocity[i][floor(Nvy/2)] << "\n";
         }

         myfile[0].close();
         myfile[1].close();
     }
     else cout << "Unable to open file";
     return 0;

}
