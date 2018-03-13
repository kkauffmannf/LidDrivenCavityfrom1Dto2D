/*
 * pressure_equation_solve.cpp
 *
 * This module solves the pressure correction equation and sets the pressure matrix values
 *
 *  Created on: Jul 17, 2017
 *      Author: Karla Kauffmann
 *      E-mail: karla.kauffmann@cchen.cl
 */

#include <iostream>
#include <string>
#include <vector>
//#include <Eigen/SVD>
//#include <Eigen/Dense>
#include "globals.h"

using namespace std;

void pressure_correction_equation_solve(vector<vector<double>> u_star, vector<vector<double>> v_star, vector<vector<double>> &pressure_prime, int i_iter)
	{

		   /* Solve the discretised pressure correction equation using the guessed values of
		    * velocity and the values of pressure. The equation has the form:
		    *
		    * a_I,J * p_I,J^(') = a_I-1,J * p_I-1,J^(') + a_I+1,J * p_I+1,J^(') +
		    *                     + a_I,J-1 * p_I,J-1^(') + a_I,J+1 * p_I,J+1^(') + Source_p
		    * where p' is the pressure corrections p'=(p-p*) in each node.
		    *
		    * */

		   vector<vector<double>> a_pressure;
		   vector<vector<double>> Source_pressure;

		   /* initialize values of the pressure residual */
		   pressure_residual_sum[i_iter] = 0.0;

	 	   /* this is the matrix A that stores the a_I,J values and the vector
	  	    * b that stores the source values to solve the system of equations for
	   	    * the pressure corrections. The system of equations is given by the pressure
	   	    * corrections equations. We solve the system A*p'=b and the solutions of p'
	   	    */
		   vector<vector<double>> A_pressure(Npy, vector<double>(Npy));
		   vector<double> b_pressure(Npy, 0.0);

		   /* For pressure nodes we assume density and Gamma=1/Reynolds_num constant */

		   double a_w;
		   double a_e;
		   double a_s;
		   double a_n;
		   double S_pr;
		   double a_p;

		   for(int i=ngc; i<(Npx - ngc); i++){

			   /* For each new i, we set the A matrix and b vector to zero */
			   for (int loopp1=0;loopp1 < Npy; loopp1++){
				   for (int loopp2=0; loopp2 < Npy; loopp2++){
					   A_pressure[loopp1][loopp2] = 0.0;
				   }
				   b_pressure[loopp1] = 0.0;
			   }

			   for(int j=ngc;j<(Npy - ngc);j++){
				   a_w = density * d_u[i][j] * Area_velocity_node_u[i][j];;
				   a_e = density * d_u[i+1][j] * Area_velocity_node_u[i+1][j];
				   a_s = density * d_v[i][j] * Area_velocity_node_v[i][j];
				   a_n = density * d_v[i][j+1] * Area_velocity_node_v[i][j+1];
				   a_p = a_w + a_e + a_s + a_n;
				   S_pr = (density * u_star[i][j] * Area_velocity_node_u[i][j]) - (density * u_star[(i+1)][j] * Area_velocity_node_u[i+1][j]) + (density * v_star[i][j] * Area_velocity_node_v[i][j]) - (density * v_star[i][j+1] * Area_velocity_node_v[i][j+1]);
				   A_pressure[j][j] = a_p;
				   A_pressure[j][j+1] = - a_n;
				   A_pressure[j][j-1] = - a_s;
				   b_pressure[j] = a_e * pressure_prime[(i+1)][j] + a_w * pressure_prime[(i-1)][j] + S_pr;
//				   pressure_residual_sum[i_iter] = pressure_residual_sum[i_iter] + abs(A_pressure(j,j)*pressure_prime(i,j) + A_pressure(j,j+1)*pressure_prime(i,j+1) + A_pressure(j,j-1)*pressure_prime(i,j-1) - b_pressure[j]);
				   pressure_residual_sum[i_iter] = pressure_residual_sum[i_iter] + abs(b_pressure[j]);

				   /* We also solve pressure_prime */
				   pressure_prime[i][j] = (b_pressure[j] - A_pressure[j][j-1]*pressure_prime[i][j-1] - A_pressure[j][j+1]*pressure_prime[i][j+1] )/A_pressure[j][j];
			   }
		   }

		   if(i_iter == 0){
			   pressure_residual_sum_norm = pressure_residual_sum[0];

			   /* to avoid dividing by zero */
			   if (pressure_residual_sum_norm == 0.0){
				   pressure_residual_sum_norm = 1.0E4/residual_threshold;
			   }
		   }

		   pressure_residual_sum[i_iter] = pressure_residual_sum[i_iter]/pressure_residual_sum_norm;
	}

