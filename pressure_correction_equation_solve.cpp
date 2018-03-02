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
#include <Eigen/SVD>
#include <Eigen/Dense>
#include "globals.h"

using namespace std;

void pressure_correction_equation_solve(MatrixXd u_star, MatrixXd v_star, MatrixXd &pressure_prime, int i_iter)
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
	   	   MatrixXd A_pressure = MatrixXd::Zero(Nx,Nx);
	   	   VectorXd b_pressure = VectorXd::Zero(Nx);

		   /* For pressure nodes we assume density and Gamma=1/Reynolds_num constant */

		   double a_w;
		   double a_e;
		   double a_s;
		   double a_n;
		   double S_pr;
		   double a_p;

		   for(int i=ngcx; i<(Nodesx+ngcx); i++){
			   A_pressure = MatrixXd::Zero(Ny,Ny);
			   b_pressure = VectorXd::Zero(Ny);

			   for(int j=ngcy;j<(Nodesy+ngcx);j++){
				   a_w = density * d_u[i][j] * Area_velocity_node_u[i][j];;
				   a_e = density * d_u[i+1][j] * Area_velocity_node_u[i+1][j];
				   a_s = density * d_v[i][j] * Area_velocity_node_v[i][j];
				   a_n = density * d_v[i][j+1] * Area_velocity_node_v[i][j+1];
				   a_p = a_w + a_e + a_s + a_n;
				   S_pr = (density * u_star(i,j) * Area_velocity_node_u[i][j]) - (density * u_star((i+1),j) * Area_velocity_node_u[i+1][j]) + (density * v_star(i,j) * Area_velocity_node_v[i][j]) - (density * v_star(i,j+1) * Area_velocity_node_v[i][j+1]);
				   A_pressure(j,j) = a_p;
				   A_pressure(j,j+1) = - a_n;
				   A_pressure(j,j-1) = - a_s;
				   b_pressure(j) = a_e * pressure_prime((i+1),j) + a_w * pressure_prime((i-1),j) + S_pr;
//				   pressure_residual_sum[i_iter] = pressure_residual_sum[i_iter] + abs(A_pressure(j,j)*pressure_prime(i,j) + A_pressure(j,j+1)*pressure_prime(i,j+1) + A_pressure(j,j-1)*pressure_prime(i,j-1) - b_pressure[j]);
				   pressure_residual_sum[i_iter] = pressure_residual_sum[i_iter] + abs(b_pressure[j]);
			   }

//			   				      	      cout << "Pressure matrix A:" << endl;
//			   				      	      cout << A_pressure << endl;
//			   				      	      cout << endl;
//			   				      	      cout << "Pressure source vector b:" << endl;
//			   				      	      cout << b_pressure << endl;
			   /* We also solve pressure_prime */
			   pressure_prime.row(i) = A_pressure.colPivHouseholderQr().solve(b_pressure);
//			   pressure_prime.row(i) = A_pressure.ldlt().solve(b_pressure);
		   }

		   if(i_iter == 0){
			   pressure_residual_sum_norm = pressure_residual_sum[0];
		   }

		   pressure_residual_sum[i_iter] = pressure_residual_sum[i_iter]/pressure_residual_sum_norm;

//		   /* Guard cells not just at the border */
//
//		   	/* west guard cells */
//		   	for (int i=0;i<(ngcx - 1);i++){
//		   		for (int j=0;j<Ny;j++){
//		   			pressure[i][j] = pressure[ngcx][j];
//		   		}
//
//		   	}
//
//		   	/* east guard cells */
//		   	for (int i=(Nodesx + ngcx);i<Nx;i++){
//		   		for (int j=0;j<Ny;j++){
//		   			pressure[i][j] = pressure[(Nodesx + ngcx - 1)][j];
//		   		}
//
//		   	}
//
//		   	/* south guard cells */
//		   	for (int i=0;i<Nx;i++){
//		   		for (int j=0;j<(ngcy - 1);j++){
//		   			pressure[i][j] = pressure[i][ngcy];
//		   		}
//
//		   	}
//
//		   	/* north guard cells */
//		   	for (int i=0;i<Nx;i++){
//		   		for (int j=(Nodesy + ngcy);j<Ny;j++){
//		   			pressure[i][j] = pressure[i][(Nodesy + ngcy - 1)];
//		   		}
//
//		   	}

		   	/* Guard cells just at the border */
		   	/* We set the guard cells adjacent to nodes that are not exactly on the boundary, equal to the  */
		   	/* negative velocity, so when they sum, it is equal to zero, thus giving zero in between. */


//		   	for (int j=ngcy;j<(Nodesy + ngcy);j++){
//		   		/* west guard cells */
//		   		pressure[(ngcx - 1)][j] = pressure[ngcx][j];
//
//		   		/* east guard cells */
//		   		pressure[(Nodesx + ngcx)][j] = pressure[(Nodesx + ngcx - 1)][j];
//		   	}
//
//		   	for (int i=ngcx;i<(Nodesx + ngcx);i++) {
//		   		/* south guard cells */
//		   		pressure[i][(ngcy - 1)] = pressure[i][ngcy];
//
//		   		/* north guard cells */
//		   		pressure[i][(Nodesy + ngcy)] = pressure[i][(Nodesy + ngcy - 1)];
//		   	}


//		   /* BC */
//		   /* west guard cells */
//		   for (int i=0;i<ngcx;i++){
//			   for (int j=0;j<Ny;j++){
//				   pressure[i][j] = pressure[ngcx][j];
//			   }
//
//		   }
//		   /* east guard cells */
//		   for (int i=(Nodesx + ngcx - 1);i<Nx;i++){
//			   for (int j=0;j<Ny;j++){
//				   pressure[i][j] = pressure[(Nodesx + ngcx - 2)][j];
//			   }
//		   }
//
//		   /* south guard cells */
//		   for (int i=0;i<Nx;i++){
//			   for (int j=0;j<ngcy;j++){
//				   pressure[i][j] = pressure[i][ngcy];
//			   }
//		   }
//
//		   /* north guard cells */
//		   for (int i=0;i<Nx;i++){
//			   for (int j=(Nodesy + ngcy - 1);j<Ny;j++){
//				   pressure[i][j] = pressure[i][(Nodesy + ngcy - 2)];
//			   }
//		   }
	}

