/*
 * momentum_equation_solve.cpp
 *
 * This module solves the momentum equation and sets the velocity matrix values.
 * It takes the vector u* and modifies it.
 *
 *  Created on: Jul 13, 2017
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

void momentum_equation_solve(MatrixXd &u_star, MatrixXd &v_star, int i_iter)
{
	   /* Solve the discretised momentum equation using the guessed values of
		    * velocity and pressure. The equations for u and v velocities have the form:
		    *
		    * a_i,J * u_i,J^(*) = a_i-1,J * u_i-1,J^(*) + a_i+1,J * u_i+1,J^(*) +
		    *                     + a_i,J-1 * u_i,J-1^(*) + a_i,J+1 * u_i,J+1^(*) + pressure term + Source_u
		    *
		    * where u is the velocity in the x direction in the position i,J.
		    *
		    * a_I,j * v_I,j^(*) = a_I+1,j * v_I+1,j^(*) + a_I-1,j * v_I-1,j^(*) +
		    *                     + a_I,j+1 * v_I,j+1^(*) + a_I,j-1 * v_I,j-1^(*) + pressure term + Source_v
		    *
		    * For the boundary conditions we set all the border velocities (except the top border, or lid) to zero.
		    * The lid velocity (or top border) we set equal to lid_velocity from globals.
		    *
		    * */

		   vector<double> x_momentum_residual_sum_prev;
		   vector<double> y_momentum_residual_sum_prev;

		   x_momentum_residual_sum_prev.resize(MAX_ITER);
		   y_momentum_residual_sum_prev.resize(MAX_ITER);

		   /* initialize values of the previous momentum residual */
		   for (int i=0;i<(MAX_ITER);i++) {
			  x_momentum_residual_sum_prev[i] = x_momentum_residual_sum[i];
			  y_momentum_residual_sum_prev[i] = y_momentum_residual_sum[i];
		   }

		   /* these are the matrices A for u and v velocities that store the a_p, a_w, a_e, a_n and a_s
		    * and the vector b that stores the S_u and S_v values to solve the system of equations for
		    * the velocities. The system of equations is given by the momentum
		    * equations. We solve the system A_u*u=b_u, A_v*v=b_v and the solutions of u and v are set to
		    * be the guessed velocities u=u* and v=v*
		    */
		   MatrixXd A_u_velocity = MatrixXd::Zero(Ny,Ny);
		   VectorXd b_u_velocity = VectorXd::Zero(Ny);
		   MatrixXd A_v_velocity = MatrixXd::Zero(Ny,Ny);
		   VectorXd b_v_velocity = VectorXd::Zero(Ny);

		   /* For velocity nodes we assume density and Gamma=1/Reynolds_num constant */

		   double F_w;
		   double F_e;
		   double F_s;
		   double F_n;
		   double D_w;
		   double D_e;
		   double D_s;
		   double D_n;
		   double a_w;
		   double a_e;
		   double a_s;
		   double a_n;
		   double S_p;
		   double S_u;
		   double a_p;
		   double delta_F;
		   double Gamma_constant = 1/Reynolds_num; /* constant diffusion coefficient */
		   double delta_x = Lx/(Nx-0.5); /* spacing of the grid in x (uniform grid) */
		   double delta_y = Ly/(Ny-0.5); /* spacing of the grid in y (uniform grid) */

		   for(int i=0; i<Nx; i++){
			   A_u_velocity = MatrixXd::Zero(Ny,Ny);
			   b_u_velocity = VectorXd::Zero(Ny);
			   A_v_velocity = MatrixXd::Zero(Ny,Ny);
			   b_v_velocity = VectorXd::Zero(Ny);
			   for(int j=0;j<Ny;j++){

				   /* The values of the variables depend on whether they are at boundaries or not
				    * so we separate for each cases. First we write the nodes at boundaries */

				   if(i==0){
					   if(j==0){
						   F_w = 0.0;
						   F_e = 0.5 * density * ( u_velocity_old[i+1][j] + 0.5 * u_velocity_old[i][j] );
						   F_s = 0.25 * density * v_velocity_old[i][j];
						   F_n = 0.5 * density * v_velocity_old[i][j+1];
						   D_w = 0.0;
						   D_e = Gamma_constant/delta_x;
						   D_s = 0.25 * Gamma_constant/delta_y;
						   D_n = 0.5 * Gamma_constant/delta_y;
						   a_w = 0.0;
						   a_e = D_e - 0.5 * F_e;
						   a_s = 0.0;
						   a_n = D_n - 0.5 * F_n;
						   S_p = - ( D_w + F_w + D_s + F_s);
						   S_u = 0.0;
						   delta_F = F_e + F_n - F_s - F_w;
						   a_p = a_w + a_e + a_s + a_n + delta_F - S_p;
						   A_u_velocity(j,j) = a_p;
						   A_u_velocity(j,j+1) = - a_n;
						   b_u_velocity(j) = a_e * u_velocity_old[i+1][j] + S_u;
					   }
					   else if (j==(Ny-1)){
						   F_w = 0.0;
						   F_e = 0.5 * density * ( u_velocity_old[i+1][j] + 0.5 * u_velocity_old[i][j] );
						   F_s = 0.5 * density * v_velocity_old[i][j];
						   F_n = 0.0;
						   D_w = 0.0;
						   D_e = Gamma_constant/delta_x;
						   D_s = 0.5 * Gamma_constant/delta_y;
						   D_n = 0.25 * Gamma_constant/delta_y;
						   a_w = 0.0;
						   a_e = D_e - 0.5 * F_e;
						   a_s = D_s + 0.5 * F_s;
						   a_n = 0.0;
						   S_p = - ( D_w + F_w + D_n - F_n );
						   S_u = ( D_n - F_n ) * lid_velocity;
						   delta_F = F_e + F_n - F_s - F_w;
						   a_p = a_w + a_e + a_s + a_n + delta_F - S_p;
						   A_u_velocity(j,j) = a_p;
						   A_u_velocity(j,j-1) = - a_s;
						   b_u_velocity(j) = a_e * u_velocity_old[i+1][j] + S_u;
					   }
					   else {
						   F_w = 0.0;
				           F_e = 0.5 * density * ( u_velocity_old[i+1][j] + 0.5 * u_velocity_old[i][j] );
				           F_s = 0.5 * density * v_velocity_old[i][j];
				           F_n = 0.5 * density * v_velocity_old[i][j+1];
				           D_w = 0.0;
				           D_e = Gamma_constant/delta_x;
				           D_s = 0.5 * Gamma_constant/delta_y;
				           D_n = 0.5 * Gamma_constant/delta_y;
				           a_w = 0.0;
				           a_e = D_e - 0.5 * F_e;
				           a_s = D_s + 0.5 * F_s;
				           a_n = D_n - 0.5 * F_n;
				           S_p = - ( D_w + F_w );
				           S_u = 0.0;
				           delta_F = F_e + F_n - F_s - F_w;
				           a_p = a_w + a_e + a_s + a_n + delta_F - S_p;
				           A_u_velocity(j,j) = a_p;
				           A_u_velocity(j,j+1) = - a_n;
				           A_u_velocity(j,j-1) = - a_s;
				           b_u_velocity(j) = a_e * u_velocity_old[i+1][j] + S_u;
					   }
//					   u_star.row(i) = A_u_velocity.colPivHouseholderQr().solve(b_u_velocity);
					   /* We set the boundary conditions for the north and west wall */
					   u_star.row(i) = VectorXd::Zero(Ny);
					   u_star(i,(Ny-1)) = lid_velocity;
				   }
				   else if (i==1){
					   if(j==0){
						   F_w = 0.5 * density * ( u_velocity_old[i][j] + 0.5 * u_velocity_old[i-1][j]);
				           F_e = 0.5 * density * ( u_velocity_old[i+1][j] + u_velocity_old[i][j] );
				           F_s = 0.25 * density * (v_velocity_old[i][j] + v_velocity_old[i-1][j]);
				           F_n = 0.5 * density * (v_velocity_old[i][j+1] + v_velocity_old[i-1][j+1]);
				           D_w = Gamma_constant/delta_x;
				           D_e = Gamma_constant/delta_x;
				           D_s = 0.5 * Gamma_constant/delta_y;
				           D_n = Gamma_constant/delta_y;
				           a_w = D_w + 0.5 * F_w;
				           a_e = D_e - 0.5 * F_e;
				           a_s = 0.0;
				           a_n = D_n - 0.5 * F_n;
				           S_p = - ( D_s + F_s );
				           S_u = 0.0;
				           delta_F = F_e + F_n - F_s - F_w;
				           a_p = a_w + a_e + a_s + a_n + delta_F - S_p;
				           A_u_velocity(j,j) = a_p;
				           A_u_velocity(j,j+1) = - a_n;
				           b_u_velocity(j) = a_e * u_velocity_old[i+1][j] + a_w * u_velocity_old[i-1][j] + ( pressure_old[i-1][j] - pressure_old[i][j] ) * Area_velocity_node_u[i][j] + S_u;
					   }
					   else if (j==(Ny-1)){
						   F_w = 0.5 * density * ( u_velocity_old[i][j] + 0.5 * u_velocity_old[i-1][j]);
				           F_e = 0.5 * density * ( u_velocity_old[i+1][j] + u_velocity_old[i][j] );
				           F_s = 0.5 * density * (v_velocity_old[i][j] + v_velocity_old[i-1][j]);
				           F_n = 0.0;
				           D_w = Gamma_constant/delta_x;
				           D_e = Gamma_constant/delta_x;
				           D_s = Gamma_constant/delta_y;
				           D_n = 0.5 * Gamma_constant/delta_y;
				           a_w = D_w + 0.5 * F_w;
				           a_e = D_e - 0.5 * F_e;
				           a_s = D_s + 0.5 * F_s;
				           a_n = 0.0;
				           S_p = - ( D_n - F_n );
				           S_u = ( D_n - F_n ) * lid_velocity;
				           delta_F = F_e + F_n - F_s - F_w;
				           a_p = a_w + a_e + a_s + a_n + delta_F - S_p;
				           A_u_velocity(j,j) = a_p;
				           A_u_velocity(j,j-1) = - a_s;
				           b_u_velocity(j) = a_e * u_velocity_old[i+1][j] + a_w * u_velocity_old[i-1][j] + ( pressure_old[i-1][j] - pressure_old[i][j] ) * Area_velocity_node_u[i][j] + S_u;
					   }
					   else {
						   F_w = 0.5 * density * ( u_velocity_old[i][j] + 0.5 * u_velocity_old[i-1][j]);
				           F_e = 0.5 * density * ( u_velocity_old[i+1][j] + u_velocity_old[i][j] );
				           F_s = 0.5 * density * (v_velocity_old[i][j] + v_velocity_old[i-1][j]);
				           F_n = 0.5 * density * (v_velocity_old[i][j+1] + v_velocity_old[i-1][j+1]);
				           D_w = Gamma_constant/delta_x;
				           D_e = Gamma_constant/delta_x;
				           D_s = Gamma_constant/delta_y;
				           D_n = Gamma_constant/delta_y;
				           a_w = D_w + 0.5 * F_w;
				           a_e = D_e - 0.5 * F_e;
				           a_s = D_s + 0.5 * F_s;
				           a_n = D_n - 0.5 * F_n;
				           S_p = 0.0;
				           S_u = 0.0;
				           delta_F = F_e + F_n - F_s - F_w;
				           a_p = a_w + a_e + a_s + a_n + delta_F - S_p;
				           A_u_velocity(j,j) = a_p;
				           A_u_velocity(j,j+1) = - a_n;
				           A_u_velocity(j,j-1) = - a_s;
				           b_u_velocity(j) = a_e * u_velocity_old[i+1][j] + a_w * u_velocity_old[i-1][j] + ( pressure_old[i-1][j] - pressure_old[i][j] ) * Area_velocity_node_u[i][j] +  S_u;
					   }
					   u_star.row(i) = A_u_velocity.colPivHouseholderQr().solve(b_u_velocity);
					   u_star(i,(Ny-1)) = lid_velocity;
				   }
				   else if (i==(Nx-1)){
					   if(j==0){
						   F_w = 0.5 * density * ( u_velocity_old[i][j] + u_velocity_old[i-1][j]);
				           F_e = 0.5 * density * u_velocity_old[i][j];
				           F_s = 0.25 * density * (v_velocity_old[i][j] + v_velocity_old[i-1][j]);
				           F_n = 0.5 * density * (v_velocity_old[i][j+1] + v_velocity_old[i-1][j+1]);
				           D_w = Gamma_constant/delta_x;
				           D_e = Gamma_constant/delta_x;
				           D_s = 0.5 * Gamma_constant/delta_y;
				           D_n = Gamma_constant/delta_y;
				           a_w = D_w + 0.5 * F_w;
				           a_e = 0.0;
				           a_s = 0.0;
				           a_n = D_n - 0.5 * F_n;
				           S_p = - ( D_e + D_s - F_e + F_s );
				           S_u = 0.0;
				           delta_F = F_e + F_n - F_s - F_w;
				           a_p = a_w + a_e + a_s + a_n + delta_F - S_p;
				           A_u_velocity(j,j) = a_p;
				           A_u_velocity(j,j+1) = - a_n;
				           b_u_velocity(j) = a_w * u_velocity_old[i-1][j] + ( pressure_old[i-1][j] - pressure_old[i][j] ) * Area_velocity_node_u[i][j] +  S_u;
					   }
					   else if (j==(Ny-1)){
						   F_w = 0.5 * density * ( u_velocity_old[i][j] + u_velocity_old[i-1][j]);
				           F_e = 0.5 * density * u_velocity_old[i][j];
				           F_s = 0.5 * density * (v_velocity_old[i][j] + v_velocity_old[i-1][j]);
				           F_n = 0.0;
				           D_w = Gamma_constant/delta_x;
				           D_e = Gamma_constant/delta_x;
				           D_s = Gamma_constant/delta_y;
				           D_n = 0.5 * Gamma_constant/delta_y;
				           a_w = D_w + 0.5 * F_w;
				           a_e = 0.0;
				           a_s = D_s + 0.5 * F_s;
				           a_n = 0.0;
				           S_p = - ( D_e + D_n - F_e - F_n );
				           S_u = ( D_n - F_n ) * lid_velocity;
				           delta_F = F_e + F_n - F_s - F_w;
				           a_p = a_w + a_e + a_s + a_n + delta_F - S_p;
				           A_u_velocity(j,j) = a_p;
				           A_u_velocity(j,j-1) = - a_s;
				           b_u_velocity(j) = a_w * u_velocity_old[i-1][j] + ( pressure_old[i-1][j] - pressure_old[i][j] ) * Area_velocity_node_u[i][j] +  S_u;
					   }
					   else {
						   F_w = 0.5 * density * ( u_velocity_old[i][j] + u_velocity_old[i-1][j]);
				           F_e = 0.5 * density * u_velocity_old[i][j];
				           F_s = 0.5 * density * (v_velocity_old[i][j] + v_velocity_old[i-1][j]);
				           F_n = 0.5 * density * (v_velocity_old[i][j+1] + v_velocity_old[i-1][j+1]);
				           D_w = Gamma_constant/delta_x;
				           D_e = Gamma_constant/delta_x;
				           D_s = Gamma_constant/delta_y;
				           D_n = Gamma_constant/delta_y;
				           a_w = D_w + 0.5 * F_w;
				           a_e = 0.0;
				           a_s = D_s + 0.5 * F_s;
				           a_n = D_n - 0.5 * F_n;
				           S_p = - ( D_e - F_e );
				           S_u = 0.0;
				           delta_F = F_e + F_n - F_s - F_w;
				           a_p = a_w + a_e + a_s + a_n + delta_F - S_p;
				           A_u_velocity(j,j) = a_p;
				           A_u_velocity(j,j+1) = - a_n;
				           A_u_velocity(j,j-1) = - a_s;
				           b_u_velocity(j) = a_w * u_velocity_old[i-1][j] + ( pressure_old[i-1][j] - pressure_old[i][j] ) * Area_velocity_node_u[i][j] +  S_u;
					   }
					   u_star.row(i) = A_u_velocity.colPivHouseholderQr().solve(b_u_velocity);
					   u_star(i,(Ny-1)) = lid_velocity;
				   }
				   else {
					   if(j==0){
						   F_w = 0.5 * density * ( u_velocity_old[i][j] + u_velocity_old[i-1][j]);
				           F_e = 0.5 * density * ( u_velocity_old[i+1][j] + u_velocity_old[i][j]);
				           F_s = 0.25 * density * (v_velocity_old[i][j] + v_velocity_old[i-1][j]);
				           F_n = 0.5 * density * (v_velocity_old[i][j+1] + v_velocity_old[i-1][j+1]);
				           D_w = Gamma_constant/delta_x;
				           D_e = Gamma_constant/delta_x;
				           D_s = 0.5 * Gamma_constant/delta_y;
				           D_n = Gamma_constant/delta_y;
				           a_w = D_w + 0.5 * F_w;
				           a_e = D_e - 0.5 * F_e;
				           a_s = 0.0;
				           a_n = D_n - 0.5 * F_n;
				           S_p = - ( D_s + F_s );
				           S_u = 0.0;
				           delta_F = F_e + F_n - F_s - F_w;
				           a_p = a_w + a_e + a_s + a_n + delta_F - S_p;
				           A_u_velocity(j,j) = a_p;
				           A_u_velocity(j,j+1) = - a_n;
				           b_u_velocity(j) = a_e * u_velocity_old[i+1][j] + a_w * u_velocity_old[i-1][j] + ( pressure_old[i-1][j] - pressure_old[i][j] ) * Area_velocity_node_u[i][j] +  S_u;
					   }
					   else if (j==(Ny-1)){
						   F_w = 0.5 * density * ( u_velocity_old[i][j] + u_velocity_old[i-1][j]);
				           F_e = 0.5 * density * ( u_velocity_old[i+1][j] + u_velocity_old[i][j]);
				           F_s = 0.5 * density * (v_velocity_old[i][j] + v_velocity_old[i-1][j]);
				           F_n = 0.0;
				           D_w = Gamma_constant/delta_x;
				           D_e = Gamma_constant/delta_x;
				           D_s = Gamma_constant/delta_y;
				           D_n = 0.5 * Gamma_constant/delta_y;
				           a_w = D_w + 0.5 * F_w;
				           a_e = D_e - 0.5 * F_e;
				           a_s = D_s + 0.5 * F_s;
				           a_n = 0.0;
				           S_p = - ( D_n - F_n );
				           S_u = ( D_n - F_n ) * lid_velocity;
				           delta_F = F_e + F_n - F_s - F_w;
				           a_p = a_w + a_e + a_s + a_n + delta_F - S_p;
				           A_u_velocity(j,j) = a_p;
				           A_u_velocity(j,j-1) = - a_s;
				           b_u_velocity(j) = a_e * u_velocity_old[i+1][j] + a_w * u_velocity_old[i-1][j] + ( pressure_old[i-1][j] - pressure_old[i][j] ) * Area_velocity_node_u[i][j] +  S_u;
					   }
					   else {
						   F_w = 0.5 * density * ( u_velocity_old[i][j] + u_velocity_old[i-1][j]);
				           F_e = 0.5 * density * ( u_velocity_old[i+1][j] + u_velocity_old[i][j]);
				           F_s = 0.5 * density * (v_velocity_old[i][j] + v_velocity_old[i-1][j]);
				           F_n = 0.5 * density * (v_velocity_old[i][j+1] + v_velocity_old[i-1][j+1]);
				           D_w = Gamma_constant/delta_x;
				           D_e = Gamma_constant/delta_x;
				           D_s = Gamma_constant/delta_y;
				           D_n = Gamma_constant/delta_y;
				           a_w = D_w + 0.5 * F_w;
				           a_e = D_e - 0.5 * F_e;
				           a_s = D_s + 0.5 * F_s;
				           a_n = D_n - 0.5 * F_n;
				           S_p = 0.0;
				           S_u = 0.0;
				           delta_F = F_e + F_n - F_s - F_w;
				           a_p = a_w + a_e + a_s + a_n + delta_F - S_p;
				           A_u_velocity(j,j) = a_p;
				           A_u_velocity(j,j+1) = - a_n;
				           A_u_velocity(j,j-1) = - a_s;
				           b_u_velocity(j) = a_e * u_velocity_old[i+1][j] + a_w * u_velocity_old[i-1][j] + ( pressure_old[i-1][j] - pressure_old[i][j] ) * Area_velocity_node_u[i][j] +  S_u;
					   }
					   u_star.row(i) = A_u_velocity.colPivHouseholderQr().solve(b_u_velocity);
					   u_star(i,(Ny-1)) = lid_velocity;
				   }
			   }
		   }


//	   /* Solve the discretised momentum equation using the guessed values of
//	    * velocity and pressure. The equation has the form:
//	    *
//	    * a_P * u_P^(*) = a_W * u_W^(*) + a_E * u_E^(*) + S_u
//	    *
//	    * where u is the velocity in the x-direction.
//	    *
//	    * */
//
//	   /* First we determine the values for the boundary nodes */
//
//
//	//   /* For velocity node 1
//	//    *
//	//    * we use the upwind scheme and the
//	//    * deferred correction approach, where the negative contribution to ap
//	//    * is placed on the RHS and we use u1_old as the nodal velocity in the
//	//    * previous iteration */
//	//
//	//   double F_e_x1;
//	//   double F_w_x1;
//	//   double a_p_x1;
//	//   double Su_x1;
//	//
//	//   F_e_x1 = density * Area_pressure_node[1] * ( velocity[0] + velocity[1]) * 0.5;
//	//   F_w_x1 = density * Area_velocity_node[0] * velocity[0];
//	//
//	//   a_p_x1 = (F_e_x1 + F_w_x1 * 0.5 * (Area_velocity_node[0]/Area_pressure_node[0]) * (Area_velocity_node[0]/Area_pressure_node[0]))/urfu;
//	//   Su_x1 = (p_in - pressure[1]) * Area_velocity_node[0] + F_w_x1 * (Area_velocity_node[0]/Area_pressure_node[0]) * velocity_old[0] + (1.0 - urfu)*a_p_x1*velocity_old[0];
//	//   dx[0] = Area_velocity_node[0]/a_p_x1;
//
//
//	   /* For velocity node 1
//	     *
//	     * we use only the upwind scheme */
//	    double F_e_x1;
//	    double F_w_x1;
//	    double a_p_x1;
//	    double Su_x1;
//
//	    F_e_x1 = density * Area_pressure_node[1] * ( velocity[0] + velocity[1]) * 0.5;
//	    //F_w_x1 = density * Area_velocity_node[0] * velocity[0]; /* BC suggested in Versteeg book */
//	    F_w_x1 = density * Area_pressure_node[0] * (velocity[0] - (velocity[1] - velocity[0]) * 0.5); /* BC extrapolated from tendency of previous node */
//
//	    a_p_x1 = (F_e_x1 - F_w_x1 * (Area_velocity_node[0]/Area_pressure_node[0]) + F_w_x1 * 0.5 * (Area_velocity_node[0]/Area_pressure_node[0]) * (Area_velocity_node[0]/Area_pressure_node[0]) )/urfu;
//	    Su_x1 = (p_in - pressure[1]) * Area_velocity_node[0] + (1.0 - urfu)*a_p_x1*velocity_old[0];
//	    dx[0] = Area_velocity_node[0]/a_p_x1;
//
//
//	   /* For velocity node Nx-1  */
//	   double F_e_xNm1;
//	   double F_w_xNm1;
//	   double a_p_xNm1;
//	   double a_w_xNm1;
//	   double Su_xNm1;
//
//	   //F_e_xNm1 = density * Area_velocity_node[Nx-2] *  velocity[Nx-2]; /* BC suggested in Versteeg book */
//	   F_e_xNm1 = density * Area_pressure_node[Nx-1] *  (velocity[Nx-2] - (velocity[Nx-3] - velocity[Nx-2]) * 0.5); /* BC extrapolated from tendency of previous node */
//	   F_w_xNm1 = density * Area_pressure_node[Nx-2] * ( velocity[Nx-3] + velocity[Nx-2]) * 0.5;
//	   a_p_xNm1 = F_e_xNm1/urfu;
//	   a_w_xNm1 = F_w_xNm1;
//	   Su_xNm1 = (pressure[Nx-2] - pressure[Nx-1]) * Area_velocity_node[Nx-2] + (1.0 - urfu)*a_p_xNm1*velocity_old[Nx-2];
//	   dx[Nx-2] = Area_velocity_node[Nx-2]/a_p_xNm1;
//
//	   /* this is the matrix A that stores the a_p, a_w and a_e and the vector
//	    * b that stores the S_u values to solve the system of equations for
//	    * the velocities. The system of equations is given by the momentum
//	    * equations. We solve the system A*u=b and the solutions of u are set to
//	    * be the guessed velocities u=u*
//	    */
//	   MatrixXd A_velocity = MatrixXd::Zero(Nx-1,Nx-1);
//	   VectorXd b_velocity = VectorXd::Zero(Nx-1);
//
//	   /* Setting the boundary nodes coefficients */
//	   A_velocity(0,0) = a_p_x1;
//	   b_velocity(0) = Su_x1;
//
//	   A_velocity(Nx-2,Nx-2) = a_p_xNm1;
//	   A_velocity(Nx-2,Nx-3) = -a_w_xNm1;
//	   b_velocity(Nx-2) = Su_xNm1;
//
//	   /* adding the values of the residuals of the momentum equations
//	    * of boundaries */
//	   x_momentum_residual_sum[i_iter] = abs(A_velocity(0,0)*velocity[0] + A_velocity(0,1)*velocity[1] - b_velocity[0]);
//	   x_momentum_residual_sum[i_iter] = x_momentum_residual_sum[i_iter] + abs(A_velocity(Nx-2,Nx-2)*velocity[Nx-2] + A_velocity(Nx-2,Nx-3)*velocity[Nx-3] - b_velocity[Nx-2]);
//
//	   /* For velocity nodes in between  */
//	   /* i=1..Nx-3 */
//	   double F_e_x;
//	   double F_w_x;
//	   double a_p_x;
//	   double a_w_x;
//	   double Su_x;
//
//	   for(int i=1; i<(Nx-2); i++){
//
//		  F_w_x = density * Area_pressure_node[i] * (velocity[i-1] + velocity[i]) * 0.5;
//		  F_e_x = density * Area_pressure_node[i+1] * (velocity[i] + velocity[i+1]) * 0.5;
//		  a_p_x = F_e_x/urfu;
//		  a_w_x = F_w_x;
//		  Su_x = (pressure[i] - pressure[i+1]) * Area_velocity_node[i] + (1.0 - urfu)*a_p_x*velocity_old[i];
//		  dx[i] = Area_velocity_node[i]/a_p_x;
//
//		  A_velocity(i,i) = a_p_x;
//		  A_velocity(i,i-1) = -a_w_x;
//		  b_velocity(i) = Su_x;
//
//		  /* adding the values of the residuals of the momentum equations
//		   * for the nodes in between */
//		  x_momentum_residual_sum[i_iter] = x_momentum_residual_sum[i_iter] + abs(A_velocity(i,i)*velocity[i] + A_velocity(i,i+1)*velocity[i+1] + A_velocity(i,i-1)*velocity[i-1] - b_velocity[i]);
//
//
//	   }
//
//
//	  /* These are the guessed velocities (u_star*) obtained by solving the system of
//	   * equations. */
//	   u_star = A_velocity.colPivHouseholderQr().solve(b_velocity);
//
//	   //    	   if (i_iter == 0){
//	   //
//	   //    	   cout << "Momentum matrix A:" << endl;
//	   //    	   cout << A_velocity << endl;
//	   //    	   cout << endl;
//	   //    	   cout << "Momentum source vector b:" << endl;
//	   //    	   cout << b_velocity << endl;
//	   //
//	   //    	   }


		   	       	   cout << "Momentum matrix A:" << endl;
		   	       	   cout << A_u_velocity << endl;
		   	       	   cout << endl;
		   	       	   cout << "Momentum source vector b:" << endl;
		   	       	   cout << b_u_velocity << endl;

		   	       		   cout << u_star << endl;
		   	       		   cout << " " << endl;
		   	       		   cout << v_star << endl;


}




