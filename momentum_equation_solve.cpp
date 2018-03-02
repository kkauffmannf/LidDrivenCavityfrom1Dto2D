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

		   /* initialize values of the momentum residual */
			x_momentum_residual_sum[i_iter] = 0.0;
			y_momentum_residual_sum[i_iter] = 0.0;

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
		   double S_v;
		   double a_p;
		   double delta_F;
		   double Gamma_constant = 1/Reynolds_num; /* constant diffusion coefficient */
		   double delta_x = Lx/(Nx-0.5); /* spacing of the grid in x (uniform grid) */
		   double delta_y = Ly/(Ny-0.5); /* spacing of the grid in y (uniform grid) */

		   /* For each new i, we set the A matrix and b vector to zero */

		   for(int i=ngcx; i<(Nodesx+ngcx); i++){
			   A_u_velocity = MatrixXd::Zero(Ny,Ny);
			   b_u_velocity = VectorXd::Zero(Ny);
			   A_v_velocity = MatrixXd::Zero(Ny,Ny);
			   b_v_velocity = VectorXd::Zero(Ny);

			   for(int j=ngcy;j<(Nodesy+ngcy);j++){

				   /* The values of the variables depend on whether they are at boundaries (north, south, east, west) or not
				    * so we separate for each cases. First we write the nodes at boundaries */
				   /* Equations for u_velocity */
				   F_w = 0.5 * density * ( u_velocity[i][j] + u_velocity[i-1][j]);
				   F_e = 0.5 * density * ( u_velocity[i+1][j] + u_velocity[i][j]);
				   F_s = 0.5 * density * (v_velocity[i][j] + v_velocity[i-1][j]);
				   F_n = 0.5 * density * (v_velocity[i][j+1] + v_velocity[i-1][j+1]);
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
				   a_p = a_p/urfu;
				   A_u_velocity(j,j) = a_p;
				   A_u_velocity(j,j+1) = - a_n;
				   A_u_velocity(j,j-1) = - a_s;
				   b_u_velocity(j) = a_e * u_velocity[i+1][j] + a_w * u_velocity[i-1][j] + ( pressure[i-1][j] - pressure[i][j] ) * Area_velocity_node_u[i][j] +  S_u + (1 - urfu) * a_p * u_velocity_old[i][j];
				   d_u[i][j] = Area_velocity_node_u[i][j]/a_p;

				   /* Equations for v_velocity */
				   F_w = 0.5 * density * ( u_velocity[i][j] + u_velocity[i][j-1]);
				   F_e = 0.5 * density * ( u_velocity[i+1][j] + u_velocity[i+1][j-1]);
				   F_s = 0.5 * density * (v_velocity[i][j-1] + v_velocity[i][j]);
				   F_n = 0.5 * density * (v_velocity[i][j] + v_velocity[i][j+1]);
				   D_w = Gamma_constant/delta_x;
				   D_e = Gamma_constant/delta_x;
				   D_s = Gamma_constant/delta_y;
				   D_n = Gamma_constant/delta_y;
				   a_w = D_w + 0.5 * F_w;
				   a_e = D_e - 0.5 * F_e;
				   a_s = D_s + 0.5 * F_s;
				   a_n = D_n - 0.5 * F_n;
				   S_p = 0.0;
				   S_v = 0.0;
				   delta_F = F_e + F_n - F_s - F_w;
				   a_p = a_w + a_e + a_s + a_n + delta_F - S_p;
				   a_p = a_p/urfv;
				   A_v_velocity(j,j) = a_p;
				   A_v_velocity(j,j+1) = - a_n;
				   A_v_velocity(j,j-1) = - a_s;
				   b_v_velocity(j) = a_e * v_velocity[i+1][j] + a_w * v_velocity[i-1][j] + ( pressure[i][j-1] - pressure[i][j] ) * Area_velocity_node_v[i][j] +  S_v + (1 - urfv) * a_p * v_velocity_old[i][j];
				   d_v[i][j] = Area_velocity_node_v[i][j]/a_p;

				   /* adding the values of the residuals of the momentum equations */
				   x_momentum_residual_sum[i_iter] = x_momentum_residual_sum[i_iter] + abs(A_u_velocity(j,j)*u_velocity[i][j] + A_u_velocity(j,j+1)*u_velocity[i][j+1] + A_u_velocity(j,j-1)*u_velocity[i][j-1] - b_u_velocity[j]);
				   y_momentum_residual_sum[i_iter] = y_momentum_residual_sum[i_iter] + abs(A_v_velocity(j,j)*v_velocity[i][j] + A_v_velocity(j,j+1)*v_velocity[i][j+1] + A_v_velocity(j,j-1)*v_velocity[i][j-1] - b_v_velocity[j]);
			   }

			   u_star.row(i) = A_u_velocity.colPivHouseholderQr().solve(b_u_velocity);
			   v_star.row(i) = A_v_velocity.colPivHouseholderQr().solve(b_v_velocity);

//			   u_star.row(i) = A_u_velocity.ldlt().solve(b_u_velocity);
//			   v_star.row(i) = A_v_velocity.ldlt().solve(b_v_velocity);
		   }

		   if(i_iter == 0){
			   x_momentum_residual_sum_norm = x_momentum_residual_sum[0];
			   y_momentum_residual_sum_norm = y_momentum_residual_sum[0];

			   /* to avoid dividing by zero */
			   if (x_momentum_residual_sum_norm == 0.0){
				   x_momentum_residual_sum_norm = 1.0E4/residual_threshold;
			   }
			   if (y_momentum_residual_sum_norm == 0.0){
				   y_momentum_residual_sum_norm = 1.0E4/residual_threshold;
			   }
		   }

		   x_momentum_residual_sum[i_iter] = x_momentum_residual_sum[i_iter]/x_momentum_residual_sum_norm;
		   y_momentum_residual_sum[i_iter] = y_momentum_residual_sum[i_iter]/y_momentum_residual_sum_norm;

			/* BC. We set the guard cells adjacent to nodes that are not exactly on the boundary, equal to the  */
			/* negative velocity, so when they sum, it is equal to zero */
			for (int j=ngcy;j<(Nodesy + ngcy);j++){
				u_star((Nodesx + ngcx),j) = -1.0 * u_star((Nodesx + ngcx - 1),j);
				v_star((ngcx - 1),j) = -1.0 * v_star(ngcx,j);
				u_star(ngcx,j) = 0.0;
				v_star((Nodesx + ngcx - 1),j) = 0.0;
			}

			for (int i=ngcx;i<(Nodesx + ngcx);i++) {
				u_star(i,(ngcy - 1)) = -1.0 * u_star(i,ngcy);
				v_star(i,(Nodesy + ngcy)) = - 1.0 * v_star(i,(Nodesy + ngcy - 1));
				u_star(i,(Nodesy + ngcy - 1)) = lid_velocity;
				v_star(i,ngcy) = 0.0;
			}

		   /* BC */
			/* west guard cells */
			for (int i=0;i<ngcx;i++){
				for (int j=0;j<Ny;j++){
					u_star(i,j) = 0.0;
					v_star(i,j) = 0.0;
				}

			}

			/* east guard cells */
			for (int i=(Nodesx + ngcx);i<Nx;i++){
				for (int j=0;j<Ny;j++){
					u_star(i,j) = 0.0;
					v_star(i,j) = 0.0;
				}

			}

			/* south guard cells */
			for (int i=0;i<Nx;i++){
				for (int j=0;j<ngcy;j++){
					u_star(i,j) = 0.0;
					v_star(i,j) = 0.0;
				}

			}

			/* north guard cells */
			for (int i=0;i<Nx;i++){
				for (int j=(Nodesy + ngcy);j<Ny;j++){
					u_star(i,j) = 0.0;
					v_star(i,j) = 0.0;
				}

			}

//			u_star(ngcx,ngcy) = 0.0;
			v_star(ngcx,ngcy) = 0.0;
			u_star((Nodesx + ngcx - 1),(Nodesy + ngcy - 1)) = lid_velocity;
//			v_star((Nodesx + ngcx - 1),(Nodesy + ngcy - 1)) = 0.0;
}




