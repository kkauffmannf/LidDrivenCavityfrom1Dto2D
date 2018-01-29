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

void momentum_equation_solve(VectorXd &u_star, int i_iter)
{
	   /* Solve the discretised momentum equation using the guessed values of
	    * velocity and pressure. The equation has the form:
	    *
	    * a_P * u_P^(*) = a_W * u_W^(*) + a_E * u_E^(*) + S_u
	    *
	    * where u is the velocity in the x-direction.
	    *
	    * */

	   /* First we determine the values for the boundary nodes */


	//   /* For velocity node 1
	//    *
	//    * we use the upwind scheme and the
	//    * deferred correction approach, where the negative contribution to ap
	//    * is placed on the RHS and we use u1_old as the nodal velocity in the
	//    * previous iteration */
	//
	//   double F_e_x1;
	//   double F_w_x1;
	//   double a_p_x1;
	//   double Su_x1;
	//
	//   F_e_x1 = density * Area_pressure_node[1] * ( velocity[0] + velocity[1]) * 0.5;
	//   F_w_x1 = density * Area_velocity_node[0] * velocity[0];
	//
	//   a_p_x1 = (F_e_x1 + F_w_x1 * 0.5 * (Area_velocity_node[0]/Area_pressure_node[0]) * (Area_velocity_node[0]/Area_pressure_node[0]))/urfu;
	//   Su_x1 = (p_in - pressure[1]) * Area_velocity_node[0] + F_w_x1 * (Area_velocity_node[0]/Area_pressure_node[0]) * velocity_old[0] + (1.0 - urfu)*a_p_x1*velocity_old[0];
	//   dx[0] = Area_velocity_node[0]/a_p_x1;


	   /* For velocity node 1
	     *
	     * we use only the upwind scheme */
	    double F_e_x1;
	    double F_w_x1;
	    double a_p_x1;
	    double Su_x1;

	    F_e_x1 = density * Area_pressure_node[1] * ( velocity[0] + velocity[1]) * 0.5;
	    //F_w_x1 = density * Area_velocity_node[0] * velocity[0]; /* BC suggested in Versteeg book */
	    F_w_x1 = density * Area_pressure_node[0] * (velocity[0] - (velocity[1] - velocity[0]) * 0.5); /* BC extrapolated from tendency of previous node */

	    a_p_x1 = (F_e_x1 - F_w_x1 * (Area_velocity_node[0]/Area_pressure_node[0]) + F_w_x1 * 0.5 * (Area_velocity_node[0]/Area_pressure_node[0]) * (Area_velocity_node[0]/Area_pressure_node[0]) )/urfu;
	    Su_x1 = (p_in - pressure[1]) * Area_velocity_node[0] + (1.0 - urfu)*a_p_x1*velocity_old[0];
	    dx[0] = Area_velocity_node[0]/a_p_x1;


	   /* For velocity node Nx-1  */
	   double F_e_xNm1;
	   double F_w_xNm1;
	   double a_p_xNm1;
	   double a_w_xNm1;
	   double Su_xNm1;

	   //F_e_xNm1 = density * Area_velocity_node[Nx-2] *  velocity[Nx-2]; /* BC suggested in Versteeg book */
	   F_e_xNm1 = density * Area_pressure_node[Nx-1] *  (velocity[Nx-2] - (velocity[Nx-3] - velocity[Nx-2]) * 0.5); /* BC extrapolated from tendency of previous node */
	   F_w_xNm1 = density * Area_pressure_node[Nx-2] * ( velocity[Nx-3] + velocity[Nx-2]) * 0.5;
	   a_p_xNm1 = F_e_xNm1/urfu;
	   a_w_xNm1 = F_w_xNm1;
	   Su_xNm1 = (pressure[Nx-2] - pressure[Nx-1]) * Area_velocity_node[Nx-2] + (1.0 - urfu)*a_p_xNm1*velocity_old[Nx-2];
	   dx[Nx-2] = Area_velocity_node[Nx-2]/a_p_xNm1;

	   /* this is the matrix A that stores the a_p, a_w and a_e and the vector
	    * b that stores the S_u values to solve the system of equations for
	    * the velocities. The system of equations is given by the momentum
	    * equations. We solve the system A*u=b and the solutions of u are set to
	    * be the guessed velocities u=u*
	    */
	   MatrixXd A_velocity = MatrixXd::Zero(Nx-1,Nx-1);
	   VectorXd b_velocity = VectorXd::Zero(Nx-1);

	   /* Setting the boundary nodes coefficients */
	   A_velocity(0,0) = a_p_x1;
	   b_velocity(0) = Su_x1;

	   A_velocity(Nx-2,Nx-2) = a_p_xNm1;
	   A_velocity(Nx-2,Nx-3) = -a_w_xNm1;
	   b_velocity(Nx-2) = Su_xNm1;

	   /* adding the values of the residuals of the momentum equations
	    * of boundaries */
	   x_momentum_residual_sum[i_iter] = abs(A_velocity(0,0)*velocity[0] + A_velocity(0,1)*velocity[1] - b_velocity[0]);
	   x_momentum_residual_sum[i_iter] = x_momentum_residual_sum[i_iter] + abs(A_velocity(Nx-2,Nx-2)*velocity[Nx-2] + A_velocity(Nx-2,Nx-3)*velocity[Nx-3] - b_velocity[Nx-2]);

	   /* For velocity nodes in between  */
	   /* i=1..Nx-3 */
	   double F_e_x;
	   double F_w_x;
	   double a_p_x;
	   double a_w_x;
	   double Su_x;

	   for(int i=1; i<(Nx-2); i++){

		  F_w_x = density * Area_pressure_node[i] * (velocity[i-1] + velocity[i]) * 0.5;
		  F_e_x = density * Area_pressure_node[i+1] * (velocity[i] + velocity[i+1]) * 0.5;
		  a_p_x = F_e_x/urfu;
		  a_w_x = F_w_x;
		  Su_x = (pressure[i] - pressure[i+1]) * Area_velocity_node[i] + (1.0 - urfu)*a_p_x*velocity_old[i];
		  dx[i] = Area_velocity_node[i]/a_p_x;

		  A_velocity(i,i) = a_p_x;
		  A_velocity(i,i-1) = -a_w_x;
		  b_velocity(i) = Su_x;

		  /* adding the values of the residuals of the momentum equations
		   * for the nodes in between */
		  x_momentum_residual_sum[i_iter] = x_momentum_residual_sum[i_iter] + abs(A_velocity(i,i)*velocity[i] + A_velocity(i,i+1)*velocity[i+1] + A_velocity(i,i-1)*velocity[i-1] - b_velocity[i]);


	   }


	  /* These are the guessed velocities (u_star*) obtained by solving the system of
	   * equations. */
	   u_star = A_velocity.colPivHouseholderQr().solve(b_velocity);

	   //    	   if (i_iter == 0){
	   //
	   //    	   cout << "Momentum matrix A:" << endl;
	   //    	   cout << A_velocity << endl;
	   //    	   cout << endl;
	   //    	   cout << "Momentum source vector b:" << endl;
	   //    	   cout << b_velocity << endl;
	   //
	   //    	   }

}




