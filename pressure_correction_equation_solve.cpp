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

void pressure_correction_equation_solve(VectorXd u_star, VectorXd &pressure_prime, int i_iter)
{

	   /* Solve the discretised pressure correction equation using the guessed values of
	    * velocity and the values of pressure. The equation has the form:
	    *
	    * a_P * p'_P = a_W * p'_W + a_E * p'_E + b'
	    *
	    * where p' is the pressure corrections p'=(p-p*) in each node.
	    *
	    * */


	   /* First we determine the values for the boundary nodes.
	    * For the first and N-1 pressure nodes the value of the
	    * correction is set to 0.0 (see reference Versteeg) */

	   /* For pressure node A
	    *
	    * the equation is p'1=0.0 */
	   double a_p_p1;
	   double b_p1;

	   a_p_p1 = 1.0;
	   b_p1 = 0.0;

	  /* For pressure node N
	   *
	   * the equation is p'N=0.0 */
	  double a_p_pN;
	  double b_pN;

	  a_p_pN = 1.0;
	  b_pN = 0.0;

	   /* this is the matrix A that stores the a_p, a_w and a_e and the vector
	    * b that stores the b' values to solve the system of equations for
	    * the pressure corrections. The system of equations is given by the pressure
	    * corrections equations. We solve the system A*p'=b and the solutions of p'
	    */
	   MatrixXd A_pressure = MatrixXd::Zero(N,N);
	   VectorXd b_pressure = VectorXd::Zero(N);

	   /* Setting the boundary nodes coefficients */
	   A_pressure(0,0) = a_p_p1;
	   b_pressure(0) = b_p1;

	   A_pressure(N-1,N-1) = a_p_pN;
	   b_pressure(N-1) = b_pN;

	   /* adding the values of the residuals of the pressure equation
	    * of boundaries */
	   pressure_residual_sum[i_iter] = abs(A_pressure(0,0)*pressure_prime[0] + A_pressure(0,1)*pressure_prime[1] - b_pressure[0]);
	   pressure_residual_sum[i_iter] = pressure_residual_sum[i_iter] + abs(A_pressure(N-1,N-1)*pressure_prime[N-1] + A_pressure(N-1,N-2)*pressure_prime[N-2] - b_pressure[N-1]);

	   /* For pressure nodes in between  */
	      /* i=1..N-2 */
		  double a_p_p;
		  double a_e_p;
		  double a_w_p;
		  double b_p;

	      for(int i=1; i<N-1; i++){

	   	     a_e_p = density * Area_velocity_node[i] * dx[i];
	   	     a_w_p = density * Area_velocity_node[i-1] * dx[i-1];
	   	     a_p_p = a_e_p + a_w_p;
	   	     b_p = density * Area_velocity_node[i-1] * u_star[i-1] - density * Area_velocity_node[i] * u_star[i];

	   	     A_pressure(i,i) = a_p_p;
	   	     A_pressure(i,i-1) = -a_w_p;
	   	     A_pressure(i,i+1) = -a_e_p;
	   	     b_pressure(i) = b_p;

		     /* adding the values of the residuals of the pressure equations
		      * for the nodes in between */
		     pressure_residual_sum[i_iter] = pressure_residual_sum[i_iter] + abs(A_pressure(i,i)*pressure_prime[i] + A_pressure(i,i+1)*pressure_prime[i+1] + A_pressure(i,i-1)*pressure_prime[i-1] - b_pressure[i]);

	      }


	      /* These are the pressure corrections obtained by solving the system of
	       * equations.
	       */
	      pressure_prime = A_pressure.colPivHouseholderQr().solve(b_pressure);

	      /* We know that the first and last values are zero, so we set them
	       * to avoid error spreading */
	      pressure_prime[0] = 0.0;
	      pressure_prime[N-1] = 0.0;

//   	      if (i_iter == 0){
//
//   	      cout << "Pressure matrix A:" << endl;
//   	      cout << A_pressure << endl;
//   	      cout << endl;
//   	      cout << "Pressure source vector b:" << endl;
//   	      cout << b_pressure << endl;
//
//          cout << endl;
//		  cout << "# iterations" << "\t" << "x momentum residual" << "\t" << "pressure residual" << endl;
//   	      }

}
