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
		   vector<double> pressure_residual_sum_prev;

		   pressure_residual_sum_prev.resize(MAX_ITER);

		   /* initialize values of the previous pressure residual */
		   for (int i=0;i<(MAX_ITER);i++) {
			   pressure_residual_sum_prev[i] = pressure_residual_sum[i];
		   }

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
				   pressure_residual_sum[i_iter] = pressure_residual_sum[i_iter] + abs(A_pressure(j,j)*pressure_prime(i,j) + A_pressure(j,j+1)*pressure_prime(i,j+1) + A_pressure(j,j-1)*pressure_prime(i,j-1) - b_pressure[j]);
			   }
			   /* keeping the highest residual, line by line (ith row) */
			   if (pressure_residual_sum[i_iter] < pressure_residual_sum_prev[i_iter]) {
				   pressure_residual_sum[i_iter]=pressure_residual_sum_prev[i_iter];
			   }
//			   				      	      cout << "Pressure matrix A:" << endl;
//			   				      	      cout << A_pressure << endl;
//			   				      	      cout << endl;
//			   				      	      cout << "Pressure source vector b:" << endl;
//			   				      	      cout << b_pressure << endl;
			   /* We also solve pressure_prime */
			   pressure_prime.row(i) = A_pressure.colPivHouseholderQr().solve(b_pressure);
		   }

		   /* BC */
		   /* west guard cells */
		   for (int i=0;i<ngcx;i++){
			   for (int j=0;j<Ny;j++){
				   pressure[i][j] = pressure[ngcx][j];
			   }

		   }
		   /* east guard cells */
		   for (int i=(Nodesx + ngcx - 1);i<Nx;i++){
			   for (int j=0;j<Ny;j++){
				   pressure[i][j] = pressure[(Nodesx + ngcx - 2)][j];
			   }
		   }

		   /* south guard cells */
		   for (int i=0;i<Nx;i++){
			   for (int j=0;j<ngcy;j++){
				   pressure[i][j] = pressure[i][ngcy];
			   }
		   }

		   /* north guard cells */
		   for (int i=0;i<Nx;i++){
			   for (int j=(Nodesy + ngcy - 1);j<Ny;j++){
				   pressure[i][j] = pressure[i][(Nodesy + ngcy - 2)];
			   }
		   }
	}




//		   for(int i=0; i<Ny; i++){
//			   for (int j=0;j<Nx;j++) {
//				   if(i==0) {
//					   if (j==0) {
//						   a_pressure[i+1][j]=0.0;
//						   a_pressure[i][j+1]=0.0;
//						   a_pressure[i][j] = 1.0;
//						   Source_pressure[i][j]= 0.0;
//					   }else if (j==(Nx-1)) {
//						   a_pressure[i+1][j]=0.0;
//						   a_pressure[i][j-1]=0.0;
//						   a_pressure[i][j] = 1.0;
//						   Source_pressure[i][j]= 0.0;
//					   } else {
//						   a_pressure[i+1][j]=0.0;
//						   a_pressure[i][j+1]=0.0;
//						   a_pressure[i][j-1]=0.0;
//						   a_pressure[i][j] = 1.0;
//						   Source_pressure[i][j]= 0.0;
//					   }
//				   } else if (i==(Ny-1)) {
//					   if (j==0) {
//						   a_pressure[i-1][j]=0.0;
//						   a_pressure[i][j+1]=0.0;
//						   a_pressure[i][j] = 1.0;
//						   Source_pressure[i][j]= 0.0;
//					   } else if (j==(Nx-1)) {
//						   a_pressure[i-1][j]=0.0;
//						   a_pressure[i][j-1]=0.0;
//						   a_pressure[i][j] = 1.0;
//						   Source_pressure[i][j]= 0.0;
//					   } else {
//						   a_pressure[i-1][j]=0.0;
//						   a_pressure[i][j+1]=0.0;
//						   a_pressure[i][j-1]=0.0;
//						   a_pressure[i][j] = 1.0;
//						   Source_pressure[i][j]= 0.0;
//					   }
//				   } else if (i==(Ny-2)) {
//					       if (j==0) {
//						   	   a_pressure[i+1][j]=0.0;
//						   	   a_pressure[i-1][j]=0.0;
//						   	   a_pressure[i][j+1]=0.0;
//						   	   a_pressure[i][j] = 1.0;
//						   	   Source_pressure[i][j]= 0.0;
//					       } else if (j==(Nx-1)) {
//						   	   a_pressure[i+1][j]=0.0;
//						   	   a_pressure[i-1][j]=0.0;
//						   	   a_pressure[i][j-1]=0.0;
//						   	   a_pressure[i][j] = 1.0;
//						   	   Source_pressure[i][j]= 0.0;
//					       } else if (j==(Nx-2)) {
//						   	   a_pressure[i+1][j]=0.0;
//						   	   a_pressure[i-1][j]=0.0;
//						   	   a_pressure[i][j+1]=0.0;
//						   	   a_pressure[i][j-1]=0.0;
//						   	   a_pressure[i][j] = 1.0;
//						   	   Source_pressure[i][j]= 0.0;
//	//					   	   a_pressure[i+1][j]=density*d_u[i][j]*Area_velocity_node_u[i][j];
//	//					   	   a_pressure[i-1][j]=density*d_u[i][j]*Area_velocity_node_u[i][j];
//	//					   	   a_pressure[i][j+1]=density*d_v[i][j]*Area_velocity_node_v[i][j];
//	//					   	   a_pressure[i][j-1]=density*d_v[i][j]*Area_velocity_node_v[i][j];
//	//					   	   a_pressure[i][j] = a_pressure[i+1][j] + a_pressure[i-1][j] + a_pressure[i][j+1] + a_pressure[i][j-1];
//	//					   	   Source_pressure[i][j]= density*u_star(i,j)*Area_velocity_node_u[i][j] - density*u_star(i,j)*Area_velocity_node_u[i][j] + density*v_star(i,j)*Area_velocity_node_v[i][j] - density*v_star(i,j)*Area_velocity_node_v[i][j];
//				   		   } else {
//						   	   a_pressure[i+1][j]=density*d_u[i][j]*Area_velocity_node_u[i][j];
//						   	   a_pressure[i-1][j]=density*d_u[i][j]*Area_velocity_node_u[i][j];
//						   	   a_pressure[i][j+1]=density*d_v[i][j+1]*Area_velocity_node_v[i][j+1];
//						   	   a_pressure[i][j-1]=density*d_v[i][j]*Area_velocity_node_v[i][j];
//						   	   a_pressure[i][j] = a_pressure[i+1][j] + a_pressure[i-1][j] + a_pressure[i][j+1] + a_pressure[i][j-1];
//						   	   Source_pressure[i][j]= density*u_star(i,j)*Area_velocity_node_u[i][j] - density*u_star(i,j)*Area_velocity_node_u[i][j] + density*v_star(i,j)*Area_velocity_node_v[i][j] - density*v_star(i,j+1)*Area_velocity_node_v[i][j+1];
//				   		   }
//				   } else {
//				       if (j==0) {
//				    	   a_pressure[i+1][j]=0.0;
//				    	   a_pressure[i-1][j]=0.0;
//				    	   a_pressure[i][j+1]=0.0;
//				    	   a_pressure[i][j] = 1.0;
//						   Source_pressure[i][j]= 0.0;
//				       } else if (j==(Nx-1)) {
//				    	   a_pressure[i+1][j]=0.0;
//				    	   a_pressure[i-1][j]=0.0;
//				    	   a_pressure[i][j-1]=0.0;
//				    	   a_pressure[i][j] = 1.0;
//				    	   Source_pressure[i][j]= 1.0;
//				       } else if (j==(Nx-2)) {
//				    	   a_pressure[i+1][j]=0.0;
//				    	   a_pressure[i-1][j]=0.0;
//				    	   a_pressure[i][j+1]=0.0;
//				    	   a_pressure[i][j-1]=0.0;
//				    	   a_pressure[i][j] = 1.0;
//				    	   Source_pressure[i][j]= 0.0;
//	//				   	   a_pressure[i+1][j]=density*d_u[i+1][j]*Area_velocity_node_u[i+1][j];
//	//				   	   a_pressure[i-1][j]=density*d_u[i][j]*Area_velocity_node_u[i][j];
//	//				   	   a_pressure[i][j+1]=density*d_v[i][j]*Area_velocity_node_v[i][j];
//	//				   	   a_pressure[i][j-1]=density*d_v[i][j]*Area_velocity_node_v[i][j];
//	//				   	   a_pressure[i][j] = a_pressure[i+1][j] + a_pressure[i-1][j] + a_pressure[i][j+1] + a_pressure[i][j-1];
//	//				   	   Source_pressure[i][j]= density*u_star(i,j)*Area_velocity_node_u[i][j] - density*u_star(i+1,j)*Area_velocity_node_u[i+1][j] + density*v_star(i,j)*Area_velocity_node_v[i][j] - density*v_star(i,j)*Area_velocity_node_v[i][j];
//				       } else {
//					   	   a_pressure[i+1][j]=density*d_u[i+1][j]*Area_velocity_node_u[i+1][j];
//					   	   a_pressure[i-1][j]=density*d_u[i][j]*Area_velocity_node_u[i][j];
//					   	   a_pressure[i][j+1]=density*d_v[i][j+1]*Area_velocity_node_v[i][j+1];
//					   	   a_pressure[i][j-1]=density*d_v[i][j]*Area_velocity_node_v[i][j];
//					   	   a_pressure[i][j] = a_pressure[i+1][j] + a_pressure[i-1][j] + a_pressure[i][j+1] + a_pressure[i][j-1];
//					   	   Source_pressure[i][j]= density*u_star(i,j)*Area_velocity_node_u[i][j] - density*u_star(i+1,j)*Area_velocity_node_u[i+1][j] + density*v_star(i,j)*Area_velocity_node_v[i][j] - density*v_star(i,j+1)*Area_velocity_node_v[i][j+1];
//				       }
//				   }
//				   for (int j=0; j<Nx; j++) {
//				   			   A_pressure(j,j) = a_pressure[i][j];
//				   			   if (j>0) {
//				   			       A_pressure(j,j-1) = -a_pressure[i][j-1];
//				   			   }
//				   			   if (j<(Nx-1)) {
//				   			   	   A_pressure(j,j+1) = -a_pressure[i][j+1];
//				   			   }
//				   			   if (i==0) {
//				   				   b_pressure(j) = a_pressure[i+1][j]*pressure_prime(i+1,j) + Source_pressure[i][j];
//				   			   } else if(i==(Ny-1)) {
//				   				   b_pressure(j) = a_pressure[i-1][j]*pressure_prime(i-1,j) + Source_pressure[i][j];
//				   			   } else {
//				   			       b_pressure(j) = a_pressure[i-1][j]*pressure_prime(i-1,j) + a_pressure[i+1][j]*pressure_prime(i+1,j) + Source_pressure[i][j];
//				   			   }
//				   		     /* adding the values of the residuals of the pressure equations */
//				   		     if (j==0){
//				   		         pressure_residual_sum[i_iter] = pressure_residual_sum[i_iter] + abs(A_pressure(j,j)*pressure_prime(i,j) + A_pressure(j,j+1)*pressure_prime(i,j+1) - b_pressure[j]);
//				   		     }else if (j==(Nx-1)){
//				   		    	 pressure_residual_sum[i_iter] = pressure_residual_sum[i_iter] + abs(A_pressure(j,j)*pressure_prime(i,j) + A_pressure(j,j-1)*pressure_prime(i,j-1) - b_pressure[j]);
//				   		     }else {
//				   		    	 pressure_residual_sum[i_iter] = pressure_residual_sum[i_iter] + abs(A_pressure(j,j)*pressure_prime(i,j) + A_pressure(j,j+1)*pressure_prime(i,j+1) + A_pressure(j,j-1)*pressure_prime(i,j-1) - b_pressure[j]);
//				   		     }
//
//				   		     /* keeping the highest residual, line by line (ith row) */
//				   		     if (pressure_residual_sum[i_iter] < pressure_residual_sum_prev[i_iter]) {
//				   		    	 pressure_residual_sum[i_iter]=pressure_residual_sum_prev[i_iter];
//				   		     }
//				   		   }
//
//				     	      /* These are the pressure corrections obtained by solving the system of
//				     	       * equations.
//				      	       */
//
//
//	//			   		   for (int i=0;i<Ny;i++) {
//	//			   		       for (int j=0;j<Nx;j++) {
//	//			   		    	   cout << A_pressure(i,j) << "\t";
//	//			   		       }
//	//			   		       cout << endl;
//	//			   		   }
//	//			   		   cout << endl;
//
//				       	   pressure_prime.row(i) = A_pressure.colPivHouseholderQr().solve(b_pressure);
//
//				   		   }
//			   }
//	   	       cout << "PressurePrime:" << endl;
//			   for (int i=0;i<Ny;i++) {
//			       for (int j=0;j<Nx;j++) {
//			    	   cout << pressure_prime(i,j) << "\t";
//			       }
//			       cout << endl;
//
//	   		   cout << "A matrix" << "\t" << "i:" << i << endl;
//	       	   cout << A_pressure << endl;
//	       	   cout << "b vector" << endl;
//	       	   cout << b_pressure << endl;
//	   		   cout << endl;
//		   }
//
//
//
//		   ////////////////////////////////////////////////////////////
//
//		   for(int i=0; i<Ny; i++){
//
//	//    	      cout << "Pressure matrix A:" << endl;
//	//    	      cout << A_pressure << endl;
//	//    	      cout << endl;
//	//    	      cout << "Pressure source vector b:" << endl;
//	//    	      cout << b_pressure << endl;
//		   }
//		   /* setting the boundary values. P' equals zero because
//		    * no corrections are used for the pressure at the boundaries. */
//		   for (int i=0; i<Ny;i++){
//			   /* left boundary */
//			   pressure_prime(i,0)=0.0;
//			   /* right boundary */
//			   pressure_prime(i,Nx-1)=0.0;
//		   }
//		   for (int j=0; j<Nx;j++){
//			   /* bottom boundary */
//			   pressure_prime(0,j)=0.0;
//			   /* top boundary */
//			   pressure_prime(Ny-1,j)=0.0;
//		   }
//
//
//
//	///////////////////////////////////////////////////////////
//	//
//	//	      /* These are the pressure corrections obtained by solving the system of
//	//	       * equations.
//	//	       */
//	//	      pressure_prime = A_pressure.colPivHouseholderQr().solve(b_pressure);
//	//
//	//	      /* We know that the first and last values are zero, so we set them
//	//	       * to avoid error spreading */
//	//	      pressure_prime[0] = 0.0;
//	//	      pressure_prime[N-1] = 0.0;
//	//
//	//   	      if (i_iter == 0){
//	//
//	//   	      cout << "Pressure matrix A:" << endl;
//	//   	      cout << A_pressure << endl;
//	//   	      cout << endl;
//	//   	      cout << "Pressure source vector b:" << endl;
//	//   	      cout << b_pressure << endl;
//	//   	      cout << "Pressure corrections:" << endl;
//	//   	      cout << pressure_prime << endl;
//	//
//	//          cout << endl;
//	//		  cout << "# iterations" << "\t" << "x momentum residual" << "\t" << "pressure residual" << endl;
//	//   	      }
























//	   /* Solve the discretised pressure correction equation using the guessed values of
//	    * velocity and the values of pressure. The equation has the form:
//	    *
//	    * a_P * p'_P = a_W * p'_W + a_E * p'_E + b'
//	    *
//	    * where p' is the pressure corrections p'=(p-p*) in each node.
//	    *
//	    * */
//
//
//	   /* First we determine the values for the boundary nodes.
//	    * For the first and Nx-1 pressure nodes the value of the
//	    * correction is set to 0.0 (see reference Versteeg) */
//
//	   /* For pressure node A
//	    *
//	    * the equation is p'1=0.0 */
//	   double a_p_p1;
//	   double b_p1;
//
//	   a_p_p1 = 1.0;
//	   b_p1 = 0.0;
//
//	  /* For pressure node Nx
//	   *
//	   * the equation is p'Nx=0.0 */
//	  double a_p_pN;
//	  double b_pN;
//
//	  a_p_pN = 1.0;
//	  b_pN = 0.0;
//
//	   /* this is the matrix A that stores the a_p, a_w and a_e and the vector
//	    * b that stores the b' values to solve the system of equations for
//	    * the pressure corrections. The system of equations is given by the pressure
//	    * corrections equations. We solve the system A*p'=b and the solutions of p'
//	    */
//	   MatrixXd A_pressure = MatrixXd::Zero(Nx,Nx);
//	   VectorXd b_pressure = VectorXd::Zero(Nx);
//
//	   /* Setting the boundary nodes coefficients */
//	   A_pressure(0,0) = a_p_p1;
//	   b_pressure(0) = b_p1;
//
//	   A_pressure(Nx-1,Nx-1) = a_p_pN;
//	   b_pressure(Nx-1) = b_pN;
//
//	   /* adding the values of the residuals of the pressure equation
//	    * of boundaries */
//	   pressure_residual_sum[i_iter] = abs(A_pressure(0,0)*pressure_prime[0] + A_pressure(0,1)*pressure_prime[1] - b_pressure[0]);
//	   pressure_residual_sum[i_iter] = pressure_residual_sum[i_iter] + abs(A_pressure(Nx-1,Nx-1)*pressure_prime[Nx-1] + A_pressure(Nx-1,Nx-2)*pressure_prime[Nx-2] - b_pressure[Nx-1]);
//
//	   /* For pressure nodes in between  */
//	      /* i=1..Nx-2 */
//		  double a_p_p;
//		  double a_e_p;
//		  double a_w_p;
//		  double b_p;
//
//	      for(int i=1; i<Nx-1; i++){
//
//	   	     a_e_p = density * Area_velocity_node[i] * dx[i];
//	   	     a_w_p = density * Area_velocity_node[i-1] * dx[i-1];
//	   	     a_p_p = a_e_p + a_w_p;
//	   	     b_p = density * Area_velocity_node[i-1] * u_star[i-1] - density * Area_velocity_node[i] * u_star[i];
//
//	   	     A_pressure(i,i) = a_p_p;
//	   	     A_pressure(i,i-1) = -a_w_p;
//	   	     A_pressure(i,i+1) = -a_e_p;
//	   	     b_pressure(i) = b_p;
//
//		     /* adding the values of the residuals of the pressure equations
//		      * for the nodes in between */
//		     pressure_residual_sum[i_iter] = pressure_residual_sum[i_iter] + abs(A_pressure(i,i)*pressure_prime[i] + A_pressure(i,i+1)*pressure_prime[i+1] + A_pressure(i,i-1)*pressure_prime[i-1] - b_pressure[i]);
//
//	      }
//
//
//	      /* These are the pressure corrections obtained by solving the system of
//	       * equations.
//	       */
//	      pressure_prime = A_pressure.colPivHouseholderQr().solve(b_pressure);
//
//	      /* We know that the first and last values are zero, so we set them
//	       * to avoid error spreading */
//	      pressure_prime[0] = 0.0;
//	      pressure_prime[Nx-1] = 0.0;
//
////   	      if (i_iter == 0){
////
////   	      cout << "Pressure matrix A:" << endl;
////   	      cout << A_pressure << endl;
////   	      cout << endl;
////   	      cout << "Pressure source vector b:" << endl;
////   	      cout << b_pressure << endl;
////
////          cout << endl;
////		  cout << "# iterations" << "\t" << "x momentum residual" << "\t" << "pressure residual" << endl;
////   	      }


