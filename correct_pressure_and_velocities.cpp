/*
 * correct_pressure_and_velocities.cpp
 *
 * Corrects pressure and velocities using the guessed velocities u* and the
 * pressure corrections p'
 *
 *  Created on: Jul 18, 2017
 *     Author: Karla Kauffmann
 *     E-mail: karla.kauffmann@cchen.cl
 */

#include <iostream>
#include <string>
#include <vector>
#include <Eigen/SVD>
#include <Eigen/Dense>
#include "globals.h"

using namespace std;

void correct_pressure_and_velocities(MatrixXd u_star, MatrixXd v_star, MatrixXd pressure_prime)
{

	  /* We set the new values for the pressures due to the corrections */

	  for(int i=0;i<Nx;i++){
		  for(int j=0;j<Ny;j++){
			  pressure[i][j] = pressure[i][j] + pressure_prime(i,j);
			  /* The corrected velocities are */
			  if (i==0){
				  if (j==0){
					  u_velocity[i][j] = 0.0;
					  v_velocity[i][j] = v_star(i,j);
				  }
				  else {
					  u_velocity[i][j] = u_star(i,j);
					  v_velocity[i][j] = v_star(i,j) + d_v[i][j] * ( pressure_prime(i,(j-1)) - pressure_prime(i,j) );
				  }
			  }
			  else {
				  if (j==0){
					  u_velocity[i][j] = u_star(i,j) + d_u[i][j] * ( pressure_prime((i-1),j) - pressure_prime(i,j) );
					  v_velocity[i][j] = v_star(i,j);
				  }
				  else {
					  u_velocity[i][j] = u_star(i,j) + d_u[i][j] * ( pressure_prime((i-1),j) - pressure_prime(i,j) );
					  v_velocity[i][j] = v_star(i,j) + d_v[i][j] * ( pressure_prime(i,(j-1)) - pressure_prime(i,j) );
				  }
			  }
		  }
	  }
}
