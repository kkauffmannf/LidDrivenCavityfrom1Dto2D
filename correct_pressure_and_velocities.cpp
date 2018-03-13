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
//#include <Eigen/SVD>
//#include <Eigen/Dense>
#include "globals.h"

using namespace std;

void correct_pressure_and_velocities(vector<vector<double>> u_star, vector<vector<double>> v_star, vector<vector<double>> pressure_prime)
{

	  /* We set the new values for the pressures due to the corrections */

	  for(int i=ngc; i<(Npx-ngc);i++){
		  for(int j=ngc; j<(Npy-ngc);j++){
			  pressure[i][j] = pressure_old[i][j] + urfp * pressure_prime[i][j];
		  }
	  }

	  for(int i=ngc; i<(Nux-ngc);i++){
		  for(int j=ngc; j<(Nuy-ngc);j++){
			  u_velocity[i][j] = u_star[i][j] + d_u[i][j] * ( pressure_prime[(i-1)][j] - pressure_prime[i][j] );
		  }
	  }
	  for(int i=ngc; i<(Nvx-ngc);i++){
		  for(int j=ngc; j<(Nvy-ngc);j++){
			  v_velocity[i][j] = v_star[i][j] + d_v[i][j] * ( pressure_prime[i][(j-1)] - pressure_prime[i][j] );
		  }
	  }


}
