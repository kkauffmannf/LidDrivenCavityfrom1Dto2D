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

	  for(int i=ngcx; i<(Nodesx+ngcx); i++){
		  for(int j=ngcy;j<(Nodesy+ngcy);j++){
			  pressure[i][j] = pressure[i][j] + urfp * pressure_prime(i,j);
			  u_velocity[i][j] = u_star(i,j) + d_u[i][j] * ( pressure_prime((i-1),j) - pressure_prime(i,j) );
			  v_velocity[i][j] = v_star(i,j) + d_v[i][j] * ( pressure_prime(i,(j-1)) - pressure_prime(i,j) );
		   }
	  }

}
