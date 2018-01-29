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

void correct_pressure_and_velocities(VectorXd u_star, VectorXd pressure_prime)
{

	   /* We set the new values for the pressures due to the corrections */
	   /* For the boundary nodes */
	      pressure[0] = p_in - 0.5*density*u_star[0]*u_star[0]*(Area_velocity_node[0]/Area_pressure_node[0]) * (Area_velocity_node[0]/Area_pressure_node[0]);
	      pressure[Nx-1] = pressure[Nx-1] + pressure_prime[Nx-1];

	      /* For the nodes in between */
	      for(int i=1;i<(Nx-1);i++){
	          pressure[i] = pressure[i] + pressure_prime[i];
	      }

	      /* The corrected velocities are */
	      for(int i=0;i<(Nx-1);i++){
	    	  velocity[i] = u_star[i] + dx[i]*(pressure_prime[i] - pressure_prime[i+1]);
	      }

	      /* The corrected nodal pressure for the first node is */
	      pressure[0] = p_in - 0.5*density*velocity[0]*velocity[0]*(Area_velocity_node[0]/Area_pressure_node[0]) * (Area_velocity_node[0]/Area_pressure_node[0]);

}
