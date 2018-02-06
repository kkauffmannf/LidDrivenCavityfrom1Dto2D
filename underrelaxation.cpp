/*
 * underrelaxation.cpp
 *
 * This module applies the underrelaxation factors for velocity (urfu) and pressure (urfp)
 *
 *  Created on: Jul 18, 2017
 *		Author: Karla Kauffmann
 *      E-mail: karla.kauffmann@cchen.cl
 *
 */

#include <iostream>
#include <string>
#include <vector>
#include <Eigen/SVD>
#include <Eigen/Dense>
#include "globals.h"

using namespace std;

void underrelaxation(MatrixXd pressure_prime)
{
	      /* The underrelaxation imposes */
	      for(int i=0;i<Nx;i++){
	    	  for(int j=0;j<Ny;j++){
	    		  u_velocity[i][j] = (1.0-urfu)*u_velocity_old[i][j] + urfu * u_velocity[i][j];
	    		  if (j!=(Ny-1)){
	    			  /* west boundary */
	                  u_velocity[0][j] = 0.0;
	    		  }

	    		  v_velocity[i][j] = (1.0-urfv)*v_velocity_old[i][j] + urfv * v_velocity[i][j];
	    		  /* south boundary */
	    		  v_velocity[i][0] = 0.0;
	    		  /* east boundary */
	    		  v_velocity[(Nx-1)][j] = 0.0;

	    		  pressure[i][j] = (1.0-urfp)*pressure_old[i][j] + urfp*pressure_prime(i,j);
	    	  }
	    	  /* north boundary */
	    	  u_velocity[i][(Ny-1)] = lid_velocity;
	      }
		  /* For the bottom west corner we set it constant */
		  pressure[0][0] = p_init;
}
