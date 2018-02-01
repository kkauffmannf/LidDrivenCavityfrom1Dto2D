/*
 * initialization.cpp
 *
 * Initializes the variables
 *

 *
 *
 *  Created on: Jul 05, 2017
 *      Author: Karla Kauffmann
 *      E-mail: karla.kauffmann@cchen.cl
 */

#include <iostream>
#include <string>
#include <vector>

#include "globals.h"

using namespace std;

void initialization()
{
	/* initialize values of the x-momentum residual */
		   for (int i=0;i<(MAX_ITER);i++) {
			       x_momentum_residual_sum[i] = 0.0;
		         }

		   /* initialize values of the x-momentum residual */
		   for (int i=0;i<(MAX_ITER);i++) {
			       y_momentum_residual_sum[i] = 0.0;
		         }

		   /* initialize values of the pressure residual */
		   for (int i=0;i<(MAX_ITER);i++) {
			       pressure_residual_sum[i] = 0.0;
		         }

		   /* initialize values of the position of the pressure and velocity nodes.
		    * It is a uniform grid. */
		   double delta_x = Lx/(Nx-0.5); /* spacing of the grid in x (uniform grid) */
		   double delta_y = Ly/(Ny-0.5); /* spacing of the grid in y (uniform grid) */

		   for (int i=0;i<Nx;i++) {
		   	   position_pressure_node_x[i] = delta_x*(i+0.5);
		   }

		   for (int j=0;j<Ny;j++) {
		   	   position_pressure_node_y[j] = delta_y*(j+0.5);
		   }

		   for (int i=0;i<Nx;i++) {
		       position_u_velocity_node_x[i] = delta_x*i;
		      }

		   for (int j=0;j<Ny;j++) {
		       position_u_velocity_node_y[j] = delta_y*(j+0.5);
		      }

		   for (int i=0;i<Nx;i++) {
		       position_v_velocity_node_x[i] = delta_x*(i+0.5);
		      }

		   for (int j=0;j<Ny;j++) {
		       position_v_velocity_node_y[j] = delta_y*j;
		      }

		   /* initialize values of Areas pressure and velocity nodes */
		   for (int i=0;i<Nx;i++) {
			   for (int j=0;j<Ny;j++) {
				   Area_velocity_node_u[i][j] = delta_x*delta_y;
			   }
		   }

		   for (int i=0;i<Nx;i++) {
			   for (int j=0;j<Ny;j++) {
				   Area_velocity_node_v[i][j] = delta_x*delta_y;
			   }
		   }

		   for (int i=0;i<Nx;i++) {
			   for (int j=0;j<Ny;j++) {
				   Area_pressure_node[i][j] = delta_x*delta_y;
			   }
		   }

		   /* initialize values of velocities */
		   for (int i=0;i<Nx;i++) {
			   for (int j=0;j<Ny;j++) {
		          u_velocity[i][j] = 0.0;
		          if (j==(Ny-1)){
		              u_velocity[i][j] = lid_velocity;
		          }
			   }
		   }

		   for (int i=0;i<Nx;i++) {
			   for (int j=0;j<Ny;j++) {
		   	      v_velocity[i][j] = 0.0;
			   }
		   }

		   /* initialize values of old velocities, in this case, the same as the actual
		    * velocities because we haven't performed any iteration yet */
		   for (int i=0;i<Nx;i++) {
			   for (int j=0;j<Ny;j++) {
		          u_velocity_old[i][j] = 0.0;
		          if (j==(Ny-1)){
			          u_velocity_old[i][j] = lid_velocity;
		          }
			   }
		   }

		   for (int i=0;i<Nx;i++) {
			   for (int j=0;j<Ny;j++) {
		   	      v_velocity_old[i][j] = 0.0;
			   }
		   }

		   /* initialize values of the pressure */
		   for (int i=0;i<Nx;i++) {
			   for (int j=0;j<Ny;j++) {
		          pressure[i][j] = p_init;
			   }
		   }

		   /* initialize values of old pressures, in this case, the same as the actual
		    * pressures because we haven't performed any iteration yet */
		   for (int i=0;i<Nx;i++) {
			   for (int j=0;j<Ny;j++) {
		          pressure_old[i][j] = p_init;
			   }
		   }

		   /* initialize values of the parameter d for the pressure correction equation */
		   for (int i=0;i<Nx;i++) {
			   for (int j=0;j<Ny;j++) {
				   d_u[i][j] = 0.0;
			   }
		   }

		   for (int i=0;i<Nx;i++) {
		  		   for (int j=0;j<Ny;j++) {
		  			   d_v[i][j] = 0.0;
		  		   }
		  	   }
}
