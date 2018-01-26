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

	   /* initialize values of the pressure residual */
	   for (int i=0;i<(MAX_ITER);i++) {
		       pressure_residual_sum[i] = 0.0;
	         }

	   /* initialize values of the position of the pressure and velocity nodes.
	    * It is a uniform grid. */
	   double delta_x = L/(N-1); /* spacing of the grid (uniform grid) */

	   for (int i=0;i<N;i++) {
	   	   position_pressure_node[i] = delta_x*i;
	      }

	   for (int i=0;i<(N-1);i++) {
	       position_velocity_node[i] = delta_x*(i+0.5);
	      }

	   /* initialize values of Areas for N pressure and (N-1) velocity nodes.
	     * It decreases linearly. */
	   for (int i=0;i<N;i++) {
		   Area_pressure_node[i] = Area_in - ((Area_in-Area_out)/(N-1)*i);
	   }

	   for (int i=0;i<(N-1);i++) {
	   	   Area_velocity_node[i] = Area_in - ((Area_in-Area_out)/(N-1)*(i+0.5));
	      }

	   /* initialize values of velocities using the mass flow */
	   for (int i=0;i<(N-1);i++) {
	      	   velocity[i] = m_dot/(density*Area_velocity_node[i]);
	         }

	   /* initialize values of old velocities, in this case, the same as the actual
	    * velocities because we haven't performed any iteration yet */
	   for (int i=0;i<(N-1);i++) {
	      	   velocity_old[i] = m_dot/(density*Area_velocity_node[i]);
	         }

	   /* initialize values of the pressure
	    * we assume a linear variation between the inlet and outlet */
	   for (int i=0;i<N;i++) {
	      	   pressure[i] = p_in - (p_in - p_out)/(N-1)*i;
	         }

	   /* initialize values of old pressures, in this case, the same as the actual
	    * pressures because we haven't performed any iteration yet */
	   for (int i=0;i<N;i++) {
	           pressure_old[i] = p_in - (p_in - p_out)/(N-1)*i;
	         }

	   /* initialize values of the parameter d for the pressure correction equation */
	   for (int i=0;i<(N-1);i++) {
			   dx[i] = 0.0;
		     }
}
