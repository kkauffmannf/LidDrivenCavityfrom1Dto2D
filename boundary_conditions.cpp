/*
 * boundary_conditions.cpp
 *
 * This module applies the boundary conditions to the guard cells.
 *
 *  Created on: Feb 7, 2018
 *		Author: Karla Kauffmann
 *      E-mail: karla.kauffmann@cchen.cl
 *
 */

#include <iostream>
#include <string>
#include <vector>
//#include <Eigen/SVD>
//#include <Eigen/Dense>
#include "globals.h"

using namespace std;

void boundary_conditions()
{
	/* For nodes that are not exactly on the boundary, we use linear interpolation, */
	/* for example (u[0][0] + u[0][1])/2 = U_BOTTOM, so that u[0][0] = 2*U_BOTTOM - u[0][1]. */
	/* In this case, U_BOTTOM is zero, so u[0][0] = - u[0][1] */

	/* u-velocity */

    /* West u-velocity */
	for(int j=ngc; j < (Nuy - ngc); j++){
		u_velocity[0][j]=0.0;
	}
	u_velocity[0][0] = -u_velocity[0][1];
	u_velocity[0][(Nuy - ngc)] = 2*lid_velocity - u_velocity[0][(Nuy - ngc - 1)];

	/* East u-velocity */
	for(int j=ngc; j < (Nuy - ngc); j++){
		u_velocity[(Nux - ngc - 1)][j]=0.0;
	}
	u_velocity[(Nux - ngc - 1)][0] = -u_velocity[(Nux - ngc - 1)][1];
	u_velocity[(Nux - ngc - 1)][(Nuy - ngc)] = 2*lid_velocity - u_velocity[(Nux - ngc - 1)][(Nuy - ngc -1)];

	/* North u-velocity */
	for(int i=(ngc + 1); i < (Nux - ngc - 1); i++){
		u_velocity[i][(Nuy - ngc)] = 2*lid_velocity - u_velocity[i][(Nuy - ngc -1)];
	}

	/* South u-velocity */
	for(int i=(ngc + 1); i < (Nux - ngc - 1); i++){
		u_velocity[i][0] = - u_velocity[i][1];
	}

	/* West ghost cells */

	for(int i=0; i<ngc;i++){
		for(int j=0; j<Nuy; j++){
			u_velocity[i][j]=0.0;
		}
	}

	/* East ghost cells */

	for(int i=(Nux - ngc); i<Nux;i++){
		for(int j=0; j<Nuy; j++){
			u_velocity[i][j]=0.0;
		}
	}

	/* North ghost cells */
	for(int i=ngc; i<(Nux-ngc);i++){
		for(int j=(Nuy - ngc + 1); j<Nuy; j++){
			u_velocity[i][j]=0.0;
		}
	}

	/* South ghost cells */
	for(int i=ngc; i<(Nux-ngc);i++){
		for(int j=0; j<(ngc - 1); j++){
			u_velocity[i][j]=0.0;
		}
	}

	/* v-velocity */

    /* South v-velocity */
	for(int i=ngc; i < (Nvx - ngc); i++){
		v_velocity[i][0]=0.0;
	}
	v_velocity[0][0] = -v_velocity[1][0];
	v_velocity[(Nvx - ngc)][0] = - v_velocity[(Nvx - ngc - 1)][0];

    /* North v-velocity */
	for(int i=ngc; i < (Nvx - ngc); i++){
		v_velocity[i][(Nvy - ngc - 1)]=0.0;
	}
	v_velocity[0][(Nvy - ngc - 1)] = -v_velocity[1][(Nvy - ngc - 1)];
	v_velocity[(Nvx - ngc)][(Nvy - ngc - 1)] = - v_velocity[(Nvx - ngc - 1)][(Nvy - ngc - 1)];

	/* West v-velocity */
	for(int j=(ngc + 1); j < (Nvy - ngc - 1); j++){
		v_velocity[0][j] = - v_velocity[1][j];
	}

	/* East v-velocity */
	for(int j=(ngc + 1); j < (Nvy- ngc - 1); j++){
		v_velocity[(Nvx - ngc)][j] = - v_velocity[(Nvx - ngc - 1)][j];
	}

	/* South ghost cells */
	for(int j=0; j<ngc;j++){
		for(int i=0; i<Nvx; i++){
			v_velocity[i][j]=0.0;
		}
	}

	/* North ghost cells */
	for(int j=(Nvy - ngc); j<Nvy;j++){
		for(int i=0; i<Nvx; i++){
			v_velocity[i][j]=0.0;
		}
	}

	/* West ghost cells */
	for(int j=ngc; j<(Nvy - ngc);j++){
		for(int i=0; i<(ngc - 1); i++){
			v_velocity[i][j]=0.0;
		}
	}

	/* East ghost cells */
	for(int j=ngc; j<(Nvy - ngc);j++){
		for(int i=(Nvx - ngc + 1); i<Nvx; i++){
			v_velocity[i][j]=0.0;
		}
	}


	/* pressure */

	pressure[ngc][ngc] = p_init;
}
