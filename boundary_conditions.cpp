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
#include <Eigen/SVD>
#include <Eigen/Dense>
#include "globals.h"

using namespace std;

void boundary_conditions()
{

	/* west guard cells */
	for (int i=0;i<ngcx;i++){
		for (int j=0;j<Ny;j++){
			u_velocity[i][j] = 0.0;
			v_velocity[i][j] = 0.0;
			pressure[i][j] = pressure[ngcx][j];
		}

	}

	/* east guard cells */
	for (int i=(Nodesx + ngcx);i<Nx;i++){
		for (int j=0;j<Ny;j++){
			u_velocity[i][j] = 0.0;
			v_velocity[i][j] = 0.0;
			pressure[i][j] = pressure[(Nodesx + ngcx - 1)][j];
		}

	}

	/* south guard cells */
	for (int i=0;i<Nx;i++){
		for (int j=0;j<ngcy;j++){
			u_velocity[i][j] = 0.0;
			v_velocity[i][j] = 0.0;
			pressure[i][j] = pressure[i][ngcy];
		}

	}

	/* north guard cells */
	for (int i=0;i<Nx;i++){
		for (int j=(Nodesy + ngcy);j<Ny;j++){
			u_velocity[i][j] = 0.0;
			v_velocity[i][j] = 0.0;
			pressure[i][j] = pressure[i][(Nodesy + ngcy - 1)];
		}

	}

	/* BC. We set the guard cells adjacent to nodes that are not exactly on the boundary, equal to the  */
	/* negative velocity, so when they sum, it is equal to zero */
	for (int j=ngcy;j<(Nodesy + ngcy);j++){
		u_velocity[(Nodesx + ngcx)][j] = - u_velocity[(Nodesx + ngcx - 1)][j];
		v_velocity[(ngcx - 1)][j] = - v_velocity[ngcx][j];
		u_velocity[ngcx][j] = 0.0;
		v_velocity[(Nodesx + ngcx - 1)][j] = 0.0;
	}

	for (int i=ngcx;i<(Nodesx + ngcx);i++) {
		u_velocity[i][(ngcy - 1)] = - u_velocity[i][ngcy];
		v_velocity[i][(Nodesy + ngcy)] = - v_velocity[i][(Nodesy + ngcy - 1)];
		u_velocity[i][(Nodesy + ngcy - 1)] = lid_velocity;
		v_velocity[i][ngcy] = 0.0;
	}
//	u_velocity[ngcx][ngcy] = 0.0;
	v_velocity[ngcx][ngcy] = 0.0;
	u_velocity[(Nodesx + ngcx - 1)][(Nodesy + ngcy - 1)] = lid_velocity;
//	v_velocity[(Nodesx + ngcx - 1)][(Nodesy + ngcy - 1)] = 0.0;
}
