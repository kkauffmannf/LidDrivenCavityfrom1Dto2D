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

void underrelaxation()
{
	      /* The underrelaxation imposes */
	      for(int i=0;i<(Nx-1);i++){
	          	velocity[i] = (1.0-urfu)*velocity_old[i] + urfu*velocity[i];
	      }

	      for(int i=0;i<Nx;i++){
	          	pressure[i] = (1.0-urfp)*pressure_old[i] + urfp*pressure[i];
	      }
}
