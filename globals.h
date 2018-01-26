/*
 * globals.h
 *
 * Declares the global variables. They are declared here
 * to be used across the different source files
 *
 *  Created on: Jul 5, 2017
 *      Author: Karla Kauffmann
 *      E-mail: karla.kauffmann@cchen.cl
 */

#include <vector>
#include <Eigen/SVD>
#include <Eigen/Dense>

using namespace Eigen;

/* Variables read in input.txt */
extern int N;                       /* number of pressure nodes for the CV computation */
extern double L;                    /* length of the rod */
extern double p_in;                 /* stagnation pressure at inlet of the nozzle */
extern double p_out;                /* static pressure at the exit of the nozzle */
extern double density;              /* density of the fluid */
extern double Area_in;              /* cross-sectional area at the inlet */
extern double Area_out;             /* cross-sectional area at the exit */
extern double m_dot;                /* guessed mass flow rate to generate the an initial velocity field */
extern double urfu;                 /* under-relaxation factor for velocity u. Cannot be 0.0 since it goes in the denominator
                                    of a_p in the momentum equations */
extern double urfp;                 /* under-relaxation factor for pressure */
extern double residual_threshold;   /* the value of the residuals to end the iterations */

/* Physical variables */
extern std::vector<double> position_pressure_node; /* position of the pressure nodes */
extern std::vector<double> position_velocity_node; /* position of the velocity nodes*/
extern std::vector<double> Area_pressure_node;     /* cross section at pressure nodes */
extern std::vector<double> Area_velocity_node;     /* cross section at velocity nodes */
extern std::vector<double> velocity;
extern std::vector<double> velocity_old;           /* the value of the velocity in the previous iteration */
extern std::vector<double> pressure;
extern std::vector<double> pressure_old;           /* the value of the pressure in the previous iteration */
extern std::vector<double> dx;                     /* parameter d for the pressure correction equation */

/* Variables used for the iterations */
extern int i_iter;                                   /* number of iterations */
extern int MAX_ITER;                                 /* set the maximum number of iterations to store in the residual vector */
extern std::vector<double> x_momentum_residual_sum;  /* sum of the residuals of the x-momentum equation per iteration*/
extern std::vector<double> pressure_residual_sum;    /* sum of the residuals of the pressure equation per iteration*/

/* Modules external to main */
void readinput();
void initialization();
void plotcolormap();
void plotresiduals();
void momentum_equation_solve(VectorXd &u_star, int i_iter);
void pressure_correction_equation_solve(VectorXd u_star, VectorXd &pressure_prime, int i_iter);
void correct_pressure_and_velocities(VectorXd u_star,VectorXd pressure_prime);
void underrelaxation();
