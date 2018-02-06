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

/* Variables in input.txt to be defined at the start by the user. */
extern int Nx; /* number of pressure nodes for the CV computation in the y direction */
extern int Ny; /* number of pressure nodes for the CV computation in the x direction */
extern double Lx; /* length of the cavity in the x direction */
extern double Ly; /* length of the cavity in the y direction */
extern double p_init;
extern double density; /* density of the fluid */
extern double lid_velocity; /* velocity of the lid */
extern double Reynolds_num; /* Reynolds number */
extern double urfu; /* under-relaxation factor for velocity u. Cannot be 0.0 since it goes in the denominator
                      of a_p in the momentum equations */
extern double urfv; /* under-relaxation factor for velocity v. Cannot be 0.0 since it goes in the denominator
                      of a_p in the momentum equations */
extern double urfp; /* under-relaxation factor for pressure */

extern double residual_threshold; /* the value of the residuals to end the iterations */

/* Physical quantities */
extern std::vector<double> position_pressure_node_x; /* position of the pressure nodes in x */
extern std::vector<double> position_pressure_node_y; /* position of the pressure nodes in y */
extern std::vector<double> position_u_velocity_node_x; /* position of the u velocity nodes in x */
extern std::vector<double> position_u_velocity_node_y; /* position of the u velocity nodes in y */
extern std::vector<double> position_v_velocity_node_x; /* position of the v velocity nodes in x */
extern std::vector<double> position_v_velocity_node_y; /* position of the v velocity nodes in y */
extern std::vector<std::vector<double>> Area_pressure_node; /* Cross-section for the pressure nodes */
extern std::vector<std::vector<double>> Area_velocity_node_u; /* Cross-section for the velocity nodes for u velocity with (i,J) indexes */
extern std::vector<std::vector<double>> Area_velocity_node_v; /* Cross-section for the velocity nodes for v velocity with (I,j) indexes */
extern std::vector<std::vector<double>> u_velocity; /* velocity in the x direction */
extern std::vector<std::vector<double>> u_velocity_old; /* the value of the u_velocity in the previous iteration */
extern std::vector<std::vector<double>> v_velocity; /* velocity in the y direction */
extern std::vector<std::vector<double>> v_velocity_old; /* the value of the v_velocity in the previous iteration */
extern std::vector<std::vector<double>> pressure;
extern std::vector<std::vector<double>> pressure_old; /* the value of the pressure in the previous iteration */
extern std::vector<std::vector<double>> d_u; /* parameter d for the pressure correction equation for u velocity with (i,J) indexes */
extern std::vector<std::vector<double>> d_v; /* parameter d for the pressure correction equation for v velocity with (I,j) indexes */

/* Variables used for the iterations */
extern int i_iter; /* number of iterations */
extern int MAX_ITER; /* set the maximum number of iterations to store in the residual vector */
extern std::vector<double> x_momentum_residual_sum; /* sum of the residuals of the x-momentum equation per iteration*/
extern std::vector<double> y_momentum_residual_sum; /* sum of the residuals of the y-momentum equation per iteration*/
extern std::vector<double> pressure_residual_sum; /* sum of the residuals of the pressure equation per iteration*/

/* Modules external to main */
void readinput();
void initialization();
void plotcolormap();
void plotresiduals();
void momentum_equation_solve(MatrixXd &u_star, MatrixXd &v_star, int i_iter);
void pressure_correction_equation_solve(MatrixXd u_star, MatrixXd v_star, MatrixXd &pressure_prime, int i_iter);
void correct_pressure_and_velocities(MatrixXd u_star, MatrixXd v_star, MatrixXd pressure_prime);
void underrelaxation(MatrixXd pressure_prime);
