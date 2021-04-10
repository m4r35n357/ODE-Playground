
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "dual.h"
#include "taylor-ode.h"

typedef enum {T = 0, R = 1, THETA = 2, PHI = 3} coordinate;

typedef struct {
    real m;
    real a;
    real v0;
} parameters;

typedef real vector4[4];

typedef real matrix4x4[4][4];

typedef real matrix4x4x4[4][4][4];

typedef series *series4;

void t_input (char **argv, long *dp, long *n, real *h, long *nsteps, real *x, real *y, real *z);

real r_ra2 (real a, real r);

real r_delta (real m, real a, real r);

real r_sigma (real a, real r, real theta);

dual d_ra2 (real a, dual r);

dual d_delta (real m, real a, dual r);

dual d_sigma (real a, dual r, dual theta);

dual g_t_t (real m, real a, dual r, dual theta);

dual g_t_r (real m, real a, dual r, dual theta);

dual g_t_theta (real m, real a, dual r, dual theta);

dual g_t_phi (real m, real a, dual r, dual theta);

dual g_r_r (real m, real a, dual r, dual theta);

dual g_r_theta (real m, real a, dual r, dual theta);

dual g_r_phi (real m, real a, dual r, dual theta);

dual g_theta_theta (real m, real a, dual r, dual theta);

dual g_theta_phi (real m, real a, dual r, dual theta);

dual g_phi_phi (real m, real a, dual r, dual theta);

real i_t_t (real m, real a, real r, real theta);

real i_t_r (real m, real a, real r, real theta);

real i_t_theta (real m, real a, real r, real theta);

real i_t_phi (real m, real a, real r, real theta);

real i_r_r (real m, real a, real r, real theta);

real i_r_theta (real m, real a, real r, real theta);

real i_r_phi (real m, real a, real r, real theta);

real i_theta_theta (real m, real a, real r, real theta);

real i_theta_phi (real m, real a, real r, real theta);

real i_phi_phi (real m, real a, real r, real theta);

void mm_mult (matrix4x4 c, matrix4x4 a, matrix4x4 b);

void real_metric (matrix4x4 metric, dual r, dual theta, parameters p);

void real_inverse (matrix4x4 inverse, real r, real theta, parameters p);

void christoffel (matrix4x4x4 symbols, matrix4x4 inverse, dual r, dual theta, parameters p);

typedef void *(*rk4_params)(int, char **);

real error (real e);

real mod2_v (series4 x, series4 v, parameters p);

void gr_output (long dp, series4 x, series4 v, real t, parameters p);

void geodesic (matrix4x4x4 symbols, series4 x_dot, series4 v_dot, series4 v, int k);

void euler (int argc, char **argv);

series *t_jet4 (long n, vector4 a);

void tsm4 (int argc, char **argv);
