
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "dual.h"

typedef struct {
    real m;
    real a;
} parameters;

typedef real vector4[4];

typedef real matrix4x4[4][4];

typedef real matrix4x4x4[4][4][4];

void t_input (char **argv, long *dp, long *n, real *h, long *nsteps, real *x, real *y, real *z);

dual ra2 (real a, dual r);

dual delta (real m, real a, dual r);

dual sigma (real a, dual r, dual theta);

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

void real_metric (matrix4x4 metric, real m, real a, real r, real theta);

void real_inverse (matrix4x4 inverse, real m, real a, real r, real theta);

void dg_dr (matrix4x4 dr, real m, real a, real r, real theta);

void dg_dtheta (matrix4x4 dtheta, real m, real a, real r, real theta);

void christoffel (matrix4x4x4 symbols, matrix4x4 inverse, real m, real a, real r, real theta);

//void geodesic (vector4 acclereration, matrix4x4 inverse, vector4 v1, vector4 v2, real m, real a, real r, real theta);

//typedef void *(*tsm_params)(int, char **, long);

//typedef void *(*tsm_inters)(long);

typedef void *(*rk4_params)(int, char **);

//typedef components (*tsm_model)(series, series, series, void *, void *, int);

//typedef void (*rk4_model)(vector4, vector4, vector4, vector4, void *);

void ode (vector4 x_dot, vector4 v_dot, vector4 x, vector4 v, parameters p);

void euler (int argc, char **argv);

//void rk4 (int argc, char **argv, rk4_model ode, rk4_params get_p);
