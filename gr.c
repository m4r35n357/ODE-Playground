/*
 * GR Geodesic equations, solved using Taylor Series integration
 *
 * (c) 2018-2021 m4r35n357@gmail.com (Ian Smith), for licencing see the LICENCE file
 */

#include "gr.h"

series *t_jet4 (long n, vector4 a) {
    series *jet4 = calloc(4, sizeof (series));
    for (int mu = 0; mu < 4; mu++) {
        jet4[mu] = t_jet_c(n + 1, a[mu]);
    }
    return jet4;
}

real r_ra2 (real a, real r) {
    return r * r + a * a;
}

real r_delta (real m, real a, real r) {
    return r * r - 2.0L * m * r + a * a;
}

real r_sigma (real a, real r, real theta) {
    return r * r + a * a * cosl(theta) * cosl(theta);
}

dual d_ra2 (real a, dual r) {
    return d_shift(d_sqr(r), a * a);
}

dual d_delta (real m, real a, dual r) {
    return d_add(d_sqr(r), d_shift(d_scale(r, - 2.0L * m), a * a));
}

dual d_sigma (real a, dual r, dual theta) {
    return d_add(d_sqr(r), d_scale(d_sqr(d_cos(theta)), a * a));
}

void gr_control (char **argv, long *dp, long *n, real *h, long *nsteps,
                real *t, real *r, real *theta, real *phi, real *t_dot, real *r_dot, real *theta_dot, real *phi_dot) {
    *dp = strtol(argv[1], NULL, 10);
    *n = strtol(argv[2], NULL, 10);
    assert(*n >= 1 && *n <= 34);
    *h = strtold(argv[3], NULL);
    assert(*h > 0.0L && *h <= 1.0L);
    *nsteps = strtol(argv[4], NULL, 10);
    assert(*nsteps >= 1 && *nsteps <= 100000);
    *t = strtold(argv[5], NULL);
    *r = strtold(argv[6], NULL);
    *theta = strtold(argv[7], NULL);
    *phi = strtold(argv[8], NULL);
    *t_dot = strtold(argv[9], NULL);
    *r_dot = strtold(argv[10], NULL);
    *theta_dot = strtold(argv[11], NULL);
    *phi_dot = strtold(argv[12], NULL);
}

real mod2_v (series4 x, series4 v, parameters p) {
    matrix4x4 g;
    real_metric (g, d_dual(x[R][0]), d_dual(x[THETA][0]), p);
    real v_dot_v = 0.0L;
    for (int mu = 0; mu < 4; mu++) {
        for (int nu = 0; nu < 4; nu++) {
            v_dot_v += g[mu][nu] * v[mu][0] * v[nu][0];
        }
    }
    return v_dot_v;
}

real error (real e) {
    return 10.0L * log10l(fabsl(e) >= 1e-36L ? fabsl(e) : 1e-36L);
}

void gr_output (long dp, series4 x, series4 v, real t, parameters p) {
    char fs[256];
    real v_dot_v = mod2_v(x, v, p);
    real ev_dot_v = error(p.v0 - v_dot_v);
    real ra = sqrtl(x[R][0] * x[R][0] + p.a * p.a);
    real sth = sinl(x[THETA][0]);
    sprintf(fs, "%%+.%ldLe %%+.%ldLe %%+.%ldLe %%+.6Le %%+.3Le %%+.3Le\n", dp, dp, dp);
    printf(fs, ra * sth * cosl(x[PHI][0]), ra * sth * sinl(x[PHI][0]), x[R][0] * cosl(x[THETA][0]), t, v_dot_v, ev_dot_v);
}

void real_metric (matrix4x4 metric, dual r, dual theta, parameters p) {
    metric[T][T] = g_t_t(p.m, p.a, r, theta).val;
    metric[T][R] = metric[R][T] = g_t_r(p.m, p.a, r, theta).val;
    metric[T][THETA] = metric[THETA][T] = g_t_theta(p.m, p.a, r, theta).val;
    metric[T][PHI] = metric[PHI][T] = g_t_phi(p.m, p.a, r, theta).val;
    metric[R][R] = g_r_r(p.m, p.a, r, theta).val;
    metric[R][THETA] = metric[THETA][R] = g_r_theta(p.m, p.a, r, theta).val;
    metric[R][PHI] = metric[PHI][R] = g_r_phi(p.m, p.a, r, theta).val;
    metric[THETA][THETA] = g_theta_theta(p.m, p.a, r, theta).val;
    metric[THETA][PHI] = metric[PHI][THETA] = g_theta_phi(p.m, p.a, r, theta).val;
    metric[PHI][PHI] = g_phi_phi(p.m, p.a, r, theta).val;
}

void real_inverse (matrix4x4 inverse, real r, real theta, parameters p) {
    inverse[T][T] = i_t_t(p.m, p.a, r, theta);
    inverse[T][R] = inverse[R][T] = i_t_r(p.m, p.a, r, theta);
    inverse[T][THETA] = inverse[THETA][T] = i_t_theta(p.m, p.a, r, theta);
    inverse[T][PHI] = inverse[PHI][T] = i_t_phi(p.m, p.a, r, theta);
    inverse[R][R] = i_r_r(p.m, p.a, r, theta);
    inverse[R][THETA] = inverse[THETA][R] = i_r_theta(p.m, p.a, r, theta);
    inverse[R][PHI] = inverse[PHI][R] = i_r_phi(p.m, p.a, r, theta);
    inverse[THETA][THETA] = i_theta_theta(p.m, p.a, r, theta);
    inverse[THETA][PHI] = inverse[PHI][THETA] = i_theta_phi(p.m, p.a, r, theta);
    inverse[PHI][PHI] = i_phi_phi(p.m, p.a, r, theta);
}

void christoffel (matrix4x4x4 symbols, matrix4x4 inverse, dual r, dual theta, parameters p) {
    matrix4x4 dg_dr;
    r = d_var(r.val);
    dg_dr[T][T] = g_t_t(p.m, p.a, r, theta).dot;
    dg_dr[T][R] = dg_dr[R][T] = g_t_r(p.m, p.a, r, theta).dot;
    dg_dr[T][THETA] = dg_dr[THETA][T] = g_t_theta(p.m, p.a, r, theta).dot;
    dg_dr[T][PHI] = dg_dr[PHI][T] = g_t_phi(p.m, p.a, r, theta).dot;
    dg_dr[R][R] = g_r_r(p.m, p.a, r, theta).dot;
    dg_dr[R][THETA] = dg_dr[THETA][R] = g_r_theta(p.m, p.a, r, theta).dot;
    dg_dr[R][PHI] = dg_dr[PHI][R] = g_r_phi(p.m, p.a, r, theta).dot;
    dg_dr[THETA][THETA] = g_theta_theta(p.m, p.a, r, theta).dot;
    dg_dr[THETA][PHI] = dg_dr[PHI][THETA] = g_theta_phi(p.m, p.a, r, theta).dot;
    dg_dr[PHI][PHI] = g_phi_phi(p.m, p.a, r, theta).dot;
    matrix4x4 dg_dtheta;
    r = d_dual(r.val);
    theta = d_var(theta.val);
    dg_dtheta[T][T] = g_t_t(p.m, p.a, r, theta).dot;
    dg_dtheta[T][R] = dg_dtheta[R][T] = g_t_r(p.m, p.a, r, theta).dot;
    dg_dtheta[T][THETA] = dg_dtheta[THETA][T] = g_t_theta(p.m, p.a, r, theta).dot;
    dg_dtheta[T][PHI] = dg_dtheta[PHI][T] = g_t_phi(p.m, p.a, r, theta).dot;
    dg_dtheta[R][R] = g_r_r(p.m, p.a, r, theta).dot;
    dg_dtheta[R][THETA] = dg_dtheta[THETA][R] = g_r_theta(p.m, p.a, r, theta).dot;
    dg_dtheta[R][PHI] = dg_dtheta[PHI][R] = g_r_phi(p.m, p.a, r, theta).dot;
    dg_dtheta[THETA][THETA] = g_theta_theta(p.m, p.a, r, theta).dot;
    dg_dtheta[THETA][PHI] = dg_dtheta[PHI][THETA] = g_theta_phi(p.m, p.a, r, theta).dot;
    dg_dtheta[PHI][PHI] = g_phi_phi(p.m, p.a, r, theta).dot;
    matrix4x4x4 dg_dx;
    for (int j = 0; j < 4; j++) {
        for (int k = 0; k < 4; k++) {
            dg_dx[j][k][T] = 0.0L;
            dg_dx[j][k][R] = dg_dr[j][k];
            dg_dx[j][k][THETA] = dg_dtheta[j][k];
            dg_dx[j][k][PHI] = 0.0L;
        }
    }
    for (int i = 0; i < 4; i++) {
        for (int k = 0; k < 4; k++) {
            for (int l = 0; l < 4; l++) {
                real sum = 0.0L;
                for (int m = 0; m < 4; m++) {
                    sum += inverse[i][m] * (dg_dx[m][k][l] + dg_dx[m][l][k] - dg_dx[k][l][m]);
                }
                symbols[i][k][l] = 0.5L * sum;
            }
        }
    }
}

void geodesic (matrix4x4x4 symbols, series4 x_dot, series4 v_dot, series4 v, int k) {
    for (int mu = 0; mu < 4; mu++) {
        real sum = 0.0L;
        for (int a = 0; a < 4; a++) {
            for (int b = 0; b < 4; b++) {
                sum -= symbols[mu][a][b] * t_prod(v[a], v[b], k);
            }
        }
        v_dot[mu][k] = sum;
        x_dot[mu][k] = v[mu][k];
    }
}

void tsm4 (int argc, char **argv) {
    (void)argc;
    long n, steps, dp;
    real h;
    vector4 x, v, xdot, vdot;
    gr_control(argv, &dp, &n, &h, &steps, &x[0], &x[1], &x[2], &x[3], &v[0], &v[1], &v[2], &v[3]);
    series4 x4 = t_jet4(n, x), v4 = t_jet4(n, v), xdot4 = t_jet4(n, xdot), vdot4 = t_jet4(n, vdot);
    parameters p = (parameters) {
        .m = strtold(argv[13], NULL),
        .a = strtold(argv[14], NULL)
    };
    p.v0 = mod2_v(x4, v4, p);
    gr_output(dp, x4, v4, 0.0, p);
    for (long step = 1; step <= steps; step++) {
        matrix4x4 inverse;
        real_inverse(inverse, x4[R][0], x4[THETA][0], p);
        matrix4x4x4 symbols;
        christoffel(symbols, inverse, d_dual(x4[R][0]), d_dual(x4[THETA][0]), p);
        for (int k = 0; k < n; k++) {
            geodesic(symbols, xdot4, vdot4, v4, k);
            for (int mu = 0; mu < 4; mu++) {
                x4[mu][k + 1] = xdot4[mu][k] / (k + 1);
                v4[mu][k + 1] = vdot4[mu][k] / (k + 1);
            }
        }
        for (int mu = 0; mu < 4; mu++) {
            *x4[mu] = t_horner(x4[mu], n, h);
            *v4[mu] = t_horner(v4[mu], n, h);
        }
        gr_output(dp, x4, v4, h * step, p);
    }
}
