
#include "gr.h"

void t_control (char **argv, long *dp, long *n, real *h, long *nsteps,
                real *t, real *r, real *theta, real *phi, real *t_dot, real *r_dot, real *theta_dot, real *phi_dot) {
    *dp = strtol(argv[1], NULL, 10);
    *n = strtol(argv[2], NULL, 10);
    assert(*n >= 2 && *n <= 34);
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

void mm_mult (matrix4x4 c, matrix4x4 a, matrix4x4 b) {
    for (int i = 0; i < 4; i++) {
        for (int j = 0; j < 4; j++) {
            real tmp = 0.0L;
            for (int k = 0; k < 4; k++) {
                tmp += a[i][k] * b[k][j];
            }
            c[i][j] = tmp;
        }
    }
}

dual ra2 (real a, dual r) {
    return d_shift(d_sqr(r), a * a);
}

dual delta (real m, real a, dual r) {
    return d_add(d_sqr(r), d_shift(d_scale(r, - 2.0L * m), a * a));
}

dual sigma (real a, dual r, dual theta) {
    return d_add(d_sqr(r), d_scale(d_sqr(d_cos(theta)), a * a));
}

void real_metric (matrix4x4 metric, real m, real a, real r, real theta) {
    metric[0][0] = g_t_t(m, a, d_dual(r), d_dual(theta)).val;
    metric[0][1] = metric[1][0] = g_t_r(m, a, d_dual(r), d_dual(theta)).val;
    metric[0][2] = metric[2][0] = g_t_theta(m, a, d_dual(r), d_dual(theta)).val;
    metric[0][3] = metric[3][0] = g_t_phi(m, a, d_dual(r), d_dual(theta)).val;
    metric[1][1] = g_r_r(m, a, d_dual(r), d_dual(theta)).val;
    metric[1][2] = metric[2][1] = g_r_theta(m, a, d_dual(r), d_dual(theta)).val;
    metric[1][3] = metric[3][1] = g_r_phi(m, a, d_dual(r), d_dual(theta)).val;
    metric[2][2] = g_theta_theta(m, a, d_dual(r), d_dual(theta)).val;
    metric[2][3] = metric[3][2] = g_theta_phi(m, a, d_dual(r), d_dual(theta)).val;
    metric[3][3] = g_phi_phi(m, a, d_dual(r), d_dual(theta)).val;
}

void real_inverse (matrix4x4 inverse, real m, real a, real r, real theta) {
    inverse[0][0] = i_t_t(m, a, r, theta);
    inverse[0][1] = inverse[1][0] = i_t_r(m, a, r, theta);
    inverse[0][2] = inverse[2][0] = i_t_theta(m, a, r, theta);
    inverse[0][3] = inverse[3][0] = i_t_phi(m, a, r, theta);
    inverse[1][1] = i_r_r(m, a, r, theta);
    inverse[1][2] = inverse[2][1] = i_r_theta(m, a, r, theta);
    inverse[1][3] = inverse[3][1] = i_r_phi(m, a, r, theta);
    inverse[2][2] = i_theta_theta(m, a, r, theta);
    inverse[2][3] = inverse[3][2] = i_theta_phi(m, a, r, theta);
    inverse[3][3] = i_phi_phi(m, a, r, theta);
}

void geodesic (vector4 acclereration, matrix4x4 inverse, vector4 v1, vector4 v2, real m, real a, real r, real theta) {
    dual r_dual = d_var(r);
    dual theta_dual = d_dual(theta);
    matrix4x4 dr;
    dr[0][0] = g_t_t(m, a, r_dual, theta_dual).dot;
    dr[0][1] = dr[1][0] = g_t_r(m, a, r_dual, theta_dual).dot;
    dr[0][2] = dr[2][0] = g_t_theta(m, a, r_dual, theta_dual).dot;
    dr[0][3] = dr[3][0] = g_t_phi(m, a, r_dual, theta_dual).dot;
    dr[1][1] = g_r_r(m, a, r_dual, theta_dual).dot;
    dr[1][2] = dr[2][1] = g_r_theta(m, a, r_dual, theta_dual).dot;
    dr[1][3] = dr[3][1] = g_r_phi(m, a, r_dual, theta_dual).dot;
    dr[2][2] = g_theta_theta(m, a, r_dual, theta_dual).dot;
    dr[2][3] = dr[3][2] = g_theta_phi(m, a, r_dual, theta_dual).dot;
    dr[3][3] = g_phi_phi(m, a, r_dual, theta_dual).dot;
    r_dual = d_dual(r);
    theta_dual = d_var(theta);
    matrix4x4 dtheta;
    dtheta[0][0] = g_t_t(m, a, r_dual, theta_dual).dot;
    dtheta[0][1] = dtheta[1][0] = g_t_r(m, a, r_dual, theta_dual).dot;
    dtheta[0][2] = dtheta[2][0] = g_t_theta(m, a, r_dual, theta_dual).dot;
    dtheta[0][3] = dtheta[3][0] = g_t_phi(m, a, r_dual, theta_dual).dot;
    dtheta[1][1] = g_r_r(m, a, r_dual, theta_dual).dot;
    dtheta[1][2] = dtheta[2][1] = g_r_theta(m, a, r_dual, theta_dual).dot;
    dtheta[1][3] = dtheta[3][1] = g_r_phi(m, a, r_dual, theta_dual).dot;
    dtheta[2][2] = g_theta_theta(m, a, r_dual, theta_dual).dot;
    dtheta[2][3] = dtheta[3][2] = g_theta_phi(m, a, r_dual, theta_dual).dot;
    dtheta[3][3] = g_phi_phi(m, a, r_dual, theta_dual).dot;
    matrix4x4x4 derivatives;
    for (int j = 0; j < 4; j++) {
        for (int k = 0; k < 4; k++) {
            derivatives[j][k][0] = 0.0L;
            derivatives[j][k][1] = dr[j][k];
            derivatives[j][k][2] = dtheta[j][k];
            derivatives[j][k][3] = 0.0L;
        }
    }
    matrix4x4x4 symbols;
    for (int i = 0; i < 4; i++) {
        for (int k = 0; k < 4; k++) {
            for (int l = 0; l < 4; l++) {
                real tmp = 0.0L;
                for (int m = 0; m < 4; m++) {
                    tmp += 0.5L * inverse[i][m] * (derivatives[m][k][l] + derivatives[m][l][k] - derivatives[k][l][m]);
                }
                symbols[i][k][l] = tmp;
            }
        }
    }
    for (int mu = 0; mu < 4; mu++) {
        real tmp = 0.0L;
        for (int a = 0; a < 4; a++) {
            for (int b = 0; b < 4; b++) {
                tmp -= symbols[mu][a][b] * v1[a] * v2[b];
            }
        }
        acclereration[mu] = tmp;
    }
}

static void t_output (long dp, vector4 x, vector4 v, real t) {
    char fs[256];
    sprintf(fs, "%%+.%ldLe %%+.%ldLe %%+.%ldLe %%+.%ldLe %%+.%ldLe %%+.%ldLe %%+.%ldLe %%+.%ldLe %%+.6Le\n",
            dp, dp, dp, dp, dp, dp, dp, dp);
    printf(fs, x[0], x[1], x[2], x[3], v[0], v[1], v[2], v[3], t);
}

void ode (vector4 x_dot, vector4 v_dot, vector4 x, vector4 v, parameters p) {
    matrix4x4 inverse;
    real_inverse(inverse, p.m, p.a, x[1], x[2]);
    geodesic(v_dot, inverse, v, v, p.m, p.a, x[1], x[2]);
    for (int i = 0; i < 4; i++) {
        x_dot[i] = v[i];
    }
}

void euler (int argc, char **argv) {
    (void)argc;
    long interval, steps, dp;
    real h;
    vector4 x, v, xdot, vdot;
    t_control(argv, &dp, &interval, &h, &steps, &x[0], &x[1], &x[2], &x[3], &v[0], &v[1], &v[2], &v[3]);
    parameters p = (parameters) {
        .m = strtold(argv[13], NULL),
        .a = strtold(argv[14], NULL)
    };
    t_output(dp, x, v, 0.0);
    for (long step = 1; step < steps + 1; step++) {
        ode(xdot, vdot, x, v, p);
        for (int i = 0; i < 4; i++) {
            x[i] += h * xdot[i];
            v[i] += h * vdot[i];
        }
        if (step % interval == 0) {
            t_output(dp, x, v, h * step);
        }
    }
}
/*
void rk4 (int argc, char **argv, rk4_model ode, rk4_params get_p) {
    long interval, steps, dp;
    real x, y, z, h, xdot = 0.0L, ydot = 0.0L, zdot = 0.0L, kx = 0.0L, ky = 0.0L, kz = 0.0L;
    t_control(argv, &dp, &interval, &h, &steps, &x, &y, &z);
    void *p = get_p(argc, argv);
    t_output(dp, x, y, z, 0.0, "_", "_", "_");
    for (long step = 1; step < steps + 1; step++) {
        components k1 = ode(x, y, z, p);
        components k2 = ode(x + 0.5 * k1.x * h, y + 0.5 * k1.y * h, z + 0.5 * k1.z * h, p);
        components k3 = ode(x + 0.5 * k2.x * h, y + 0.5 * k2.y * h, z + 0.5 * k2.z * h, p);
        components k4 = ode(x + k3.x * h, y + k3.y * h, z + k3.z * h, p);
        x += h * (kx = (k1.x + 2.0 * (k2.x + k3.x) + k4.x) / 6.0);
        y += h * (ky = (k1.y + 2.0 * (k2.y + k3.y) + k4.y) / 6.0);
        z += h * (kz = (k1.z + 2.0 * (k2.z + k3.z) + k4.z) / 6.0);
        if (isnan(x) || isinf(x) || isnan(y) || isinf(y) || isnan(z) || isinf(z)) exit(2);
        if (step % interval == 0) {
            t_output(dp, x, y, z, h * step,
                     kx * xdot < 0.0L ? "X" : "_",
                     ky * ydot < 0.0L ? "Y" : "_",
                     kz * zdot < 0.0L ? "Z" : "_");
            xdot = kx, ydot = ky, zdot = kz;
        }
    }
}
*/
