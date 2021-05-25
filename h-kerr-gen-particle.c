
/*
 * Parameter (E, L, Q) generation for Kerr spacetime
 *
 * Example:  ./h-kerr-gen-particle-dbg 1e-9 3 12 63 1 1 .8 >/tmp/$USER/data
 * Example:  ./h-kerr-gen-particle-dbg 1e-9 12 -1 63 1 1 .8 >/tmp/$USER/data
 *
 * Potential plots:
 *
 * Example:  gnuplot -p -e "set terminal wxt size 600,450; set yrange [*:10]; plot 0.0, '/tmp/$USER/data' using 1:2 with lines"
 * Example:  gnuplot -p -e "set terminal wxt size 600,450; set yrange [*:1]; plot 0.0, '/tmp/$USER/data' using 3:4 with lines"
 *
 * (c) 2018-2021 m4r35n357@gmail.com (Ian Smith), for licencing see the LICENCE file
 */

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>
#include "dual.h"

typedef struct {
    real a; real b; real c;
} vector3;

typedef struct {
    real a; real b; real c;
    real d; real e; real f;
    real g; real h; real i;
} matrix3x3;

static matrix3x3 invert (matrix3x3 m) {
    matrix3x3 c = (matrix3x3) {
        .a =  (m.e * m.i - m.f * m.h), .b = -(m.d * m.i - m.f * m.g), .c =  (m.d * m.h - m.e * m.g),
        .d = -(m.b * m.i - m.c * m.h), .e =  (m.a * m.i - m.c * m.g), .f = -(m.a * m.h - m.b * m.g),
        .g =  (m.b * m.f - m.c * m.e), .h = -(m.a * m.f - m.c * m.d), .i =  (m.a * m.e - m.b * m.d)
    };
    real d = m.a * c.a + m.b * c.b + m.c * c.c;
    assert(d != 0.0L);
    return (matrix3x3) {
        .a = c.a / d, .b = c.d / d, .c = c.g / d,
        .d = c.b / d, .e = c.e / d, .f = c.h / d,
        .g = c.c / d, .h = c.f / d, .i = c.i / d
    };
}

static vector3 mv_mult (matrix3x3 m, vector3 v) {
    return (vector3) {
        .a = m.a * v.a + m.b * v.b + m.c * v.c,
        .b = m.d * v.a + m.e * v.b + m.f * v.c,
        .c = m.g * v.a + m.h * v.b + m.i * v.c
    };
}

static vector3 v_sub (vector3 u, vector3 v) {
    return (vector3) {
        .a = u.a - v.a,
        .b = u.b - v.b,
        .c = u.c - v.c
    };
}

static _Bool converged (vector3 v, real epsilon) {
    return fabsl(v.a) < epsilon && fabsl(v.b) < epsilon && fabsl(v.c) < epsilon;
}

static dual R (real r, dual E, dual L, dual Q, real M, real a, real mu2) {
    real ra2 = r * r + a * a;
    return d_sub(d_sqr(d_sub(d_scale(E, ra2), d_scale(L, a))),
                 d_scale(d_add(d_sqr(d_sub(L, d_scale(E, a))), d_shift(Q, mu2 * r * r)), ra2 - 2.0L * M * r));
}

static dual dR_dr (real r, dual E, dual L, dual Q, real M, real a, real mu2) {
    real ra2 = r * r + a * a;
    return d_sub(d_scale(d_mul(E, d_sub(d_scale(E, ra2), d_scale(L, a))), 4.0L * r),
                 d_shift(d_scale(d_shift(d_add(Q, d_sqr(d_sub(L, d_scale(E, a)))), mu2 * r * r), 2.0L * r - 2.0L * M),
                         2.0L * mu2 * r * (ra2 - 2.0L * M * r)));
}

static dual THETA (real theta, dual E, dual L, dual Q, real mu2, real a) {
    real sth2 = sinl(theta) * sinl(theta);
    return d_sub(Q, d_scale(d_add(d_scale(d_shift(d_sqr(E), - mu2), - a * a), d_scale(d_sqr(L), 1.0L / sth2)), 1.0L - sth2));
}

typedef struct {
    real epsilon;  // precision
    real bh_mass;  // central mass
    real pmass2;  // particle mass (squared)
    real E, L, Q;  // constants of motion
    real spin;  // global constants
    real rmin, rmax, thmax;  // constraints
} parameters;

static parameters *get_p (char **argv) {
    parameters *p = malloc(sizeof (parameters));
    p->epsilon = strtold(argv[1], NULL);
    p->rmin = strtold(argv[2], NULL);
    p->rmax = strtold(argv[3], NULL);
    p->thmax = elevation_to_colatitude(strtold(argv[4], NULL));
    p->bh_mass = strtold(argv[5], NULL);
    p->pmass2 = strtold(argv[6], NULL);
    p->spin = strtold(argv[7], NULL);
    p->E = 1.0L;
    p->L = 5.0L;
    p->Q = 0.0L;
    return p;
}

int main (int argc, char **argv) {
    assert(argc == 8);
    parameters *p = get_p(argv);
    matrix3x3 J;
    vector3 x = (vector3) {
        .a = p->E,
        .b = p->L,
        .c = p->Q
    };
    vector3 f = (vector3) {
        .a = 1.0L,
        .b = 1.0L,
        .c = 1.0L
    };
    fprintf(stderr, "\n");
    long count = 0L;
    _Bool circular = p->rmin * p->rmax < 0.0L;
    while (! converged(f, p->epsilon)) {
        if (! circular) {
            J = (matrix3x3) {
                .a = R(p->rmin,  d_var(p->E), d_dual(p->L), d_dual(p->Q), p->bh_mass, p->spin, p->pmass2).dot,
                .b = R(p->rmin, d_dual(p->E),  d_var(p->L), d_dual(p->Q), p->bh_mass, p->spin, p->pmass2).dot,
                .c = R(p->rmin, d_dual(p->E), d_dual(p->L),  d_var(p->Q), p->bh_mass, p->spin, p->pmass2).dot,
                .d = R(p->rmax,  d_var(p->E), d_dual(p->L), d_dual(p->Q), p->bh_mass, p->spin, p->pmass2).dot,
                .e = R(p->rmax, d_dual(p->E),  d_var(p->L), d_dual(p->Q), p->bh_mass, p->spin, p->pmass2).dot,
                .f = R(p->rmax, d_dual(p->E), d_dual(p->L),  d_var(p->Q), p->bh_mass, p->spin, p->pmass2).dot,
                .g = THETA(p->thmax,  d_var(p->E), d_dual(p->L), d_dual(p->Q), p->pmass2, p->spin).dot,
                .h = THETA(p->thmax, d_dual(p->E),  d_var(p->L), d_dual(p->Q), p->pmass2, p->spin).dot,
                .i = THETA(p->thmax, d_dual(p->E), d_dual(p->L),  d_var(p->Q), p->pmass2, p->spin).dot
            };
            f = (vector3) {
                .a = R(p->rmin, d_dual(p->E), d_dual(p->L), d_dual(p->Q), p->bh_mass, p->spin, p->pmass2).val,
                .b = R(p->rmax, d_dual(p->E), d_dual(p->L), d_dual(p->Q), p->bh_mass, p->spin, p->pmass2).val,
                .c = THETA(p->thmax, d_dual(p->E), d_dual(p->L), d_dual(p->Q), p->pmass2, p->spin).val
            };
        } else {
            J = (matrix3x3) {
                .a = R(p->rmin,  d_var(p->E), d_dual(p->L), d_dual(p->Q), p->bh_mass, p->spin, p->pmass2).dot,
                .b = R(p->rmin, d_dual(p->E),  d_var(p->L), d_dual(p->Q), p->bh_mass, p->spin, p->pmass2).dot,
                .c = R(p->rmin, d_dual(p->E), d_dual(p->L),  d_var(p->Q), p->bh_mass, p->spin, p->pmass2).dot,
                .d = dR_dr(p->rmin,  d_var(p->E), d_dual(p->L), d_dual(p->Q), p->bh_mass, p->spin, p->pmass2).dot,
                .e = dR_dr(p->rmin, d_dual(p->E),  d_var(p->L), d_dual(p->Q), p->bh_mass, p->spin, p->pmass2).dot,
                .f = dR_dr(p->rmin, d_dual(p->E), d_dual(p->L),  d_var(p->Q), p->bh_mass, p->spin, p->pmass2).dot,
                .g = THETA(p->thmax,  d_var(p->E), d_dual(p->L), d_dual(p->Q), p->pmass2, p->spin).dot,
                .h = THETA(p->thmax, d_dual(p->E),  d_var(p->L), d_dual(p->Q), p->pmass2, p->spin).dot,
                .i = THETA(p->thmax, d_dual(p->E), d_dual(p->L),  d_var(p->Q), p->pmass2, p->spin).dot
            };
            f = (vector3) {
                .a = R(p->rmin, d_dual(p->E), d_dual(p->L), d_dual(p->Q), p->bh_mass, p->spin, p->pmass2).val,
                .b = dR_dr(p->rmin, d_dual(p->E), d_dual(p->L), d_dual(p->Q), p->bh_mass, p->spin, p->pmass2).val,
                .c = THETA(p->thmax, d_dual(p->E), d_dual(p->L), d_dual(p->Q), p->pmass2, p->spin).val
            };
        }
        x = v_sub(x, mv_mult(invert(J), f));
        fprintf(stderr, "%.18Lf %.18Lf %.18Lf\n", x.a, x.b, x.c);
        p->E = x.a;
        p->L = x.b;
        p->Q = x.c;
        count += 1;
    }
    _Bool valid = 1;
    if (! circular) {
        valid = ! (dR_dr(p->rmin, d_dual(p->E), d_dual(p->L), d_dual(p->Q), p->bh_mass, p->spin, p->pmass2).val < 0.0L &&
                   dR_dr(p->rmax, d_dual(p->E), d_dual(p->L), d_dual(p->Q), p->bh_mass, p->spin, p->pmass2).val > 0.0L);
    }
    fprintf(stderr, "%.ld iterations, precision %.1Le %s\n",
            count, p->epsilon, valid ? (p->spin * p->L < 0.0L ? "RETROGRADE" : "PROGRADE") : "INVALID");
    fprintf(stderr, "\n");
    fprintf(stderr, "Simulate:\n");
    fprintf(stderr, "./h-kerr-sd-dbg 6 8 .01 10000 0 %.3Lf 1.0 1.0 %La %La 1.0 %La %.3Lf %.3Lf >/tmp/$USER/data\n",
            p->spin, p->E, p->L, p->Q, circular ? p->rmin : 0.5L * (p->rmin + p->rmax), 0.0L);
    fprintf(stderr, "\n");
    fprintf(stderr, "Generate ICs:\n");
    fprintf(stderr, "./h-kerr-sd-dbg 15 8 .01 0 2 %.3Lf 1.0 1.0 %La %La 1.0 %La %.3Lf %.3Lf\n",
            p->spin, p->E, p->L, p->Q, circular ? p->rmin : 0.5L * (p->rmin + p->rmax), 0.0L);
    fprintf(stderr, "\n");

    real r_range = (circular ? p->rmin + 1.0L : p->rmax + 1.0L);
    real theta_range = get_PI();
    for (int i = 1; i < 1000; i += 1) {
        real r_plot = r_range * i / 1000.0L;
        real R_plot = R(r_plot, d_dual(p->E), d_dual(p->L), d_dual(p->Q), p->bh_mass, p->spin, p->pmass2).val;
        real theta_plot = theta_range * i / 1000.0L;
        real THETA_plot = THETA(theta_plot, d_dual(p->E), d_dual(p->L), d_dual(p->Q), p->pmass2, p->spin).val;
        fprintf(stdout, "%.6Lf %.12Lf %.6Lf %.12Lf\n", r_plot, -0.5L * R_plot, theta_plot, -0.5L * THETA_plot);
    }
    return 0;
}
