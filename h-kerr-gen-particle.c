/*
 * Command generation for Kerr spacetime particle simulations
 *
 * Example:  ./h-kerr-gen-particle-dbg 1e-9 4 12 63 1 1 .8 >/tmp/$USER/data
 * Example:  ./h-kerr-gen-particle-dbg 1e-9 12 -1 63 1 1 .8 >/tmp/$USER/data
 *
 * (c) 2018-2022 m4r35n357@gmail.com (Ian Smith), for licencing see the LICENCE file
 */

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>
#include "h-kerr.h"

typedef struct Vector3 { real a, b, c; } vector3;

typedef struct Matrix3x3 { real a, b, c, d, e, f, g, h, i; } matrix3x3;

static real det2x2 (real a, real d, real b, real c) {
    real w = b * c;
    return fmal(a, d, -w) + fmal(-b, c, w);
}

static matrix3x3 m_invert (matrix3x3 m) {
    matrix3x3 c = (matrix3x3){
        .a =  det2x2(m.e, m.i, m.f, m.h), .b = -det2x2(m.d, m.i, m.f, m.g), .c =  det2x2(m.d, m.h, m.e, m.g),
        .d = -det2x2(m.b, m.i, m.c, m.h), .e =  det2x2(m.a, m.i, m.c, m.g), .f = -det2x2(m.a, m.h, m.b, m.g),
        .g =  det2x2(m.b, m.f, m.c, m.e), .h = -det2x2(m.a, m.f, m.c, m.d), .i =  det2x2(m.a, m.e, m.b, m.d)
    };
    real d = m.a * c.a + m.b * c.b + m.c * c.c;
    assert(d != 0.0L);
    return (matrix3x3){
        .a = c.a / d, .b = c.d / d, .c = c.g / d,
        .d = c.b / d, .e = c.e / d, .f = c.h / d,
        .g = c.c / d, .h = c.f / d, .i = c.i / d
    };
}

static vector3 mv_mult (matrix3x3 m, vector3 v) {
    return (vector3){
        .a = m.a * v.a + m.b * v.b + m.c * v.c,
        .b = m.d * v.a + m.e * v.b + m.f * v.c,
        .c = m.g * v.a + m.h * v.b + m.i * v.c
    };
}

static vector3 v_sub (vector3 u, vector3 v) {
    return (vector3){
        .a = u.a - v.a,
        .b = u.b - v.b,
        .c = u.c - v.c
    };
}

static _Bool converged (vector3 v, real epsilon) {
    return fabsl(v.a) < epsilon && fabsl(v.b) < epsilon && fabsl(v.c) < epsilon;
}

static dual R (real r, dual E, dual L, dual Q, real a) {
    real ra2 = r * r + a * a;
    return d_sub(d_sqr(d_sub(d_scale(E, ra2), d_scale(L, a))),
                 d_scale(d_add(d_sqr(d_sub(L, d_scale(E, a))), d_shift(Q, r * r)), ra2 - 2.0L * r));
}

static dual dR_dr (real r, dual E, dual L, dual Q, real a) {
    real ra2 = r * r + a * a;
    return d_sub(d_scale(d_mul(E, d_sub(d_scale(E, ra2), d_scale(L, a))), 4.0L * r),
                 d_shift(d_scale(d_shift(d_add(Q, d_sqr(d_sub(L, d_scale(E, a)))), r * r), 2.0L * r - 2.0L),
                         2.0L * r * (ra2 - 2.0L * r)));
}

static dual THETA (real theta, dual E, dual L, dual Q, real a) {
    real sth2 = sinl(theta) * sinl(theta);
    return d_sub(Q, d_scale(d_add(d_scale(d_shift(d_sqr(E), - 1.0L), - a * a), d_scale(d_sqr(L), 1.0L / sth2)), 1.0L - sth2));
}

static kerr *get_p_gen (char **argv) {
    kerr *p = malloc(sizeof (kerr));
    p->epsilon = strtold(argv[1], NULL);
    p->rmin = strtold(argv[2], NULL);
    p->rmax = strtold(argv[3], NULL);
    p->thmax = elevation_to_colatitude(strtold(argv[4], NULL));
    p->a = strtold(argv[5], NULL);
    p->E = 1.0L;
    p->L = 5.0L;
    p->Q = 0.0L;
    return p;
}

int main (int argc, char **argv) {
    assert(argc == 6);
    kerr *p = get_p_gen(argv);
    matrix3x3 J;
    vector3 x = (vector3){.a = p->E, .b = p->L, .c = p->Q};
    vector3 f = (vector3){.a = 1.0L, .b = 1.0L, .c = 1.0L};
    fprintf(stderr, "\n");
    long count = 0L;
    _Bool circular = p->rmin * p->rmax < 0.0L;
    while (! converged(f, p->epsilon)) {
        J = (matrix3x3){
            .a = R(p->rmin,  d_var(p->E), d_dual(p->L), d_dual(p->Q), p->a).dot,
            .b = R(p->rmin, d_dual(p->E),  d_var(p->L), d_dual(p->Q), p->a).dot,
            .c = R(p->rmin, d_dual(p->E), d_dual(p->L),  d_var(p->Q), p->a).dot,
            .g = THETA(p->thmax,  d_var(p->E), d_dual(p->L), d_dual(p->Q), p->a).dot,
            .h = THETA(p->thmax, d_dual(p->E),  d_var(p->L), d_dual(p->Q), p->a).dot,
            .i = THETA(p->thmax, d_dual(p->E), d_dual(p->L),  d_var(p->Q), p->a).dot
        };
        f = (vector3){
            .a = R(p->rmin, d_dual(p->E), d_dual(p->L), d_dual(p->Q), p->a).val,
            .c = THETA(p->thmax, d_dual(p->E), d_dual(p->L), d_dual(p->Q), p->a).val
        };
        if (! circular) {
            J.d = R(p->rmax,  d_var(p->E), d_dual(p->L), d_dual(p->Q), p->a).dot;
            J.e = R(p->rmax, d_dual(p->E),  d_var(p->L), d_dual(p->Q), p->a).dot;
            J.f = R(p->rmax, d_dual(p->E), d_dual(p->L),  d_var(p->Q), p->a).dot;
            f.b = R(p->rmax, d_dual(p->E), d_dual(p->L), d_dual(p->Q), p->a).val;
        } else {
            J.d = dR_dr(p->rmin,  d_var(p->E), d_dual(p->L), d_dual(p->Q), p->a).dot;
            J.e = dR_dr(p->rmin, d_dual(p->E),  d_var(p->L), d_dual(p->Q), p->a).dot;
            J.f = dR_dr(p->rmin, d_dual(p->E), d_dual(p->L),  d_var(p->Q), p->a).dot;
            f.b = dR_dr(p->rmin, d_dual(p->E), d_dual(p->L), d_dual(p->Q), p->a).val;
        }
        x = v_sub(x, mv_mult(m_invert(J), f));
        fprintf(stderr, "%.18Lf %.18Lf %.18Lf\n", x.a, x.b, x.c);
        p->E = x.a;
        p->L = x.b;
        p->Q = x.c;
        count += 1;
    }
    _Bool valid = 1;
    if (! circular) {
        valid = ! (dR_dr(p->rmin, d_dual(p->E), d_dual(p->L), d_dual(p->Q), p->a).val < 0.0L &&
                   dR_dr(p->rmax, d_dual(p->E), d_dual(p->L), d_dual(p->Q), p->a).val > 0.0L);
    }
    fprintf(stderr, "%.ld iterations, precision %.1Le %s\n",
            count, p->epsilon, valid ? (p->a * p->L < 0.0L ? "RETROGRADE" : "PROGRADE") : "INVALID");
    fprintf(stderr, "\n");
    fprintf(stderr, "Simulate:\n");
    fprintf(stderr, "./h-kerr-dbg 6 8 .01 10000 0 %.3Lf 1.0 %.9Le %.9Le 1.0 %.9Le %.3Lf 0.0 >/tmp/$USER/data\n",
            p->a, p->E, p->L, p->Q, circular ? p->rmin : 0.5L * (p->rmin + p->rmax));
    fprintf(stderr, "\n");
    fprintf(stderr, "./h-kerr-dbg 6 8 .01 10000 0 %.3Lf 1.0 %La %La 1.0 %La %.3Lf 0.0 >/tmp/$USER/data\n",
            p->a, p->E, p->L, p->Q, circular ? p->rmin : 0.5L * (p->rmin + p->rmax));
    fprintf(stderr, "\n");
    fprintf(stderr, "./h-kerr-gl $(yad --columns=2 --title='Kerr Particle Orbit GL' --form --separator=' ' --align=right ");
    fprintf(stderr, "--field='Display Mode':CB --field='Order':NUM --field='Step Size':NUM --field='Steps':NUM ");
    fprintf(stderr, "--field='Track Length':NUM --field='BH spin':NUM --field='particle mass':RO ");
    fprintf(stderr, "--field='energy' --field='momentum' --field='momentum factor' --field='Carter constant' ");
    fprintf(stderr, "--field='r0' --field='theta0' ");
    fprintf(stderr, "-- '0!1!2' '4!2..10!2' '.01!0.001..0.1!0.001!3' '10000!1..1000000!1000' '1000!1..100000!1' "),
    fprintf(stderr, "'%.3Lf!-1.0..1.0!0.1!1' 1.0 %.9Le %.9Le 1.0 %.9Le %.3Lf 0.0)\n",
            p->a, p->E, p->L, p->Q, circular ? p->rmin : 0.5L * (p->rmin + p->rmax));
    fprintf(stderr, "\n");
    fprintf(stderr, "Generate ICs:\n");
    fprintf(stderr, "./h-kerr-dbg 15 8 .01 0 2 %.3Lf 1.0 %La %La 1.0 %La %.3Lf 0.0\n",
            p->a, p->E, p->L, p->Q, circular ? p->rmin : 0.5L * (p->rmin + p->rmax));
    fprintf(stderr, "\n");

    real r_range = (circular ? p->rmin + 1.0L : p->rmax + 1.0L);
    for (int i = 1; i < 1000; i++) {
        real r_plot = r_range * i / 1000.0L;
        real R_plot = R(r_plot, d_dual(p->E), d_dual(p->L), d_dual(p->Q), p->a).val;
        real theta_plot = MY_PI * i / 1000.0L;
        real THETA_plot = THETA(theta_plot, d_dual(p->E), d_dual(p->L), d_dual(p->Q), p->a).val;
        fprintf(stdout, "%.6Lf %.12Lf %.6Lf %.12Lf\n", r_plot, -0.5L * R_plot, theta_plot, -0.5L * THETA_plot);
    }
    return 0;
}
