/*
 * Command generation for Kerr spacetime particle simulations
 *
 * (c) 2018-2023 m4r35n357@gmail.com (Ian Smith), for licencing see the LICENCE file
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "h-kerr.h"

typedef struct Vector3 { real a, b, c; } vector3;

typedef struct Matrix3x3 { real a, b, c, d, e, f, g, h, i; } matrix3x3;

static real det2x2 (real a, real d, real b, real c) {
    real w = b * c;
    return fmal(a, d, -w) + fmal(-b, c, w);
}

static matrix3x3 m_invert (matrix3x3 m) {
    matrix3x3 c = {
        .a =  det2x2(m.e, m.i, m.f, m.h), .b = -det2x2(m.d, m.i, m.f, m.g), .c =  det2x2(m.d, m.h, m.e, m.g),
        .d = -det2x2(m.b, m.i, m.c, m.h), .e =  det2x2(m.a, m.i, m.c, m.g), .f = -det2x2(m.a, m.h, m.b, m.g),
        .g =  det2x2(m.b, m.f, m.c, m.e), .h = -det2x2(m.a, m.f, m.c, m.d), .i =  det2x2(m.a, m.e, m.b, m.d)
    };
    real d = m.a * c.a + m.b * c.b + m.c * c.c;
    CHECK(d != 0.0L);
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

static bool converged (vector3 v, real epsilon) {
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

static parameters *get_p_gen (char **argv) {
    parameters *_ = malloc(sizeof (parameters)); CHECK(_);
    _->epsilon = strtold(argv[1], NULL);
    _->rmin = strtold(argv[2], NULL);
    _->rmax = strtold(argv[3], NULL);
    _->thmax = elevation_to_colatitude(strtold(argv[4], NULL));
    _->a = strtold(argv[5], NULL);
    _->E = 1.0L;
    _->L = 5.0L;
    _->Q = 0.0L;
    return _;
}

int main (int argc, char **argv) { (void)argc;
    PRINT_ARGS(argc, argv);
    CHECK(argc == 6);
    parameters *k = get_p_gen(argv);
    matrix3x3 J;
    vector3 x = {k->E, k->L, k->Q}, f = {1.0L, 1.0L, 1.0L};
    long count = 0L;
    bool circular = k->rmin * k->rmax < 0.0L;
    while (!converged(f, k->epsilon)) {
        J = (matrix3x3){
            .a = R(k->rmin,  d_var(k->E), d_dual(k->L), d_dual(k->Q), k->a).dot,
            .b = R(k->rmin, d_dual(k->E),  d_var(k->L), d_dual(k->Q), k->a).dot,
            .c = R(k->rmin, d_dual(k->E), d_dual(k->L),  d_var(k->Q), k->a).dot,
            .g = THETA(k->thmax,  d_var(k->E), d_dual(k->L), d_dual(k->Q), k->a).dot,
            .h = THETA(k->thmax, d_dual(k->E),  d_var(k->L), d_dual(k->Q), k->a).dot,
            .i = THETA(k->thmax, d_dual(k->E), d_dual(k->L),  d_var(k->Q), k->a).dot
        };
        f = (vector3){
            .a = R(k->rmin, d_dual(k->E), d_dual(k->L), d_dual(k->Q), k->a).val,
            .c = THETA(k->thmax, d_dual(k->E), d_dual(k->L), d_dual(k->Q), k->a).val
        };
        if (!circular) {
            J.d = R(k->rmax,  d_var(k->E), d_dual(k->L), d_dual(k->Q), k->a).dot;
            J.e = R(k->rmax, d_dual(k->E),  d_var(k->L), d_dual(k->Q), k->a).dot;
            J.f = R(k->rmax, d_dual(k->E), d_dual(k->L),  d_var(k->Q), k->a).dot;
            f.b = R(k->rmax, d_dual(k->E), d_dual(k->L), d_dual(k->Q), k->a).val;
        } else {
            J.d = dR_dr(k->rmin,  d_var(k->E), d_dual(k->L), d_dual(k->Q), k->a).dot;
            J.e = dR_dr(k->rmin, d_dual(k->E),  d_var(k->L), d_dual(k->Q), k->a).dot;
            J.f = dR_dr(k->rmin, d_dual(k->E), d_dual(k->L),  d_var(k->Q), k->a).dot;
            f.b = dR_dr(k->rmin, d_dual(k->E), d_dual(k->L), d_dual(k->Q), k->a).val;
        }
        x = v_sub(x, mv_mult(m_invert(J), f));
        fprintf(stderr, "%.18Lf %.18Lf %.18Lf\n", x.a, x.b, x.c);
        k->E = x.a;
        k->L = x.b;
        k->Q = x.c;
        count++;
    }
    bool valid = true;
    if (!circular) {
        valid = ! (dR_dr(k->rmin, d_dual(k->E), d_dual(k->L), d_dual(k->Q), k->a).val < 0.0L &&
                   dR_dr(k->rmax, d_dual(k->E), d_dual(k->L), d_dual(k->Q), k->a).val > 0.0L);
    }
    fprintf(stderr, "%.ld iterations, precision %.1Le %s\n",
            count, k->epsilon, valid ? (k->a * k->L < 0.0L ? "RETROGRADE" : "PROGRADE") : "INVALID");
    fprintf(stderr, "\nSimulate:\n");
    fprintf(stderr, "./h-kerr-std 6 8 .01 10000 0 %.3Lf %.9Le %.9Le 1.0 %.9Le %.3Lf 0.0 >/tmp/$USER/data\n",
            k->a, k->E, k->L, k->Q, circular ? k->rmin : 0.5L * (k->rmin + k->rmax));
    fprintf(stderr, "\n./h-kerr-std 6 8 .01 10000 0 %.3Lf %La %La 1.0 %La %.3Lf 0.0 >/tmp/$USER/data\n",
            k->a, k->E, k->L, k->Q, circular ? k->rmin : 0.5L * (k->rmin + k->rmax));
    fprintf(stderr, "\n./h-kerr-gl $(yad --columns=2 --title='Kerr Particle Orbit GL' --form --separator=' ' --align=right ");
    fprintf(stderr, "--field='Trail Length':NUM --field='Order':NUM --field='Step Size':NUM --field='Steps':NUM ");
    fprintf(stderr, "--field='BH spin':NUM --field='particle mass':RO ");
    fprintf(stderr, "--field='energy' --field='momentum' --field='momentum factor' --field='Carter constant' ");
    fprintf(stderr, "--field='r0' --field='theta0' ");
    fprintf(stderr, "-- '2000!1000..10000!1000' '4!2..10!2' '.01!0.001..0.1!0.001!3' '10000!1..1000000!1000' "),
    fprintf(stderr, "'%.3Lf!-1.0..1.0!0.1!1' 1.0 %.9Le %.9Le 1.0 %.9Le %.3Lf 0.0)\n",
            k->a, k->E, k->L, k->Q, circular ? k->rmin : 0.5L * (k->rmin + k->rmax));

    real r_range = (circular ? k->rmin + 1.0L : k->rmax + 1.0L);
    real PI = acosl(-1.0L);
    for (int i = 1; i < 1000; i++) {
        real r = r_range * i / 1000.0L;
        real theta = PI * i / 1000.0L;
        fprintf(stdout, "%.6Lf %.12Lf %.6Lf %.12Lf\n",
                r, -0.5L * R(r, d_dual(k->E), d_dual(k->L), d_dual(k->Q), k->a).val,
                theta, -0.5L * THETA(theta, d_dual(k->E), d_dual(k->L), d_dual(k->Q), k->a).val);
    }
    return 0;
}
