/*
 * Newtonian central value problem using Hamilton's equations with automatic differentiation
 *
 * (c) 2018-2023 m4r35n357@gmail.com (Ian Smith), for licencing see the LICENCE file
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "symplectic.h"
#include "dual.h"

static dual H (real GM, real m, dual r, dual p_r, dual p_phi) {
    return d_sub(d_scale(d_add(d_sqr(p_r), d_sqr(d_div(p_phi, r))), 0.5L / m), d_scale(d_inv(r), GM * m));
}

struct Parameters {
    real GM, m, r, p_r, phi, p_phi, h0;
};

model *symp_init_p (int argc, char **argv) { (void)argc;
    CHECK(argc == 9);
    model *_ = malloc(sizeof (model)); CHECK(_);
    _->GM = strtold(argv[5], NULL);
    _->m = strtold(argv[6], NULL);
    _->r = strtold(argv[7], NULL);
    _->phi = _->p_r = 0.0L;
    _->p_phi = strtold(argv[8], NULL) * _->m * sqrtl(_->GM * _->r); // arg 8 = 1.0 gives a circular orbit
    _->h0 = H(_->GM, _->m, d_dual(_->r), d_dual(_->p_r), d_dual(_->p_phi)).val;
    return _;
}

void update_q (model *_, real c) {
    _->r   += c * H(_->GM, _->m, d_dual(_->r),  d_var(_->p_r), d_dual(_->p_phi)).dot;
    _->phi += c * H(_->GM, _->m, d_dual(_->r), d_dual(_->p_r),  d_var(_->p_phi)).dot;
}

void update_p (model *_, real d) {
    _->p_r -= d * H(_->GM, _->m, d_var(_->r), d_dual(_->p_r), d_dual(_->p_phi)).dot;
}

static void plot (int dp, model *_, real t) {
    real h = H(_->GM, _->m, d_dual(_->r), d_dual(_->p_r), d_dual(_->p_phi)).val;
    printf("% .*Le % .*Le % .*Le % .*Le %.6Le % .*Le % .*Le\n",
           dp, _->r * sinl(_->phi), dp, _->r * cosl(_->phi), dp, _->r, dp, _->p_r, t, dp, error(h - _->h0), dp, h);
}

int main (int argc, char **argv) {
    solve(symp_get_c(argc, argv), symp_init_p(argc, argv), plot);
    return 0;
}
