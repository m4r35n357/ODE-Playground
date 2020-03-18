/*
 * Double Pendulum System
 *
 * Example: ./tsm-double-pendulum-dbg 16 10 0.01 10000 1 1 1 1 3 -1 3 -1
 *
 * (c) 2018-2020 m4r35n357@gmail.com (Ian Smith), for licencing see the LICENCE file
 */

#include <assert.h>
#include <mpfr.h>
#include "taylor-ode.h"

static void polar_to_rectangular (mpfr_t *xa, mpfr_t *ya, mpfr_t *xb, mpfr_t *yb, mpfr_t la, mpfr_t lb, mpfr_t tt, mpfr_t ta, mpfr_t tb, mpfr_t wa, mpfr_t wb, mpfr_t *s, mpfr_t *c) {
    mpfr_sin_cos(*s, *c, ta, RND);
    mpfr_mul(*xa, la, *s, RND);
    mpfr_mul(*ya, la, *c, RND);
    mpfr_neg(*ya, *ya, RND);
    mpfr_sin_cos(*s, *c, tb, RND);
    mpfr_fma(*xb, lb, *s, *xa, RND);
    mpfr_fms(*yb, lb, *c, *ya, RND);
    mpfr_neg(*yb, *yb, RND);
    mpfr_printf("%+.9RNe %+.9RNe %+.9RNe %+.9RNe %+.5RNe %+.9RNe %+.9RNe %+.9RNe %+.9RNe\n", xa, ya, xb, yb, tt, ta, tb, wa, wb);
}

int main (int argc, char **argv) {
    long n, nsteps;
    mpfr_t t, t10, w10, t20, w20, g, m1, m2, l1, l2, h, _, __, x1, y1, x2, y2, *d, *st1, *ct1, *w1_2, *w2_2;
    mpfr_t *_t1_t2, *_t1_2t2, *_2t1_2t2, *st1_t2, *ct1_t2, *st1_2t2, *ct1_2t2, *s2t1_2t2, *c2t1_2t2, *n1, *n1_, *n2, *n2_, *q1, *q2;

    // initialize from command arguments
    assert(argc == 13);
    t_stepper(argv, &n, &t, &h, &nsteps);
    t_args(argv, argc, &m1, &m2, &l1, &l2, &t10, &w10, &t20, &w20);
    mpfr_inits(_, __, x1, y1, x2, y2, NULL);
    mpfr_init_set_str(g, "9.80665", 10, RND);

    // initialize the derivative and temporary jets
    mpfr_t *t1 = t_jet_c(n + 1, t10);
    mpfr_t *w1 = t_jet_c(n + 1, w10);
    mpfr_t *t2 = t_jet_c(n + 1, t20);
    mpfr_t *w2 = t_jet_c(n + 1, w20);
    st1 = t_jet(n);
    ct1 = t_jet(n);
    _t1_t2 = t_jet(n);
    _t1_2t2 = t_jet(n);
    _2t1_2t2 = t_jet(n);
    st1_t2 = t_jet(n);
    ct1_t2 = t_jet(n);
    st1_2t2 = t_jet(n);
    ct1_2t2 = t_jet(n);
    s2t1_2t2 = t_jet(n);
    c2t1_2t2 = t_jet(n);
    w1_2 = t_jet(n);
    w2_2 = t_jet(n);
    n1 = t_jet(n);
    n2 = t_jet(n);
    n1_ = t_jet(n);
    n2_ = t_jet(n);
    q1 = t_jet(n);
    q2 = t_jet(n);
    d = t_jet(n);

    // main loop
    polar_to_rectangular(&x1, &y1, &x2, &y2, l1, l2, t, t1[0], t2[0], w1[0], w2[0], &_, &__);
    for (long step = 1; step < nsteps + 1; step++) {
        // compute the taylor coefficients
        for (int k = 0; k < n; k++) {
            t_sin_cos(st1, ct1, t1, k, TRIG);

            mpfr_sub(_t1_t2[k], t1[k], t2[k], RND);
            t_sin_cos(st1_t2, ct1_t2, _t1_t2, k, TRIG);

            mpfr_mul_2ui(__, t2[k], 1, RND);
            mpfr_sub(_t1_2t2[k], t1[k], __, RND);
            t_sin_cos(st1_2t2, ct1_2t2, _t1_2t2, k, TRIG);

            mpfr_mul_2ui(_, t1[k], 1, RND);
            mpfr_sub(_2t1_2t2[k], _, __, RND);
            t_sin_cos(s2t1_2t2, c2t1_2t2, _2t1_2t2, k, TRIG);

            mpfr_set(w1_2[k], *t_sqr(w1, k), RND);
            mpfr_set(w2_2[k], *t_sqr(w2, k), RND);

            // TODO l & m set to 1 for now
            mpfr_add(n1_[k], w2_2[k], *t_prod(w1_2, ct1_t2, k), RND);
            mpfr_mul_2ui(_, *t_prod(st1_t2, n1_, k), 1, RND);
            mpfr_fma(_, g, st1_2t2[k], _, RND);
            mpfr_mul_ui(__, g, 3, RND);
            mpfr_fma(n1[k], __, st1[k], _, RND);
            mpfr_neg(n1[k], n1[k], RND);

            // TODO l & m set to 1 for now
            mpfr_fma(__, g, ct1[k], w1_2[k], RND);
            mpfr_mul_2ui(__, __, 1, RND);
            mpfr_add(n2_[k], __, *t_prod(w2_2, ct1_t2, k), RND);
            mpfr_mul_2ui(n2[k], *t_prod(st1_t2, n2_, k), 1, RND);

            mpfr_mul_2ui(_, m1, 1, RND);
            mpfr_fms(__, c2t1_2t2[k], m2, m2, RND);
            mpfr_sub(d[k], _, __, RND);

            //  t1' = w1
            mpfr_div_ui(t1[k + 1], w1[k], k + 1, RND);
            //  t2' = w2
            mpfr_div_ui(t2[k + 1], w2[k], k + 1, RND);
            //  w1' = dfgsdfg
            mpfr_div(_, *t_quot(q1, n1, d, k), l1, RND);
            mpfr_div_ui(w1[k + 1], _, k + 1, RND);
            //  w2' = fgdfghdfg
            mpfr_div(_, *t_quot(q2, n2, d, k), l2, RND);
            mpfr_div_ui(w2[k + 1], _, k + 1, RND);
        }

        // sum the series using Horner's method and advance one step
        t_horner(&t10, t1, n, h);
        t_horner(&t20, t2, n, h);
        t_horner(&w10, w1, n, h);
        t_horner(&w20, w2, n, h);
        mpfr_mul_ui(t, h, step, RND);
        polar_to_rectangular(&x1, &y1, &x2, &y2, l1, l2, t, t1[0], t2[0], w1[0], w2[0], &_, &__);
    }
    return 0;
}
