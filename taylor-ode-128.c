/*
 * Compile using c99 -g -Og -o tsm-lorenz-128 taylor-ode-128.c -lquadmath
 * or, using Clang:
 * clang -Wall -g -Og -o tsm-lorenz-128 taylor-ode-128.c -lquadmath -I'/usr/lib/gcc/x86_64-linux-gnu/7/include'
 */

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <stdarg.h>
#include <quadmath.h>

#define KNRM "\x1B[0;37m"
#define KCYN "\x1B[36m"
#define KGRY "\x1B[2;37m"
#define KBLD "\x1B[1;37m"

typedef __float128 quad;

typedef struct { quad *x, *y; } jet_pair;

typedef struct { quad x, y; } quad_pair;

typedef enum {TRIG, HYP} geometry;

typedef enum {VARIABLE, CONSTANT} ad_status;

long order = 10, nsteps = 1001, n;

quad t, x, y, z, s, r, b, h, tmp, *cx, *cy, *cz, *cx0, *cx1, *cx9, *PI_3, *PI_4, *we, *wl, *wsq, *wsum, *wprod, *wquot, *wpwr;

jet_pair ts_pair;

void t_line_output (quad t, int count, ...) {
    char buf[18];
    va_list data;
    va_start(data, count);
    for (int i = 0; i < count; i++) {
        quadmath_snprintf(buf, sizeof buf, "%.34Qf", va_arg(data, quad));
        printf("%s ", buf);
    }
    va_end(data);
    quadmath_snprintf(buf, sizeof buf, "%.34Qf", t);
    printf("%s\n", buf);
}

void t_stepper (int argc, char **argv, long *n, quad *t, quad *h, long *nsteps) {
    *n = strtol(argv[2], NULL, 10);
    assert(*n > 1);
    *t = 0.0q;
    *h = strtoflt128(argv[3], NULL);
    *nsteps = strtol(argv[4], NULL, 10);
}

quad *t_jet (int n) {
    return malloc(sizeof (quad) * n);
}

quad *t_jet_constant (int n, quad value) {
    assert(n > 0);
    quad *jet = t_jet(n);
    jet[0] = value;
    for (int i = 1; i < n; i++) {
         jet[i] = 0.0q;
    }
    return jet;
}

quad t_horner (const quad *jet, long n, quad h) {
    assert(n > 0);
    quad sum = jet[n];
    for (int j = n - 1; j > - 1; j--) {
        sum = sum * h + jet[j];
    }
    return sum;
}

quad t_square (const quad *U, int k) {
    assert(k >= 0);
    if (k == 0) {
        return U[0] * U[0];
    } else {
        quad S = 0.0q;
        if (k % 2 == 1) {
            for (int j = 0; j < (k - 1) / 2 + 1; j++) {
                S += U[j] * U[k - j];
            }
            return 2.0q * S;
        } else {
            for (int j = 0; j < (k - 2) / 2 + 1; j++) {
                S += U[j] * U[k - j];
            }
            return 2.0q * S + U[k / 2] * U[k / 2];
        }
    }
}

quad t_product (const quad *V, const quad *U, int k) {
    assert(sizeof *U == sizeof *V);
    assert(k >= 0);
    quad P = 0.0q;
    for (int j = 0; j < k + 1; j++) {
        P += U[j] * V[k - j];
    }
    return P;
}

quad t_quotient (const quad *Q, const quad *U, const quad *V, int k) {
    assert(V[0] != 0.0q);
    assert(Q != U && Q != V && U != V);
    assert(sizeof *U == sizeof *Q && sizeof *V == sizeof *Q);
    assert(k >= 0);
    quad tmp = 0.0q;
    for (int j = 1; j < k + 1; j++) {
        tmp += V[j] * Q[k - j];
    }
    return (U[k] - tmp) / V[0];
}

static quad ddot (const quad *V, const quad *U, int k) {
    assert(sizeof *U == sizeof *V);
    quad dd = 0.0q;
    for (int j = 1; j < k; j++) {
        dd += j * U[j] * V[k - j];
    }
    return dd / k;
}

quad t_exp (const quad *E, const quad *U, int k) {
    assert(E != U);
    assert(sizeof *U == sizeof *E);
    assert(k >= 0);
    if (k == 0) {
        return expq(U[0]);
    } else {
        return E[0] * U[k] + ddot(E, U, k);
    }
}

quad_pair t_sin_cos (const quad *S, const quad *C, const quad *U, int k, geometry g) {
    assert(S != C && S != U && C != U);
    assert(sizeof *U == sizeof *S && sizeof *U == sizeof *C);
    assert(k >= 0);
    if (k == 0) {
        if (g == TRIG) {
            quad sin, cos;
            sincosq(U[0], &sin, &cos);
            return (quad_pair) {sin, cos};
        } else {
            return (quad_pair) {sinhq(U[0]), coshq(U[0])};
        }
    } else {
        quad sin = C[0] * U[k] + ddot(C, U, k);
        quad cos = S[0] * U[k] + ddot(S, U, k);
        if (g == TRIG) {
            return (quad_pair) {sin, - cos};
        } else {
            return (quad_pair) {sin, cos};
        }
    }
}

quad_pair t_tan_sec2 (const quad *T, const quad *S2, const quad *U, int k, geometry g) {
    assert(T != S2 && T != U && S2 != U);
    assert(sizeof *U == sizeof *T && sizeof *U == sizeof *S2);
    assert(k >= 0);
    if (k == 0) {
        if (g == TRIG) {
            quad tan = tanq(U[0]);
            return (quad_pair) {tan, 1.0q + tan * tan};
        } else {
            quad tanh = tanhq(U[0]);
            return (quad_pair) {tanh, 1.0q - tanh * tanh};
        }
    } else {
        quad tan = S2[0] * U[k] + ddot(S2, U, k);
        quad sec2 = 2.0q * (T[0] * tan + ddot(T, T, k));
        if (g == TRIG) {
            return (quad_pair) {tan, sec2};
        } else {
            return (quad_pair) {tan, - sec2};
        }
    }
}

quad t_power (const quad *P, const quad *U, const quad a, int k) {
    assert(U[0] > 0.0q);
    assert(P != U);
    assert(sizeof *U == sizeof *P);
    assert(k >= 0);
    if (k == 0) {
        return powq(U[0], a);
    } else {
        return (a * (P[0] * U[k] + ddot(P, U, k)) - ddot(U, P, k)) / U[0];
    }
}

void set_ad_status (quad *jet, ad_status status) {
    if (status == VARIABLE) {
        jet[1] = 1.0q;
    } else {
        jet[1] = 0.0q;
    }
}

quad *ad_scale (const quad *U, quad a, int n) {
    quad *scaled = t_jet(n);
    for (int k = 0; k < n; k++) {
        scaled[k] = a * U[k];
    }
    return scaled;
}

quad *ad_plus (const quad *U, const quad *V, int n) {
    quad *sum = t_jet(n);
    for (int k = 0; k < n; k++) {
        sum[k] = U[k] + V[k];
    }
    return sum;
}

quad *ad_minus (const quad *U, const quad *V, int n) {
    quad *difference = t_jet(n);
    for (int k = 0; k < n; k++) {
        difference[k] = U[k] - V[k];
    }
    return difference;
}

quad *ad_square (const quad *U, int n) {
    quad *square = t_jet(n);
    for (int k = 0; k < n; k++) {
        square[k] = t_square(U, k);
    }
    return square;
}

quad *ad_product (const quad *V, const quad *U, int n) {
    quad *product = t_jet(n);
    for (int k = 0; k < n; k++) {
        product[k] = t_product(V, U, k);
    }
    return product;
}

quad *ad_quotient (const quad *U, const quad *V, int n) {
    quad *quotient = t_jet(n);
    for (int k = 0; k < n; k++) {
        quotient[k] = t_quotient(quotient, U, V, k);
    }
    return quotient;
}

quad *ad_exp (const quad *U, int n) {
    quad *exp = t_jet(n);
    for (int k = 0; k < n; k++) {
        exp[k] = t_exp(exp, U, k);
    }
    return exp;
}

jet_pair ad_sin_cos (const quad *U, int n, geometry g) {
    quad *sin = t_jet(n);
    quad *cos = t_jet(n);
    for (int k = 0; k < n; k++) {
        quad_pair result = t_sin_cos(sin, cos, U, k, g);
        sin[k] = result.x;
        cos[k] = result.y;
    }
    return (jet_pair) {sin, cos};
}

jet_pair ad_tan_sec2 (const quad *U, int n, geometry g) {
    quad *tan = t_jet(n);
    quad *sec2 = t_jet(n);
    for (int k = 0; k < n; k++) {
        quad_pair result = t_tan_sec2(tan, sec2, U, k, g);
        tan[k] = result.x;
        sec2[k] = result.y;
    }
    return (jet_pair) {tan, sec2};
}

quad *ad_power (const quad *U, quad a, int n) {
    quad *power = t_jet(n);
    for (int k = 0; k < n; k++) {
        power[k] = t_power(power, U, a, k);
    }
    return power;
}

void jet_output (const quad *jet, long n, char* f_colour, char *fk_colour) {
    char buf[15];
    quadmath_snprintf(buf, sizeof buf, "%9.6Qf", jet[0]);
    printf("%s%s", f_colour, buf);
    for (int i = 1; i < n; i++) {
        quadmath_snprintf(buf, sizeof buf, "%9.6Qf", jet[i]);
        printf(" %s%s", fk_colour, buf);
    }
    printf("%s\n", f_colour);
}

void jet_to_derivs (quad *jet, long n) {
    for (int i = 1, fac = 1; i < n; i++) {
        fac *= i;
        jet[i] *= fac;
    }
}

void derivative_output (quad *jet, long n, char* f_colour, char *fk_colour) {
    jet_to_derivs(jet, n);
    jet_output(jet, n, f_colour, fk_colour);
}

int ad_test (int argc, char **argv) {
    assert(argc == 4);

    n = strtol(argv[1], NULL, 10);
    assert(n > 1);
    x = strtoflt128(argv[2], NULL);
    y = strtoflt128(argv[3], NULL);

    cx = t_jet_constant(n, x);
    cy = t_jet_constant(n, y);

    cx0 = t_jet_constant(n, 0.0q);
    cx0[1] = 1.0q;

    cx1 = t_jet_constant(n, 1.0q);

    cx9 = t_jet_constant(n, 9.0q);
    cx9[1] = 1.0q;

    PI_3 = t_jet_constant(n, M_PIq / 3.0);
    PI_3[1] = 1.0q;

    PI_4 = t_jet_constant(n, M_PI_4q);
    PI_4[1] = 1.0q;

    wsq = t_jet(n);
    wsum = t_jet(n);
    wprod = t_jet(n);
    wpwr = t_jet(n);
    wl = t_jet(n);

    printf("\n%sx = %s, y = %s, order = %ld%s\n\n", KBLD, argv[2], argv[3], n - 1, KNRM);

    set_ad_status(cx, VARIABLE);
    set_ad_status(cy, CONSTANT);

    printf("%s%s%s\n", KCYN, "f(x) = x", KNRM);
    jet_output(cx, n, KNRM, KGRY);
    derivative_output(cx, n, KBLD, KGRY);
    printf("%s\n", KNRM);

    printf("%s%s%s\n", KCYN, "f(x) = x^2", KNRM);
    wprod = ad_product(cx, cx, n);
    jet_output(wprod, n, KNRM, KGRY);
    wsq = ad_square(cx, n);
    jet_output(wsq, n, KNRM, KGRY);
    derivative_output(wsq, n, KBLD, KGRY);
    printf("%s\n", KNRM);

    printf("%s%s%s\n", KCYN, "f(x) = x^3", KNRM);
    wsq = ad_square(cx, n);
    wprod = ad_product(wsq, cx, n);
    jet_output(wprod, n, KNRM, KGRY);
    derivative_output(wprod, n, KBLD, KGRY);
    printf("%s\n", KNRM);

    printf("%s%s%s\n", KCYN, "f(x) = x^4", KNRM);
    wsq = ad_square(cx, n);
    wprod = ad_product(wsq, wsq, n);
    jet_output(wprod, n, KNRM, KGRY);
    wprod = ad_product(cx, cx, n);
    wsq = ad_square(wprod, n);
    jet_output(wsq, n, KNRM, KGRY);
    derivative_output(wsq, n, KBLD, KGRY);
    printf("%s\n", KNRM);

    printf("%s%s%s\n", KCYN, "f(x) = 1 / x", KNRM);
    wpwr = ad_power(cx, - 1.0q, n);
    jet_output(wpwr, n, KNRM, KGRY);
    derivative_output(wpwr, n, KBLD, KGRY);
    wquot = ad_quotient(cx1, cx, n);
    jet_output(wquot, n, KNRM, KGRY);
    derivative_output(wquot, n, KBLD, KGRY);
    printf("%s\n", KNRM);

    printf("%s%s%s\n", KCYN, "f(x) = x^a, x = 9, a = -0.5", KNRM);
    wpwr = ad_power(cx9, - 0.5q, n);
    jet_output(wpwr, n, KNRM, KGRY);
    derivative_output(wpwr, n, KBLD, KGRY);
    printf("%s\n", KNRM);

    printf("%s%s%s\n", KCYN, "f(x) = e^x", KNRM);
    set_ad_status(cx1, VARIABLE);
    we = ad_exp(cx1, n);
    set_ad_status(cx1, CONSTANT);
    jet_output(we, n, KNRM, KGRY);
    derivative_output(we, n, KBLD, KGRY);
    printf("%s\n", KNRM);

    printf("%s%s%s\n", KCYN, "f(x) = sin(x), f(x) = cos(x), x = 0", KNRM);
    ts_pair = ad_sin_cos(cx0, n, TRIG);
    jet_output(ts_pair.x, n, KNRM, KGRY);
    derivative_output(ts_pair.x, n, KBLD, KGRY);
    printf("%s", KNRM);
    jet_output(ts_pair.y, n, KNRM, KGRY);
    derivative_output(ts_pair.y, n, KBLD, KGRY);
    printf("%s\n", KNRM);

    printf("%s%s%s\n", KCYN, "f(x) = sin(x), f(x) = cos(x), x = pi / 3", KNRM);
    ts_pair = ad_sin_cos(PI_3, n, TRIG);
    jet_output(ts_pair.x, n, KNRM, KGRY);
    derivative_output(ts_pair.x, n, KBLD, KGRY);
    printf("%s", KNRM);
    jet_output(ts_pair.y, n, KNRM, KGRY);
    derivative_output(ts_pair.y, n, KBLD, KGRY);
    printf("%s\n", KNRM);

    printf("%s%s%s\n", KCYN, "f(x) = sinh(x), f(x) = cosh(x), x = pi / 3", KNRM);
    ts_pair = ad_sin_cos(PI_3, n, HYP);
    jet_output(ts_pair.x, n, KNRM, KGRY);
    derivative_output(ts_pair.x, n, KBLD, KGRY);
    printf("%s", KNRM);
    jet_output(ts_pair.y, n, KNRM, KGRY);
    derivative_output(ts_pair.y, n, KBLD, KGRY);
    printf("%s\n", KNRM);

    printf("%s%s%s\n", KCYN, "f(x) = tan(x), f(x) = sec(x)^2, x = pi / 4", KNRM);
    ts_pair = ad_tan_sec2(PI_4, n, TRIG);
    jet_output(ts_pair.x, n, KNRM, KGRY);
    derivative_output(ts_pair.x, n, KBLD, KGRY);
    printf("%s", KNRM);
    jet_output(ts_pair.y, n, KNRM, KGRY);
    derivative_output(ts_pair.y, n, KBLD, KGRY);
    printf("%s\n", KNRM);

    printf("%s%s%s\n", KCYN, "f(x) = tanh(x), f(x) = sech(x)^2, x = pi / 4", KNRM);
    ts_pair = ad_tan_sec2(PI_4, n, HYP);
    jet_output(ts_pair.x, n, KNRM, KGRY);
    derivative_output(ts_pair.x, n, KBLD, KGRY);
    printf("%s", KNRM);
    jet_output(ts_pair.y, n, KNRM, KGRY);
    derivative_output(ts_pair.y, n, KBLD, KGRY);
    printf("%s\n", KNRM);

    printf("%s%s%s\n", KCYN, "f(x, y) = x * y, d/dx", KNRM);
    wprod = ad_product(cx, cy, n);
    jet_output(wprod, n, KNRM, KGRY);
    derivative_output(wprod, n, KBLD, KGRY);
    printf("%s\n", KNRM);

    printf("%s%s%s\n", KCYN, "f(x, y) = e^x / (x - y * sin(x^2), d/dx", KNRM);
    wsq = ad_square(cx, n);
    ts_pair = ad_sin_cos(wsq, n, TRIG);
    wprod = ad_product(cy, ts_pair.x, n);
    wsum = ad_minus(cx, wprod, n);
    we = ad_exp(cx, n);
    wquot = ad_quotient(we, wsum, n);
    jet_output(wquot, n, KNRM, KGRY);
    derivative_output(wquot, n, KBLD, KGRY);
    printf("%s\n", KNRM);

    set_ad_status(cx, CONSTANT);
    set_ad_status(cy, VARIABLE);

    printf("%s%s%s\n", KCYN, "f(x, y) = x * y, d/dy", KNRM);
    wprod = ad_product(cx, cy, n);
    jet_output(wprod, n, KNRM, KGRY);
    derivative_output(wprod, n, KBLD, KGRY);
    printf("%s\n", KNRM);

    printf("%s%s%s\n", KCYN, "f(x, y) = e^x / (x - y * sin(x^2), d/dy", KNRM);
    wsq = ad_square(cx, n);
    ts_pair = ad_sin_cos(wsq, n, TRIG);
    wprod = ad_product(cy, ts_pair.x, n);
    wsum = ad_minus(cx, wprod, n);
    we = ad_exp(cx, n);
    wquot = ad_quotient(we, wsum, n);
    jet_output(wquot, n, KNRM, KGRY);
    derivative_output(wquot, n, KBLD, KGRY);
    printf("%s\n", KNRM);

    return 0;
}

int t_test(int argc, char **argv) {
    assert(argc == 12);
    // initialize from command arguments
    t_stepper(argc, argv, &order, &t, &h, &nsteps);
    x = strtoflt128(argv[5], NULL);
    y = strtoflt128(argv[6], NULL);
    z = strtoflt128(argv[7], NULL);
    s = strtoflt128(argv[8], NULL);
    r = strtoflt128(argv[9], NULL);
    b = strtoflt128(argv[10], NULL) / strtoflt128(argv[11], NULL);
    // initialize the derivative and temporary jets
    cx = t_jet(order + 1);
    cy = t_jet(order + 1);
    cz = t_jet(order + 1);
    // main loop
    for (long step = 1; step < nsteps + 1; step++) {
        // print a line of output
        t_line_output(t, 3, x, y, z);
        // compute the taylor coefficients
        cx[0] = x;
        cy[0] = y;
        cz[0] = z;
        for (int k = 0; k < order; k++) {
            cx[k + 1] = s * (cy[k] - cx[k]) / (k + 1);
            cy[k + 1] = (r * cx[k] - t_product(cx, cz, k) - cy[k]) / (k + 1);
            cz[k + 1] = (t_product(cx, cy, k) - b * cz[k]) / (k + 1);
        }
        // sum the series using Horner's method and advance one step
        x = t_horner(cx, order, h);
        y = t_horner(cy, order, h);
        z = t_horner(cz, order, h);
        t = h * step;
    }
    return 0;
}

int main(int argc, char **argv) {
    if (argc == 4) {
        // ./tsm-lorenz-128 7 2 1
        ad_test(argc, argv);
    } else if (argc == 12) {
        // ./tsm-lorenz-128 16 10 .01 10001 -15.8 -17.48 35.64 10 28 8 3
        t_test(argc, argv);
    }
}

