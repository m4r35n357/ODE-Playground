/*
 * Compile using c99 -g -Og -o tsm-lorenz-128 taylor-ode-128.c -lquadmath
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

typedef struct { quad x, y; } pair;

typedef enum {TRIG, HYP} geometry;

long order, n;
quad x, y, tmp, *cx, *cy, *cx0, *cx1, *cx9, *PI_3, *PI_4, *we, *wl, *ws, *wc, *wt, *ws2, *wsq, *wsum, *wprod, *wquot, *wpwr;

//long order = 10, nsteps = 1001;
//quad t, x, y, z, s, r, b, h, *cx, *cy, *cz;

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

quad t_horner (quad *jet, long n, quad h) {
    assert(n > 0);
    quad sum = jet[n];
    for (int j = n - 1; j > - 1; j--) {
        sum = sum * h + jet[j];
    }
    return sum;
}

quad t_square (quad *U, int k) {
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

quad t_product (quad *V, quad *U, int k) {
    assert(sizeof *U == sizeof *V);
    assert(k >= 0);
    quad P = 0.0q;
    for (int j = 0; j < k + 1; j++) {
        P += U[j] * V[k - j];
    }
    return P;
}

quad t_quotient (quad *Q, quad *U, quad *V, int k) {
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

quad ddot (quad *V, quad *U, int k) {
    assert(sizeof *U == sizeof *V);
    quad dd = 0.0q;
    for (int j = 1; j < k; j++) {
        dd += j * U[j] * V[k - j];
    }
    return dd / k;
}

quad t_exp (quad *E, quad *U, int k) {
    assert(E != U);
    assert(sizeof *U == sizeof *E);
    assert(k >= 0);
    if (k == 0) {
        return expq(U[0]);
    } else {
        return E[0] * U[k] + ddot(E, U, k);
    }
}

pair t_sin_cos (quad *S, quad *C, quad *U, int k, geometry g) {
    assert(S != C && S != U && C != U);
    assert(sizeof *U == sizeof *S && sizeof *U == sizeof *C);
    assert(k >= 0);
    if (k == 0) {
        if (g == TRIG) {
            sincosq(U[0], S, C);
            return (pair) {S[0], C[0]};
        } else {
            return (pair) {sinhq(U[0]), coshq(U[0])};
        }
    } else {
        quad sin = C[0] * U[k] + ddot(C, U, k);
        quad cos = S[0] * U[k] + ddot(S, U, k);
        if (g == TRIG) {
            return (pair) {sin, - cos};
        } else {
            return (pair) {sin, cos};
        }
    }
}

pair t_tan_sec2 (quad *T, quad *S2, quad *U, int k, geometry g) {
    assert(T != S2 && T != U && S2 != U);
    assert(sizeof *U == sizeof *T && sizeof *U == sizeof *S2);
    assert(k >= 0);
    if (k == 0) {
        if (g == TRIG) {
            quad tan = tanq(U[0]);
            return (pair) {tan, 1.0q + tan * tan};
        } else {
            quad tanh = tanhq(U[0]);
            return (pair) {tanh, 1.0q - tanh * tanh};
        }
    } else {
        quad tan = S2[0] * U[k] + ddot(S2, U, k);
        quad sec2 = 2.0q * (T[0] * tan + ddot(T, T, k));
        if (g == TRIG) {
            return (pair) {tan, sec2};
        } else {
            return (pair) {tan, - sec2};
        }
    }
}

quad t_power (quad *P, quad *U, quad a, int k) {
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

void ad_scale (quad *S, quad *U, quad a, int n) {
    for (int k = 0; k < n; k++) {
        S[k] = a * U[k];
    }
}

void ad_plus (quad *P, quad *U, quad *V, int n) {
    for (int k = 0; k < n; k++) {
        P[k] = U[k] + V[k];
    }
}

void ad_minus (quad *M, quad *U, quad *V, int n) {
    for (int k = 0; k < n; k++) {
        M[k] = U[k] - V[k];
    }
}

void ad_square (quad *S, quad *U, int n) {
    for (int k = 0; k < n; k++) {
        S[k] = t_square(U, k);
    }
}

void ad_product (quad *P, quad *V, quad *U, int n) {
    assert(P != U && P != V);
    for (int k = 0; k < n; k++) {
        P[k] = t_product(V, U, k);
    }
}

void ad_quotient (quad *Q, quad *U, quad *V, int n) {
    for (int k = 0; k < n; k++) {
        Q[k] = t_quotient(Q, U, V, k);
    }
}

void ad_exp (quad *E, quad *U, int n) {
    for (int k = 0; k < n; k++) {
        E[k] = t_exp(E, U, k);
    }
}

void ad_sin_cos (quad *S, quad *C, quad *U, int n) {
    for (int k = 0; k < n; k++) {
        pair result = t_sin_cos(S, C, U, k, TRIG);
        S[k] = result.x;
        C[k] = result.y;
    }
}

void ad_sinh_cosh (quad *S, quad *C, quad *U, int n) {
    for (int k = 0; k < n; k++) {
        pair result = t_sin_cos(S, C, U, k, HYP);
        S[k] = result.x;
        C[k] = result.y;
    }
}

void ad_tan_sec2 (quad *T, quad *S2, quad *U, int n) {
    for (int k = 0; k < n; k++) {
        pair result = t_tan_sec2(T, S2, U, k, TRIG);
        T[k] = result.x;
        S2[k] = result.y;
    }
}

void ad_tanh_sech2 (quad *T, quad *S2, quad *U, int n) {
    for (int k = 0; k < n; k++) {
        pair result = t_tan_sec2(T, S2, U, k, HYP);
        T[k] = result.x;
        S2[k] = result.y;
    }
}

void ad_power (quad *P, quad *U, quad a, int n) {
    for (int k = 0; k < n; k++) {
        P[k] = t_power(P, U, a, k);
    }
}

void jet_output (quad *jet, long n, char* f_colour, char *fk_colour) {
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

int main (int argc, char **argv) {
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
    wquot = t_jet(n);
    wpwr = t_jet(n);
    we = t_jet(n);
    wl = t_jet(n);
    ws = t_jet(n);
    wc = t_jet(n);
    wt = t_jet(n);
    ws2 = t_jet(n);

    printf("\n%sx = %s, y = %s, order = %ld%s\n\n", KBLD, argv[2], argv[3], n - 1, KNRM);

    cx[1] = 1.0q;
    cy[1] = 0.0q;

    printf("%s%s%s\n", KCYN, "f(x) = x", KNRM);
    jet_output(cx, n, KNRM, KGRY);
    derivative_output(cx, n, KBLD, KGRY);
    printf("%s\n", KNRM);

    printf("%s%s%s\n", KCYN, "f(x) = x^2", KNRM);
    ad_product(wprod, cx, cx, n);
    jet_output(wprod, n, KNRM, KGRY);
    ad_square(wsq, cx, n);
    jet_output(wsq, n, KNRM, KGRY);
    derivative_output(wsq, n, KBLD, KGRY);
    printf("%s\n", KNRM);

    printf("%s%s%s\n", KCYN, "f(x) = x^3", KNRM);
    ad_square(wsq, cx, n);
    ad_product(wprod, wsq, cx, n);
    jet_output(wprod, n, KNRM, KGRY);
    derivative_output(wprod, n, KBLD, KGRY);
    printf("%s\n", KNRM);

    printf("%s%s%s\n", KCYN, "f(x) = x^4", KNRM);
    ad_square(wsq, cx, n);
    ad_product(wprod, wsq, wsq, n);
    jet_output(wprod, n, KNRM, KGRY);
    ad_product(wprod, cx, cx, n);
    ad_square(wsq, wprod, n);
    jet_output(wsq, n, KNRM, KGRY);
    derivative_output(wsq, n, KBLD, KGRY);
    printf("%s\n", KNRM);

    printf("%s%s%s\n", KCYN, "f(x) = 1 / x", KNRM);
    ad_quotient(wquot, cx1, cx, n);
    jet_output(wquot, n, KNRM, KGRY);
    derivative_output(wquot, n, KBLD, KGRY);
    printf("%s\n", KNRM);

    printf("%s%s%s\n", KCYN, "f(x) = x^a, x = 9, a = -0.5", KNRM);
    ad_power(wpwr, cx9, - 0.5q, n);
    jet_output(wpwr, n, KNRM, KGRY);
    derivative_output(wpwr, n, KBLD, KGRY);
    printf("%s\n", KNRM);

    printf("%s%s%s\n", KCYN, "f(x) = e^x", KNRM);
    cx1[1] = 1.0q;
    ad_exp(we, cx1, n);
    cx1[1] = 0.0q;
    jet_output(we, n, KNRM, KGRY);
    derivative_output(we, n, KBLD, KGRY);
    printf("%s\n", KNRM);

    printf("%s%s%s\n", KCYN, "f(x) = sin(x), f(x) = cos(x), x = 0", KNRM);
    ad_sin_cos(ws, wc, cx0, n);
    jet_output(ws, n, KNRM, KGRY);
    derivative_output(ws, n, KBLD, KGRY);
    printf("%s", KNRM);
    jet_output(wc, n, KNRM, KGRY);
    derivative_output(wc, n, KBLD, KGRY);
    printf("%s\n", KNRM);

    printf("%s%s%s\n", KCYN, "f(x) = sin(x), f(x) = cos(x), x = pi / 3", KNRM);
    ad_sin_cos(ws, wc, PI_3, n);
    jet_output(ws, n, KNRM, KGRY);
    derivative_output(ws, n, KBLD, KGRY);
    printf("%s", KNRM);
    jet_output(wc, n, KNRM, KGRY);
    derivative_output(wc, n, KBLD, KGRY);
    printf("%s\n", KNRM);

    printf("%s%s%s\n", KCYN, "f(x) = sinh(x), f(x) = cosh(x), x = pi / 3", KNRM);
    ad_sinh_cosh(ws, wc, PI_3, n);
    jet_output(ws, n, KNRM, KGRY);
    derivative_output(ws, n, KBLD, KGRY);
    printf("%s", KNRM);
    jet_output(wc, n, KNRM, KGRY);
    derivative_output(wc, n, KBLD, KGRY);
    printf("%s\n", KNRM);

    printf("%s%s%s\n", KCYN, "f(x) = tan(x), f(x) = sec(x)^2, x = pi / 4", KNRM);
    ad_tan_sec2(wt, ws2, PI_4, n);
    jet_output(wt, n, KNRM, KGRY);
    derivative_output(wt, n, KBLD, KGRY);
    printf("%s", KNRM);
    jet_output(ws2, n, KNRM, KGRY);
    derivative_output(ws2, n, KBLD, KGRY);
    printf("%s\n", KNRM);

    printf("%s%s%s\n", KCYN, "f(x) = tanh(x), f(x) = sech(x)^2, x = pi / 4", KNRM);
    ad_tanh_sech2(wt, ws2, PI_4, n);
    jet_output(wt, n, KNRM, KGRY);
    derivative_output(wt, n, KBLD, KGRY);
    printf("%s", KNRM);
    jet_output(ws2, n, KNRM, KGRY);
    derivative_output(ws2, n, KBLD, KGRY);
    printf("%s\n", KNRM);

    printf("%s%s%s\n", KCYN, "f(x, y) = x * y, d/dx", KNRM);
    ad_product(wprod, cx, cy, n);
    jet_output(wprod, n, KNRM, KGRY);
    derivative_output(wprod, n, KBLD, KGRY);
    printf("%s\n", KNRM);

    printf("%s%s%s\n", KCYN, "f(x, y) = e^x / (x - y * sin(x^2), d/dx", KNRM);
    ad_square(wsq, cx, n);
    ad_sin_cos(ws, wc, wsq, n);
    ad_product(wprod, cy, ws, n);
    ad_minus(wsum, cx, wprod, n);
    ad_exp(we, cx, n);
    ad_quotient(wquot, we, wsum, n);
    jet_output(wquot, n, KNRM, KGRY);
    derivative_output(wquot, n, KBLD, KGRY);
    printf("%s\n", KNRM);

    cx[1] = 0.0q;
    cy[1] = 1.0q;

    printf("%s%s%s\n", KCYN, "f(x, y) = x * y, d/dy", KNRM);
    ad_product(wprod, cx, cy, n);
    jet_output(wprod, n, KNRM, KGRY);
    derivative_output(wprod, n, KBLD, KGRY);
    printf("%s\n", KNRM);

    printf("%s%s%s\n", KCYN, "f(x, y) = e^x / (x - y * sin(x^2), d/dy", KNRM);
    ad_square(wsq, cx, n);
    ad_sin_cos(ws, wc, wsq, n);
    ad_product(wprod, cy, ws, n);
    ad_minus(wsum, cx, wprod, n);
    ad_exp(we, cx, n);
    ad_quotient(wquot, we, wsum, n);
    jet_output(wquot, n, KNRM, KGRY);
    derivative_output(wquot, n, KBLD, KGRY);
    printf("%s\n", KNRM);

    return 0;
}

/*
int main(int argc, char **argv) {
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
*/
