/*
 * (c) 2018-2023 m4r35n357@gmail.com (Ian Smith), for licencing see the LICENCE file
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "symplectic.h"
#include "h-nbody.h"

model *get_p_nbody (int argc, char **argv) {
    CHECK((argc - 6) % 7 == 0);
    model *_ = malloc(sizeof (model)); CHECK(_);
    _->G = strtold(argv[5], NULL); CHECK(_->G > 0.0L);
    _->n = (argc - 6) / 7;
    _->bodies = malloc((size_t)_->n * sizeof (body)); CHECK(_->bodies);
    for (int i = 0; i < _->n; i++) {
        _->bodies[i].m = strtold(argv[7 * i + 6], NULL); CHECK(_->bodies[i].m > 0.0L);
        _->bodies[i].r = (float)powl(_->bodies[i].m, 1.0L / 3.0L);
        _->bodies[i].x = strtold(argv[7 * i + 7], NULL);
        _->bodies[i].y = strtold(argv[7 * i + 8], NULL);
        _->bodies[i].z = strtold(argv[7 * i + 9], NULL);
        _->bodies[i].px = strtold(argv[7 * i + 10], NULL);
        _->bodies[i].py = strtold(argv[7 * i + 11], NULL);
        _->bodies[i].pz = strtold(argv[7 * i + 12], NULL);
    }
    reset_cog(_);
    _->h0 = H(_);
    return _;
}

void reset_cog (model *p) {
    body *b = p->bodies;
    real X = 0.0L, Y = 0.0L, Z = 0.0L, M = 0.0L;
    for (int i = 0; i < p->n; i++) {
        X += b[i].x * b[i].m;
        Y += b[i].y * b[i].m;
        Z += b[i].z * b[i].m;
        M += b[i].m;
    }
    for (int i = 0; i < p->n; i++) {
        b[i].x -= X / M;
        b[i].y -= Y / M;
        b[i].z -= Z / M;
    }
}

static real distance (real x, real y, real z, real X, real Y, real Z) {
    return sqrtl(SQR(x - X) + SQR(y - Y) + SQR(z - Z));
}

real H (model *p) {
    body *b = p->bodies;
    real e = 0.0L;
    for (int i = 0; i < p->n; i++) {
        e += 0.5L * (SQR(b[i].px) + SQR(b[i].py) + SQR(b[i].pz)) / b[i].m;
        for (int j = 0; j < i; j++) {
            e -= p->G * b[i].m * b[j].m / distance(b[i].x, b[i].y, b[i].z, b[j].x, b[j].y, b[j].z);
        }
    }
    return e;
}

void update_q (model *p, real c) {
    body *b = p->bodies;
    for (int i = 0; i < p->n; i++) {
        real _ = c / b[i].m;
        b[i].x += b[i].px * _;
        b[i].y += b[i].py * _;
        b[i].z += b[i].pz * _;
    }
}

void update_p (model *p, real c) {
    body *b = p->bodies;
    for (int i = 0; i < p->n; i++) {
        for (int j = 0; j < i; j++) {
            real d = distance(b[i].x, b[i].y, b[i].z, b[j].x, b[j].y, b[j].z);
            real _ = c * p->G * b[i].m * b[j].m / (d * d * d);
            real dPx = (b[j].x - b[i].x) * _;
            real dPy = (b[j].y - b[i].y) * _;
            real dPz = (b[j].z - b[i].z) * _;
            b[j].px -= dPx; b[i].px += dPx;
            b[j].py -= dPy; b[i].py += dPy;
            b[j].pz -= dPz; b[i].pz += dPz;
        }
    }
}
