/*
 * (c) 2018-2023 m4r35n357@gmail.com (Ian Smith), for licencing see the LICENCE file
 */

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>
#include "symplectic.h"
#include "h-nbody.h"

nbody *get_p_nbody (int argc, char **argv, int n_bodies) {
    CHECK(argc == 6 + 7 * n_bodies);
    nbody *nb = malloc(sizeof (nbody));
    nb->g = strtold(argv[5], NULL);
    nb->n = n_bodies;
    nb->bodies = malloc((size_t)nb->n * sizeof (body));
    for (int i = 0; i < nb->n; i++) {
        nb->bodies[i].m = strtold(argv[7 * i + 6], NULL);
        nb->bodies[i].r = (float)powl(nb->bodies[i].m, 1.0L / 3.0L);
        nb->bodies[i].x = strtold(argv[7 * i + 7], NULL);
        nb->bodies[i].y = strtold(argv[7 * i + 8], NULL);
        nb->bodies[i].z = strtold(argv[7 * i + 9], NULL);
        nb->bodies[i].px = strtold(argv[7 * i + 10], NULL);
        nb->bodies[i].py = strtold(argv[7 * i + 11], NULL);
        nb->bodies[i].pz = strtold(argv[7 * i + 12], NULL);
    }
    reset_cog(nb);
    nb->h0 = hamiltonian(nb);
    return nb;
}

void reset_cog (nbody *nb) {
    body *b = nb->bodies;
    real X = 0.0L, Y = 0.0L, Z = 0.0L, M = 0.0L;
    for (int i = 0; i < nb->n; i++) {
        X += b[i].x * b[i].m;
        Y += b[i].y * b[i].m;
        Z += b[i].z * b[i].m;
        M += b[i].m;
    }
    for (int i = 0; i < nb->n; i++) {
        b[i].x -= X / M;
        b[i].y -= Y / M;
        b[i].z -= Z / M;
    }
}

static real distance (real x, real y, real z, real X, real Y, real Z) {
    return sqrtl((x - X) * (x - X) + (y - Y) * (y - Y) + (z - Z) * (z - Z));
}

real hamiltonian (nbody *nb) {
    body *b = nb->bodies;
    real e = 0.0L;
    for (int i = 0; i < nb->n; i++) {
        e += 0.5L * (b[i].px * b[i].px + b[i].py * b[i].py + b[i].pz * b[i].pz) / b[i].m;
        for (int j = 0; j < i; j++) {
            e -= nb->g * b[i].m * b[j].m / distance(b[i].x, b[i].y, b[i].z, b[j].x, b[j].y, b[j].z);
        }
    }
    return e;
}

void update_q (void *n_body, real c) {
    nbody *nb = (nbody *)n_body;
    body *b = nb->bodies;
    for (int i = 0; i < nb->n; i++) {
        real _ = c / b[i].m;
        b[i].x += b[i].px * _;
        b[i].y += b[i].py * _;
        b[i].z += b[i].pz * _;
    }
}

void update_p (void *n_body, real c) {
    nbody *nb = (nbody *)n_body;
    body *b = nb->bodies;
    for (int i = 0; i < nb->n; i++) {
        for (int j = 0; j < i; j++) {
            real d = distance(b[i].x, b[i].y, b[i].z, b[j].x, b[j].y, b[j].z);
            real _ = c * nb->g * b[i].m * b[j].m / (d * d * d);
            real dPx = (b[j].x - b[i].x) * _;
            real dPy = (b[j].y - b[i].y) * _;
            real dPz = (b[j].z - b[i].z) * _;
            b[j].px -= dPx; b[i].px += dPx;
            b[j].py -= dPy; b[i].py += dPy;
            b[j].pz -= dPz; b[i].pz += dPz;
        }
    }
}
