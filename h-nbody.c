/*
 * N-body problem using Hamilton's equations
 *
 * Example:  ./h-nbody-gl  6 8 .01 10000
 *
 * (c) 2018-2022 m4r35n357@gmail.com (Ian Smith), for licencing see the LICENCE file
 */

#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <assert.h>
#include <math.h>
#include "symplectic.h"
#include "h-nbody.h"

void *get_p (int argc, char **argv, int n_bodies) {
    assert((argc - 6) % 7 == 0);
    rgb colours[8];
    colours[0] = (rgb) { .r = 1.0F, .g = 1.0F, .b = 0.0F };
    colours[1] = (rgb) { .r = 0.0F, .g = 1.0F, .b = 1.0F };
    colours[2] = (rgb) { .r = 1.0F, .g = 0.0F, .b = 1.0F };
    colours[3] = (rgb) { .r = 1.0F, .g = 0.0F, .b = 0.0F };
    colours[4] = (rgb) { .r = 0.0F, .g = 1.0F, .b = 0.0F };
    colours[5] = (rgb) { .r = 0.0F, .g = 0.0F, .b = 1.0F };
    colours[6] = (rgb) { .r = 0.3F, .g = 0.3F, .b = 0.3F };
    colours[7] = (rgb) { .r = 0.7F, .g = 0.7F, .b = 0.7F };
    body *bodies = calloc((size_t)n_bodies, sizeof (body));
    for (int i = 0; i < n_bodies; i += 1) {
        bodies[i].m = strtold(argv[7 * i + 6], NULL);
        bodies[i].x = strtold(argv[7 * i + 7], NULL);
        bodies[i].y = strtold(argv[7 * i + 8], NULL);
        bodies[i].z = strtold(argv[7 * i + 9], NULL);
        bodies[i].px = strtold(argv[7 * i + 10], NULL);
        bodies[i].py = strtold(argv[7 * i + 11], NULL);
        bodies[i].pz = strtold(argv[7 * i + 12], NULL);
        bodies[i].colour = colours[i];
    }
    nbody *nb = malloc(sizeof (nbody));
    nb->n = n_bodies;
    nb->bodies = bodies;
    nb->g = strtold(argv[5], NULL);
    nb->h0 = h(nb);
    nb->ball_scale = 0.01F;
    nb->view_radius = 20.0F;
    nb->view_longitude = 0.0F;
    nb->view_latitude = 90.0F;
    return nb;
}

void cog (nbody *nb) {
    real X = 0.0L, Y = 0.0L, Z = 0.0L, mT = 0.0L;
    for (int i = 0; i < nb->n; i += 1) {
        body *a = &nb->bodies[i];
        X += a->x * a->m;
        Y += a->y * a->m;
        Z += a->z * a->m;
        mT += a->m;
    }
    nb->centre = (components) { .x = X / mT, .y = Y / mT, .z = Z / mT };
    for (int i = 0; i < nb->n; i += 1) {
        body *a = &nb->bodies[i];
        a->x -= nb->centre.x;
        a->y -= nb->centre.y;
        a->z -= nb->centre.z;
    }
}

static real distance (real xA, real yA, real zA, real xB, real yB, real zB) {
    return sqrtl((xB - xA) * (xB - xA) + (yB - yA) * (yB - yA) + (zB - zA) * (zB - zA));
}

real h (nbody *nb) {
    real energy = 0.0L;
    for (int i = 0; i < nb->n; i += 1) {
        body *a = &nb->bodies[i];
        energy += 0.5L * (a->px * a->px + a->py * a->py + a->pz * a->pz) / a->m;
        for (int j = 0; j < nb->n; j += 1) {
            if (i > j) {
                body *b = &nb->bodies[j];
                energy -= nb->g * a->m * b->m / distance(a->x, a->y, a->z, b->x, b->y, b->z);
            }
        }
    }
    return energy;
}

void update_q (void *n_body, real c) {
    nbody *nb = (nbody *)n_body;
    real tmp;
    for (int i = 0; i < nb->n; i += 1) {
        body *a = &nb->bodies[i];
        tmp = c / a->m;
        a->x += a->px * tmp;
        a->y += a->py * tmp;
        a->z += a->pz * tmp;
    }
}

void update_p (void *n_body, real c) {
    nbody *nb = (nbody *)n_body;
    for (int i = 0; i < nb->n; i += 1) {
        body *a = &nb->bodies[i];
        for (int j = 0; j < nb->n; j += 1) {
            if (i > j) {
                body *b = &nb->bodies[j];
                real d = distance(a->x, a->y, a->z, b->x, b->y, b->z);
                real tmp = - c * nb->g * a->m * b->m / (d * d * d);
                real dPx = (b->x - a->x) * tmp;
                real dPy = (b->y - a->y) * tmp;
                real dPz = (b->z - a->z) * tmp;
                a->px -= dPx;
                a->py -= dPy;
                a->pz -= dPz;
                b->px += dPx;
                b->py += dPy;
                b->pz += dPz;
            }
        }
    }
}
