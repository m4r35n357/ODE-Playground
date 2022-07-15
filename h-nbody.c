/*
 * N-body problem using Hamilton's equations
 *
 * Example:  ./h-nbody-gl  6 8 .01 10000
 *
 * (c) 2018-2022 m4r35n357@gmail.com (Ian Smith), for licencing see the LICENCE file
 */

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>
#include "symplectic.h"
#include "h-nbody.h"

void cog (nbody *nb) {
    real X = 0.0L, Y = 0.0L, Z = 0.0L, mT = 0.0L;
    for (int i = 0; i < nb->n; i += 1) {
        body *a = &nb->bodies[i];
        X += a->q_x * a->m;
        Y += a->q_y * a->m;
        Z += a->q_z * a->m;
        mT += a->m;
    }
    nb->centre = (components) { .x = X / mT, .y = Y / mT, .z = Z / mT };
    for (int i = 0; i < nb->n; i += 1) {
        body *a = &nb->bodies[i];
        a->q_x -= nb->centre.x;
        a->q_y -= nb->centre.y;
        a->q_z -= nb->centre.z;
    }
}

static real distance (real xA, real yA, real zA, real xB, real yB, real zB) {
    return sqrtl((xB - xA) * (xB - xA) + (yB - yA) * (yB - yA) + (zB - zA) * (zB - zA));
}

static real h (nbody *nb) {
    real energy = 0.0L;
    for (int i = 0; i < nb->n; i += 1) {
        body *a = &nb->bodies[i];
        energy += 0.5L * (a->p_x * a->p_x + a->p_y * a->p_y + a->p_z * a->p_z) / a->m;
        for (int j = 0; j < nb->n; j += 1) {
            if (i > j) {
                body *b = &nb->bodies[j];
                energy -= nb->g * a->m * b->m / distance(a->q_x, a->q_y, a->q_z, b->q_x, b->q_y, b->q_z);
            }
        }
    }
    return energy;
}

void *get_p (int argc, char **argv, int va_begin) { (void)argc; (void)argv; (void)va_begin;
    int n_bodies = 8;
    body *bodies = calloc((size_t)n_bodies, sizeof (body));
    bodies[0] = (body){ .m = 100.0L, .colour = (components) { 1.0F, 1.0F, 0.0F },
                        .q_x = 0.0L, .q_y = 0.0L, .q_z = 0.0L, .p_x = 0.0L, .p_y = 0.0L, .p_z = 0.0L };
    bodies[1] = (body){ .m = 2.0L, .colour = (components) { .x = 0.0F, .y = 1.0F, .z = 1.0F },
                        .q_x = 0.0L, .q_y = 4.5L, .q_z = 0.4L, .p_x = -0.2L, .p_y = 0.0L, .p_z = 1.8L };
    bodies[2] = (body){ .m = 3.0L, .colour = (components) { .x = 1.0F, .y = 0.0F, .z = 1.0F },
                        .q_x = -6.0L, .q_y = 0.0L, .q_z = -0.4L, .p_x = 0.0L, .p_y = -2.0L, .p_z = 1.0L };
    bodies[3] = (body){ .m = 5.0L, .colour = (components) { .x = 1.0F, .y = 0.0F, .z = 0.0F },
                        .q_x = 3.0L, .q_y = 0.0L, .q_z = -0.2L, .p_x = 0.0L, .p_y = 5.8L, .p_z = -0.2L };
    bodies[4] = (body){ .m = 4.0L, .colour = (components) { .x = 0.0F, .y = 1.0F, .z = 0.0F },
                        .q_x = 0.0L, .q_y = -4.0L, .q_z = 0.1L, .p_x = -3.6L, .p_y = 0.0L, .p_z = 0.2L};
    bodies[5] = (body){ .m = 3.0L, .colour = (components) { .x = 0.0F, .y = 0.0F, .z = 1.0F },
                        .q_x = -4.0L, .q_y = 0.0L, .q_z = -0.1L, .p_x = 0.0L, .p_y = -0.2L, .p_z = -2.6L };
    bodies[6] = (body){ .m = 3.0L, .colour = (components) { .x = 0.3F, .y = 0.3F, .z = 0.3F },
                        .q_x = 8.0L, .q_y = 0.0L, .q_z = -0.3L, .p_x = 0.0L, .p_y = 2.0L, .p_z = -0.2L };
    bodies[7] = (body){ .m = 4.0L, .colour = (components) { .x = 0.6F, .y = 0.6F, .z = 0.6F },
                        .q_x = 0.0L, .q_y = 4.0L, .q_z = -0.2L, .p_x = -4.8L, .p_y = 0.0L, .p_z = -0.2L };
    nbody *nb = malloc(sizeof (nbody));
    nb->n = n_bodies;
    nb->bodies = bodies;
    nb->g = 0.05L;
    nb->h0 = h(nb);
    nb->radius = 20.0L;
    nb->view_longitude = 0.0L;
    nb->view_latitude = 90.0L;
    return nb;
}

void update_q (void *n_body, real c) {
    nbody *nb = (nbody *)n_body;
    real tmp;
    for (int i = 0; i < nb->n; i += 1) {
        body *a = &nb->bodies[i];
        tmp = c / a->m;
        a->q_x += a->p_x * tmp;
        a->q_y += a->p_y * tmp;
        a->q_z += a->p_z * tmp;
    }
}

void update_p (void *n_body, real c) {
    nbody *nb = (nbody *)n_body;
    for (int i = 0; i < nb->n; i += 1) {
        body *a = &nb->bodies[i];
        for (int j = 0; j < nb->n; j += 1) {
            if (i > j) {
                body *b = &nb->bodies[j];
                real d = distance(a->q_x, a->q_y, a->q_z, b->q_x, b->q_y, b->q_z);
                real tmp = - c * nb->g * a->m * b->m / (d * d * d);
                real dPx = (b->q_x - a->q_x) * tmp;
                real dPy = (b->q_y - a->q_y) * tmp;
                real dPz = (b->q_z - a->q_z) * tmp;
                a->p_x -= dPx;
                a->p_y -= dPy;
                a->p_z -= dPz;
                b->p_x += dPx;
                b->p_y += dPy;
                b->p_z += dPz;
            }
        }
    }
}
