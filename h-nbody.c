/*
 * N-body problem using Hamilton's equations
 *
 * Example:  ./h-3body-dbg  6 8 1 10000  >/tmp/$USER/data
 * 
 * (c) 2018-2022 m4r35n357@gmail.com (Ian Smith), for licencing see the LICENCE file
 */

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>
#include "symplectic.h"
#include "h-nbody.h"

static real distance (real xA, real yA, real zA, real xB, real yB, real zB) {
    return sqrtl((xB - xA) * (xB - xA) + (yB - yA) * (yB - yA) + (zB - zA) * (zB - zA));
}

static real h (nbody *nb) {
    real energy = 0.0L;
    for (int i = 0; i < nb->n; i += 1) {
        body a = nb->bodies[i];
        energy += 0.5L * (a.p_x * a.p_x + a.p_y * a.p_y + a.p_z * a.p_z) / a.m;
        for (int j = 0; j < nb->n; j += 1) {
            if (i > j) {
                body b = nb->bodies[j];
                energy -= nb->g * a.m * b.m / distance(a.q_x, a.q_y, a.q_z, b.q_x, b.q_y, b.q_z);
            }
        }
    }
    return energy;
}

void *get_p (int argc, char **argv, int va_begin) { (void)argc; (void)argv; (void)va_begin;
	int n_bodies = 8;
    nbody *nb = calloc((size_t)n_bodies, sizeof (body));
	nb->n = n_bodies;
    nb->bodies[0] = (body){.m = 1.0L, .q_x = 0.0L, .q_y = 0.0L, .q_z = 0.0L, .p_x = 0.0L, .p_y = 0.0L, .p_z = 0.0L};
    nb->bodies[1] = (body){.m = 1.0L, .q_x = 0.0L, .q_y = 0.0L, .q_z = 0.0L, .p_x = 0.0L, .p_y = 0.0L, .p_z = 0.0L};
    nb->bodies[2] = (body){.m = 1.0L, .q_x = 0.0L, .q_y = 0.0L, .q_z = 0.0L, .p_x = 0.0L, .p_y = 0.0L, .p_z = 0.0L};
    nb->bodies[3] = (body){.m = 1.0L, .q_x = 0.0L, .q_y = 0.0L, .q_z = 0.0L, .p_x = 0.0L, .p_y = 0.0L, .p_z = 0.0L};
    nb->bodies[4] = (body){.m = 1.0L, .q_x = 0.0L, .q_y = 0.0L, .q_z = 0.0L, .p_x = 0.0L, .p_y = 0.0L, .p_z = 0.0L};
    nb->bodies[5] = (body){.m = 1.0L, .q_x = 0.0L, .q_y = 0.0L, .q_z = 0.0L, .p_x = 0.0L, .p_y = 0.0L, .p_z = 0.0L};
    nb->bodies[6] = (body){.m = 1.0L, .q_x = 0.0L, .q_y = 0.0L, .q_z = 0.0L, .p_x = 0.0L, .p_y = 0.0L, .p_z = 0.0L};
    nb->bodies[7] = (body){.m = 1.0L, .q_x = 0.0L, .q_y = 0.0L, .q_z = 0.0L, .p_x = 0.0L, .p_y = 0.0L, .p_z = 0.0L};
	nb->g = 6.674e-11L;
	nb->h0 = h(nb);
    return nb;
}

void update_q (void *n_body, real c) {
    nbody *nb = (nbody *)n_body;
    real tmp;
    for (int i = 0; i < nb->n; i += 1) {
        body a = nb->bodies[i];
        tmp = c / a.m;
        a.q_x += a.p_x * tmp;
        a.q_y += a.p_y * tmp;
        a.q_z += a.p_z * tmp;
    }
}


void update_p (void *n_body, real c) {
    nbody *nb = (nbody *)n_body;
    for (int i = 0; i < nb->n; i += 1) {
        body a = nb->bodies[i];
        for (int j = 0; j < nb->n; j += 1) {
            if (i > j) {
                body b = nb->bodies[j];
                real d = distance(a.q_x, a.q_y, a.q_z, b.q_x, b.q_y, b.q_z);
                real tmp = - c * nb->g * a.m * b.m / (d * d * d);
                real dPx = (b.q_x - a.q_x) * tmp;
                real dPy = (b.q_y - a.q_y) * tmp;
                real dPz = (b.q_z - a.q_z) * tmp;
                a.p_x -= dPx;
                a.p_y -= dPy;
                a.p_z -= dPz;
                b.p_x += dPx;
                b.p_y += dPy;
                b.p_z += dPz;
            }
        }
    }
}
