/*
 * N-body problem using Hamilton's equations
 *
 * (c) 2018-2022 m4r35n357@gmail.com (Ian Smith), for licencing see the LICENCE file
 */

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>
#include "symplectic.h"
#include "h-nbody.h"

void *get_p (int argc, char **argv, int n_bodies) {
    fprintf(stderr, "[ "); for (int i = 0; i < argc; i++) fprintf(stderr, "%s ", argv[i]); fprintf(stderr, "]\n");
    assert(argc == 7 + 7 * n_bodies);
    rgb colours[] = {
        (rgb){1.0F, 1.0F, 0.0F}, (rgb){0.0F, 1.0F, 1.0F}, (rgb){1.0F, 0.0F, 1.0F},
        (rgb){1.0F, 0.0F, 0.0F}, (rgb){0.0F, 1.0F, 0.0F}, (rgb){0.0F, 0.0F, 1.0F},
        (rgb){0.2F, 0.2F, 0.2F}, (rgb){0.8F, 0.8F, 0.8F}, (rgb){0.5F, 0.5F, 0.5F}
    };
    nbody *nb = malloc(sizeof (nbody));
    nb->max_points = (int)strtol(argv[5], NULL, BASE);
    nb->oldest = nb->current = nb->buffers_full = 0;
    nb->g = strtold(argv[6], NULL);
    nb->n = n_bodies;
    nb->bodies = calloc((size_t)nb->n, sizeof (body));
    for (int i = 0; i < nb->n; i += 1) {
        nb->bodies[i].m = strtold(argv[7 * i + 7], NULL);
        nb->bodies[i].r = (float)powl(nb->bodies[i].m, 1.0L / 3.0L);
        nb->bodies[i].x = strtold(argv[7 * i + 8], NULL);
        nb->bodies[i].y = strtold(argv[7 * i + 9], NULL);
        nb->bodies[i].z = strtold(argv[7 * i + 10], NULL);
        nb->bodies[i].px = strtold(argv[7 * i + 11], NULL);
        nb->bodies[i].py = strtold(argv[7 * i + 12], NULL);
        nb->bodies[i].pz = strtold(argv[7 * i + 13], NULL);
        nb->bodies[i].colour = colours[i % 9];
    }
    cog(nb);
    for (int i = 0; i < nb->n; i += 1) {
        nb->bodies[i].track = calloc((size_t)nb->max_points, sizeof (components));
        nb->bodies[i].track[nb->current] = (point){(float)nb->bodies[i].x, (float)nb->bodies[i].y, (float)nb->bodies[i].z};
    }
    nb->h = nb->h0 = h(nb);
    nb->ball_scale = 0.1F;
    nb->view_radius = 20.0F;
    nb->view_longitude = 0.0F;
    nb->view_latitude = 90.0F;
    return nb;
}

void cog (nbody *nb) {
    body *b = nb->bodies;
    real X = 0.0L, Y = 0.0L, Z = 0.0L, M = 0.0L;
    for (int i = 0; i < nb->n; i += 1) {
        X += b[i].x * b[i].m;
        Y += b[i].y * b[i].m;
        Z += b[i].z * b[i].m;
        M += b[i].m;
    }
    nb->centre = (components){X / M, Y / M, Z / M};
    for (int i = 0; i < nb->n; i += 1) {
        b[i].x -= nb->centre.x;
        b[i].y -= nb->centre.y;
        b[i].z -= nb->centre.z;
    }
}

static real distance (real xA, real yA, real zA, real xB, real yB, real zB) {
    return sqrtl((xB - xA) * (xB - xA) + (yB - yA) * (yB - yA) + (zB - zA) * (zB - zA));
}

real h (nbody *nb) {
    body *b = nb->bodies;
    real e = 0.0L;
    for (int i = 0; i < nb->n; i += 1) {
        e += 0.5L * (b[i].px * b[i].px + b[i].py * b[i].py + b[i].pz * b[i].pz) / b[i].m;
        for (int j = 0; j < i; j += 1) {
            e -= nb->g * b[i].m * b[j].m / distance(b[i].x, b[i].y, b[i].z, b[j].x, b[j].y, b[j].z);
        }
    }
    return e;
}

void update_q (void *n_body, real c) {
    nbody *nb = (nbody *)n_body;
    body *b = nb->bodies;
    for (int i = 0; i < nb->n; i += 1) {
        real _ = c / b[i].m;
        b[i].x += b[i].px * _;
        b[i].y += b[i].py * _;
        b[i].z += b[i].pz * _;
    }
}

void update_p (void *n_body, real c) {
    nbody *nb = (nbody *)n_body;
    body *b = nb->bodies;
    for (int i = 0; i < nb->n; i += 1) {
        for (int j = 0; j < i; j += 1) {
            real d = distance(b[i].x, b[i].y, b[i].z, b[j].x, b[j].y, b[j].z);
            real _ = - c * nb->g * b[i].m * b[j].m / (d * d * d);
            real dPx = (b[j].x - b[i].x) * _;
            real dPy = (b[j].y - b[i].y) * _;
            real dPz = (b[j].z - b[i].z) * _;
            b[i].px -= dPx;
            b[i].py -= dPy;
            b[i].pz -= dPz;
            b[j].px += dPx;
            b[j].py += dPy;
            b[j].pz += dPz;
        }
    }
}
