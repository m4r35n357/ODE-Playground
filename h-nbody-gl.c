/*
 *  N-Body OpenGL display
 *
 * (c) 2018-2023 m4r35n357@gmail.com (Ian Smith), for licencing see the LICENCE file
 */

#include <stdio.h>
#include <assert.h>
#include <GL/freeglut.h>
#include "symplectic.h"
#include "opengl.h"
#include "h-nbody.h"

static nbody *m;  // the model

point point_from_model (void *model) {
    body *b = (body *)model;
    return (point){(float)b->x, (float)b->y, (float)b->z};
}

void Animate () {
    SetupView();

    body *b = m->bodies;
    if (mode == BOTH || mode == TRAIL) {
        for (int j = 0; j < m->n; j++) {
            line_trail(&t[j]);
        }
    }

    if (mode == BOTH || mode == POSITION) {
        for (int j = 0; j < m->n; j++) {
            line_position(t[j].points[newest], t[j].colour, b[j].r);
        }
    }

    if (osd_active) {
        glColor3f(0.0F, 0.5F, 0.5F);
        real h = hamiltonian(m);
        sprintf(hud, "t: %.1Lf  h: %.6Le  ~sf: %.1Lf", c->step * c->step_size, h, error(h - m->h0));
        osd(10, glutGet(GLUT_WINDOW_HEIGHT) - 20, hud);
        osd_summary();
    }

    if (!finished && !paused) {
        if (generate(c, m)) {
            reset_cog(m);
            buffer_point();
            for (int j = 0; j < m->n; j++) {
                t[j].points[newest] = point_from_model(&b[j]);
            }
        } else finished = 1;
        if (stepping) paused = 1;
    }

    ReDraw();
}

void CloseWindow () {
    fprintf(stderr, "H : % .18Le\n", hamiltonian(m));
}

int main (int argc, char **argv) {
    since = clock();
    c = get_c_symp(argc, argv);
    m = get_p_nbody(argc, argv, (argc - 6) / 7);
    fprintf(stderr, "\nH0: % .18Le\n", hamiltonian(m));

    length = (int)strtol(argv[1], NULL, BASE); assert(length >= 0 && length <= c->steps);
    t = malloc((size_t)m->n * sizeof (trail));
    for (int j = 0; j < m->n; j++) {
        t[j].colour = get_colour(j);
        t[j].points = malloc((size_t)length * sizeof (point));
        t[j].points[0] = point_from_model(&m->bodies[j]);
    }

    ApplicationInit(argc, argv, "N-Body Plotter");
    glutCloseFunc(CloseWindow);
    glutMainLoop();     // Start the main loop.  glutMainLoop never returns.
    return 0 ;          // Compiler requires this to be here. (Never reached)
}
