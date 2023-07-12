/*
 *  N-Body OpenGL display
 *
 * (c) 2018-2023 m4r35n357@gmail.com (Ian Smith), for licencing see the LICENCE file
 */

#include <stdio.h>
#include <GL/freeglut.h>
#include "symplectic.h"
#include "opengl.h"
#include "h-nbody.h"

static nbody *nb;  // the model

point get_current_point (void *data) {
    body *b = (body *)data;
    return (point){(float)b->x, (float)b->y, (float)b->z};
}

void Animate () {
    SetupView();

    if (mode == BOTH || mode == TRAIL) {
        for (int j = 0; j < nb->n; j++) {
            line_trail(&t[j]);
        }
    }

    if (mode == BOTH || mode == POSITION) {
        for (int j = 0; j < nb->n; j++) {
            line_position(t[j].points[newest], t[j].colour, nb->bodies[j].r);
        }
    }

    if (osd_active) {
        glColor3f(0.0F, 0.5F, 0.5F);
        real h = hamiltonian(nb);
        sprintf(hud, "t: %.1Lf  h: %.6Le  ~sf: %.1Lf", c->step * c->step_size, h, error(h - nb->h0));
        osd(10, glutGet(GLUT_WINDOW_HEIGHT) - 20, hud);
        osd_summary();
    }

    if (!finished && !paused) {
        if (generate(c, nb)) {
            reset_cog(nb);
            buffer_point();
            for (int j = 0; j < nb->n; j++) {
                t[j].points[newest] = get_current_point(&nb->bodies[j]);
            }
        } else finished = true;
        if (stepping) paused = true;
    }

    ReDraw();
}

void CloseWindow () {
    fprintf(stderr, "H : % .18Le\n", hamiltonian(nb));
}

int main (int argc, char **argv) {
    since = clock();
    c = symp_get_c(argc, argv);
    nb = get_p_nbody(argc, argv);
    fprintf(stderr, "\nH0: % .18Le\n", hamiltonian(nb));

    length = (int)strtol(argv[1], NULL, BASE); CHECK(length >= 0 && length <= c->steps);
    t = malloc((size_t)nb->n * sizeof (trail)); CHECK(t);
    for (int j = 0; j < nb->n; j++) {
        t[j].colour = get_colour(j);
        t[j].points = malloc((size_t)length * sizeof (point)); CHECK(t[j].points);
        t[j].points[0] = get_current_point(&nb->bodies[j]);
    }

    ApplicationInit(argc, argv, "N-Body Plotter");
    glutCloseFunc(CloseWindow);
    glutMainLoop();     // Start the main loop.  glutMainLoop never returns.
    return 0 ;          // Compiler requires this to be here. (Never reached)
}
