/*
 *  ODE OpenGL display
 *
 * (c) 2018-2023 m4r35n357@gmail.com (Ian Smith), for licencing see the LICENCE file
 */

#include <stdio.h>
#include <assert.h>
#include <GL/freeglut.h>
#include "taylor-ode.h"
#include "opengl.h"

static void *m;  // the model
static series3 *jets;

point get_current_point (void *data) {
    series3 *j = (series3 *)data;
    return (point){(float)j->x[0], (float)j->y[0], (float)j->z[0]};
}

void Animate () {
    SetupView();

    t->colour = get_colour(colour_index);

    if (mode == BOTH || mode == TRAIL) line_trail(t);

    point p = t->points[newest];
    if (mode == BOTH || mode == POSITION) line_position(p, t->colour, 1.0F);

    if (osd_active) {
        glColor3f(0.0F, 0.5F, 0.5F);
        sprintf(hud, "t: %.1Lf  x: % .1lf  y: % .1lf  z: % .1lf  ", c->step * c->step_size, p.a, p.b, p.c);
        osd(10, glutGet(GLUT_WINDOW_HEIGHT) - 20, hud);
        osd_summary();
    }

    if (!finished && !paused) {
        if (tsm_gen(c, jets, m)) {
            buffer_point();
            t->points[newest] = get_current_point(jets);
        } else finished = 1;
        if (stepping) paused = 1;
    }

    ReDraw();
}

int main (int argc, char **argv) {
    since = clock();
    c = get_c_tsm(argc, argv);
    m = get_p(argc, argv, c->order);

    jets = initial_values(argv, c->order);

    length = (int)strtol(argv[1], NULL, BASE); assert(length >= 0 && length <= c->steps);
    t = malloc(sizeof (trail));
    t->points = malloc((size_t)length * sizeof (point));
    t->points[newest] = get_current_point(jets);

    ApplicationInit(argc, argv, "ODE Plotter");
    glutMainLoop();     // Start the main loop.  glutMainLoop never returns.
    return 0 ;          // Compiler requires this to be here. (Never reached)
}
