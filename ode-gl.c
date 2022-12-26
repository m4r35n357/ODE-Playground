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

point point_from_model (void *model) {
    series3 *j = (series3 *)model;
    return (point){(float)j->x[0], (float)j->y[0], (float)j->z[0]};
}

void Animate () {
    SetupView();

    t->colour = get_colour(colour_index);
    glColor3f(t->colour.a, t->colour.b, t->colour.c);

    if (mode == BOTH || mode == TRAIL) {
        glBegin(GL_LINE_STRIP);
        for (int i = oldest; i != newest; i = (i + 1) % length) {  // read buffers
            glVertex3f(t->points[i].a, t->points[i].b, t->points[i].c);
        }
        glVertex3f(t->points[newest].a, t->points[newest].b, t->points[newest].c);
        glEnd();
    }

    point p = t->points[newest];
    if (mode == BOTH || mode == POSITION) {
        glColor3f(0.3F, 0.3F, 0.3F);
        glBegin(GL_LINES);
        glVertex3f(0.0F, 0.0F, 0.0F);
        glVertex3f(p.a, p.b, p.c);
        glEnd();
        glTranslatef(p.a, p.b, p.c);
        glColor3f(t->colour.a, t->colour.b, t->colour.c);
        solid ? glutSolidSphere(ball_scale, mesh, mesh) : glutWireSphere(ball_scale, mesh, mesh);
    }

    if (osd_active) {
        glColor3f(0.0F, 0.5F, 0.5F);

        sprintf(hud, "t: %.1Lf  x: % .1lf  y: % .1lf  z: % .1lf  ", c->step * c->step_size, p.a, p.b, p.c);
        osd(10, glutGet(GLUT_WINDOW_HEIGHT) - 20, hud);

        sprintf(hud, "Elapsed: %.1fs  CPU: %.1fs  %.0f %%",
                      elapsed = finished ? elapsed : 0.001F * (float)glutGet(GLUT_ELAPSED_TIME),
                      cpu = finished ? cpu : (float)(clock() - since) / CLOCKS_PER_SEC,
                      (float)(100.0L * c->step / c->steps));
        osd(10, 10, hud);
    }

    if (!finished && !paused) {
        if (tsm_gen(c, jets, m)) {
            buffer_point();
            t->points[newest] = point_from_model(jets);
        } else {
            finished = 1;
        }
        if (stepping) {
            paused = 1;
        }
    }

    ReDraw();
}

int main (int argc, char** argv) {
    c = get_c_tsm(argc, argv);
    m = get_p(argc, argv, c->order);
    since = clock();

    jets = malloc(sizeof (series3));
    jets->x = t_jet(c->order + 1); jets->x[0] = strtold(argv[5], NULL);
    jets->y = t_jet(c->order + 1); jets->y[0] = strtold(argv[6], NULL);
    jets->z = t_jet(c->order + 1); jets->z[0] = strtold(argv[7], NULL);

    length = (int)strtol(argv[1], NULL, BASE); assert(length >= 0 && length <= c->steps);
    t = malloc(sizeof (trail));
    t->points = calloc((size_t)length, sizeof (point));
    t->points[newest] = point_from_model(jets);

    ApplicationInit(argc, argv, "ODE Plotter");
    glutMainLoop();     // Start the main loop.  glutMainLoop never returns.
    return 0 ;          // Compiler requires this to be here. (Never reached)
}
