/*
 *  ODE OpenGL display
 *
 * (c) 2018-2022 m4r35n357@gmail.com (Ian Smith), for licencing see the LICENCE file
 */

#include <stdio.h>
#include <math.h>
#include <GL/freeglut.h>    // OpenGL Graphics Utility Library
#include "taylor-ode.h"
#include "opengl.h"

/*
 * Particle/Body tracks
 */
typedef struct Particle {
    series3 *jets;
    struct triple_f colour, *track;
} particle;

static particle *p;
static void *m;

point point_from_model (void *model) {
    series3 *j = (series3 *)model;
    return (point){(float)j->x[0], (float)j->y[0], (float)j->z[0]};
}

void Animate (void) {
    SetupView(view_radius, view_latitude, view_longitude, light_position);

    point q = p->track[newest];
    glBegin(GL_LINES);
    glColor3f(0.3F, 0.3F, 0.3F);
    glVertex3f(0.0F, 0.0F, 0.0F);
    glVertex3f(q.a, q.b, q.c);
    glEnd();

    p->colour = get_colour(colour_index);
    glColor3f(p->colour.a, p->colour.b, p->colour.c);

    if (mode == BOTH || mode == LINES) {
        glBegin(GL_LINE_STRIP);
        for (int i = oldest; i != newest; i = (i + 1) % max_points) {  // read buffers
            glVertex3f(p->track[i].a, p->track[i].b, p->track[i].c);
        }
        glEnd();
    }

    if (mode == BOTH || mode == BALLS) {
        glTranslatef(q.a, q.b, q.c);
        glutSolidSphere(ball_scale, 10, 10);
    }

    glColor3f(0.0F, 0.5F, 0.5F);
    sprintf(hud, "t: %.1Lf  x: % .1lf  y: % .1lf  z: % .1lf  ", c->step * c->step_size, q.a, q.b, q.c);
    osd(10, glutGet(GLUT_WINDOW_HEIGHT) - 20, hud);

    sprintf(hud, "Elapsed: %.1fs  CPU: %.1fs  %.1f %%",
                  elapsed = finished ? elapsed : ((float)(glutGet(GLUT_ELAPSED_TIME)) / 1000.0F),
                  cpu = finished ? cpu : (double)(clock() - since) / CLOCKS_PER_SEC,
                  (float)(100 * c->step) / (float)c->steps);
    osd(10, 10, hud);

    if (!finished && !stopped) {
        if (tsm_gen(c, p->jets, m)) {
            buffer_point(max_points, &oldest, &newest, &buffers_full);
            p->track[newest] = point_from_model(p->jets);
        } else {
            finished = 1;
        }
        if (stepping) {
            stopped = 1;
        }
    }

    ReDraw();
}

int main (int argc, char** argv) {
    mode = (display)strtol(argv[1], NULL, BASE);
    c = get_c_tsm(argv);
    m = get_p(argc, argv, c->order);
    since = clock();
    colour_index = 13;

    p = malloc(sizeof (particle));
    p->jets = malloc(sizeof (series3));
    p->jets->x = t_jet(c->order + 1); p->jets->x[0] = strtold(argv[5], NULL);
    p->jets->y = t_jet(c->order + 1); p->jets->y[0] = strtold(argv[6], NULL);
    p->jets->z = t_jet(c->order + 1); p->jets->z[0] = strtold(argv[7], NULL);
    max_points = c->steps / 2;
    oldest = newest = buffers_full = 0;
    p->track = calloc((size_t)max_points, sizeof (components));
    p->track[newest] = point_from_model(p->jets);

    ApplicationInit(argc, argv, "ODE Plottter");
    glutMainLoop();     // Start the main loop.  glutMainLoop never returns.
    return(0);          // Compiler requires this to be here. (Never reached)
}
