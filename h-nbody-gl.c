/*
 *  N-Body OpenGL display
 *
 * (c) 2018-2022 m4r35n357@gmail.com (Ian Smith), for licencing see the LICENCE file
 */

#include <stdio.h>
#include <GL/freeglut.h>
#include "symplectic.h"
#include "opengl.h"
#include "h-nbody.h"

static nbody *m;     // the model
static track *t;

point point_from_model (void *model) {
    body *b = (body *)model;
    return (point){(float)b->x, (float)b->y, (float)b->z};
}

void Animate (void) {
    SetupView(view_radius, view_latitude, view_longitude, light_position);

    body *b = m->bodies;
    if (mode == BOTH || mode == LINES) {
        for (int i = 0; i < m->n; i += 1) {
            glBegin(GL_LINE_STRIP);
            for (int k = oldest; k != newest; k = (k + 1) % max_points) {  // read buffers
                glColor3f(t[i].colour.a, t[i].colour.b, t[i].colour.c);
                glVertex3f(t[i].points[k].a, t[i].points[k].b, t[i].points[k].c);
            }
            glEnd();
        }
    }

    if (mode == BOTH || mode == BALLS) {
        glTranslatef((float)(b[0].x - m->centre.x), (float)(b[0].y - m->centre.y), (float)(b[0].z - m->centre.z));
        glColor3f(t[0].colour.a, t[0].colour.b, t[0].colour.c);
        glutSolidSphere(ball_scale * b[0].r, 10, 10);
        for (int i = 1; i < m->n; i += 1) {
            glTranslatef((float)(b[i].x - b[i - 1].x), (float)(b[i].y - b[i - 1].y), (float)(b[i].z - b[i - 1].z));
            glColor3f(t[i].colour.a, t[i].colour.b, t[i].colour.c);
            glutSolidSphere(ball_scale * b[i].r, 10, 10);
        }
    }

    glColor3f(0.0F, 0.5F, 0.5F);
    sprintf(hud, "t: %.1Lf  h: %.6Le  ~sf: %.1Lf", c->step * c->step_size, m->h, error(m->h - m->h0));
    osd(10, glutGet(GLUT_WINDOW_HEIGHT) - 20, hud);

    sprintf(hud, "Elapsed: %.1fs  CPU: %.1fs  %.0f %%",
                  elapsed = finished ? elapsed : ((float)(glutGet(GLUT_ELAPSED_TIME)) / 1000.0F),
                  cpu = finished ? cpu : (double)(clock() - since) / CLOCKS_PER_SEC,
                  (float)(100 * c->step / c->steps));
    osd(10, 10, hud);

    if (!finished && !stopped) {
        if (generate(c, m)) {
            cog(m);
            m->h = h(m) > m->h ? h(m) : m->h;
            buffer_point(max_points, &oldest, &newest, &buffers_full);
            for (int i = 0; i < m->n; i += 1) {
                t[i].points[newest] = point_from_model(&b[i]);
            }
        } else {
            finished = 1;
        }
        if (stepping) {
            stopped = 1;
        }
    }

    ReDraw();
}

void CloseWindow (void) {
    fprintf(stderr, "H : % .18Le\n", m->h);
}

int main (int argc, char** argv) {
    mode = (display)strtol(argv[1], NULL, BASE);
    c = get_c_symp(argv);
    m = get_p_nbody(argc, argv, (argc - 7) / 7);
    fprintf(stderr, "\n");
    fprintf(stderr, "H0: % .18Le\n", m->h);
    since = clock();

    max_points = (int)strtol(argv[5], NULL, BASE);
    oldest = newest = buffers_full = 0;
    t = calloc((size_t)m->n, sizeof (track));
    for (int i = 0; i < m->n; i += 1) {
        t[i].colour = get_colour(i);
        t[i].points = calloc((size_t)max_points, sizeof (point));
        t[i].points[0] = point_from_model(&m->bodies[i]);
    }

    ApplicationInit(argc, argv, "N-Body Plotter");
    glutCloseFunc(CloseWindow);
    glutMainLoop();     // Start the main loop.  glutMainLoop never returns.
    return(0);          // Compiler requires this to be here. (Never reached)
}
