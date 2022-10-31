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
            glColor3f(t[j].colour.a, t[j].colour.b, t[j].colour.c);
            glBegin(GL_LINE_STRIP);
            for (int i = oldest; i != newest; i = (i + 1) % length) {  // read buffers
                glVertex3f(t[j].points[i].a, t[j].points[i].b, t[j].points[i].c);
            }
            glVertex3f(t[j].points[newest].a, t[j].points[newest].b, t[j].points[newest].c);
            glEnd();
        }
    }

    if (mode == BOTH || mode == POSITION) {
        point o = (point){(float)m->centre.x, (float)m->centre.y, (float)m->centre.z};
        for (int j = 0; j < m->n; j++) {
            glColor3f(0.3F, 0.3F, 0.3F);
            point p = t[j].points[newest];
            glBegin(GL_LINES);
            glVertex3f(o.a, o.b, o.c);
            glVertex3f(p.a, p.b, p.c);
            glEnd();
        }
        glTranslatef((float)(b[0].x - m->centre.x), (float)(b[0].y - m->centre.y), (float)(b[0].z - m->centre.z));
        glColor3f(t[0].colour.a, t[0].colour.b, t[0].colour.c);
        solid ? glutSolidSphere(ball_scale * b[0].r, mesh, mesh) : glutWireSphere(ball_scale * b[0].r, mesh, mesh);
        for (int j = 1; j < m->n; j++) {
            glTranslatef((float)(b[j].x - b[j - 1].x), (float)(b[j].y - b[j - 1].y), (float)(b[j].z - b[j - 1].z));
            glColor3f(t[j].colour.a, t[j].colour.b, t[j].colour.c);
            solid ? glutSolidSphere(ball_scale * b[j].r, mesh, mesh) : glutWireSphere(ball_scale * b[j].r, mesh, mesh);
        }
    }

    if (osd_active) {
        glColor3f(0.0F, 0.5F, 0.5F);

        sprintf(hud, "t: %.1Lf  h: %.6Le  ~sf: %.1Lf", c->step * c->step_size, m->h, error(m->h - m->h0));
        osd(10, glutGet(GLUT_WINDOW_HEIGHT) - 20, hud);

        sprintf(hud, "Elapsed: %.1fs  CPU: %.1fs  %.0f %%",
                      elapsed = finished ? elapsed : 0.001F * (float)glutGet(GLUT_ELAPSED_TIME),
                      cpu = finished ? cpu : (double)(clock() - since) / CLOCKS_PER_SEC,
                      (float)(100.0L * c->step / c->steps));
        osd(10, 10, hud);
    }

    if (!finished && !paused) {
        if (generate(c, m)) {
            cog(m);
            m->h = h(m) > m->h ? h(m) : m->h;
            buffer_point();
            for (int j = 0; j < m->n; j++) {
                t[j].points[newest] = point_from_model(&b[j]);
            }
        } else {
            finished = 1;
        }
        if (stepping) {
            paused = 1;
        }
    }

    ReDraw();
}

void CloseWindow () {
    fprintf(stderr, "H : % .18Le\n", m->h);
}

int main (int argc, char** argv) {
    c = get_c_symp(argv);
    m = get_p_nbody(argc, argv, (argc - 7) / 7);
    since = clock();
    fprintf(stderr, "\nH0: % .18Le\n", m->h);

    length = (int)strtol(argv[5], NULL, BASE);
    t = calloc((size_t)m->n, sizeof (trail));
    for (int j = 0; j < m->n; j++) {
        t[j].colour = get_colour(j);
        t[j].points = calloc((size_t)length, sizeof (point));
        t[j].points[0] = point_from_model(&m->bodies[j]);
    }

    ApplicationInit(argc, argv, "N-Body Plotter");
    glutCloseFunc(CloseWindow);
    glutMainLoop();     // Start the main loop.  glutMainLoop never returns.
    return(0);          // Compiler requires this to be here. (Never reached)
}
