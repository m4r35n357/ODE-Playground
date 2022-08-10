/*
 *  N-Body Plotter
 *
 * (c) 2018-2022 m4r35n357@gmail.com (Ian Smith), for licencing see the LICENCE file
 */

#include <stdio.h>
#include <time.h>
#include <GL/freeglut.h>
#include "symplectic.h"
#include "h-nbody.h"
#include "opengl.h"

static controls *c;  // integrator controls
static nbody *m;     // the model

static display d;
static char hud[128];
static clock_t since;
static double elapsed, cpu;

static float light_pos[] = { -100.0F, 100.0F, -100.0F, 0.0F };

static GLenum finished = GL_FALSE;
static GLenum stopped = GL_FALSE;
static GLenum running = GL_TRUE;
static GLenum stepping = GL_FALSE;

void SpecialKeyFunc (int Key, int x, int y) { (void)x; (void)y;
    switch (Key) {
        case GLUT_KEY_UP: m->view_latitude += 1.0F; break;
        case GLUT_KEY_DOWN: m->view_latitude -= 1.0F; break;
        case GLUT_KEY_LEFT: m->view_longitude += 1.0F; break;
        case GLUT_KEY_RIGHT: m->view_longitude -= 1.0F; break;
    }
}

void KeyPressFunc (unsigned char Key, int x, int y) { (void)x; (void)y;
    switch (Key) {
        case 'B': case 'b': m->ball_scale /= 1.1F; break;
        case 'G': case 'g': m->ball_scale *= 1.1F; break;
        case 'A': case 'a': m->view_radius -= 0.1F; break;
        case 'Z': case 'z': m->view_radius += 0.1F; break;
        case 'R': case 'r': running = !running; stopped = GL_FALSE; break;
        case 'S': case 's': stepping = !stepping; stopped = GL_FALSE; break;
        case 'F': case 'f': glutFullScreenToggle(); break;
        case 27: exit(1); // Escape key
    }
}

void Animate (void) {
    SetupView(m->view_radius, m->view_latitude, m->view_longitude, light_pos);

    body *b = m->bodies;
    if (d == BOTH || d == LINES) {
        for (int i = 0; i < m->n; i += 1) {
            glBegin(GL_LINE_STRIP);
            for (int k = m->oldest; k != m->newest; k = (k + 1) % m->max_points) {  // read buffers
                glColor3f(b[i].colour.a, b[i].colour.b, b[i].colour.c);
                glVertex3f(b[i].track[k].a, b[i].track[k].b, b[i].track[k].c);
            }
            glEnd();
        }
    }

    if (d == BOTH || d == BALLS) {
        glTranslatef((float)(b[0].x - m->centre.x), (float)(b[0].y - m->centre.y), (float)(b[0].z - m->centre.z));
        glColor3f(b[0].colour.a, b[0].colour.b, b[0].colour.c);
        glutSolidSphere(m->ball_scale * b[0].r, 10, 10);
        for (int i = 1; i < m->n; i += 1) {
            glTranslatef((float)(b[i].x - b[i - 1].x), (float)(b[i].y - b[i - 1].y), (float)(b[i].z - b[i - 1].z));
            glColor3f(b[i].colour.a, b[i].colour.b, b[i].colour.c);
            glutSolidSphere(m->ball_scale * b[i].r, 10, 10);
        }
    }

    sprintf(hud, "t: %.1Lf  h: %.6Le  ~sf: %.1Lf", c->step * c->step_size, m->h, error(m->h - m->h0));
    osd(10, glutGet(GLUT_WINDOW_HEIGHT) - 20, 0.0F, 0.5F, 0.5F, hud);

    sprintf(hud, "Elapsed: %.1fs  CPU: %.1fs  %.1f %%",
                  elapsed = finished ? elapsed : ((float)(glutGet(GLUT_ELAPSED_TIME)) / 1000.0F),
                  cpu = finished ? cpu : (double)(clock() - since) / CLOCKS_PER_SEC,
                  (float)(100 * c->step) / (float)c->steps);
    osd(10, 10, 0.0F, 0.5F, 0.5F, hud);

    if (!finished && !stopped) {
        if (generate(c, m)) {
            cog(m);
            m->h = h(m) > m->h ? h(m) : m->h;
            if (d == BOTH || d == LINES) {  // write buffers
                buffer_point(m->max_points, &m->oldest, &m->newest, &m->buffers_full);
                for (int i = 0; i < m->n; i += 1) {
                    b[i].track[m->newest] = (point){(float)b[i].x, (float)b[i].y, (float)b[i].z};
                }
            }
        } else {
            finished = GL_TRUE;
        }
        if (stepping) {
            stopped = GL_TRUE;
        }
    }

    ReDraw();
}

void CloseWindow (void) {
    fprintf(stderr, "H : % .18Le\n", m->h);
}

int main (int argc, char** argv) {
    d = (display)strtol(argv[1], NULL, BASE);
    c = get_c(argv);
    m = (nbody *)get_p(argc, argv, (argc - 7) / 7);

    fprintf(stderr, "\n");
    fprintf(stderr, "H0: % .18Le\n", m->h);
    since = clock();

    ApplicationInit(argc, argv, "N-Body Plotter");
    glutCloseFunc(CloseWindow);
    glutMainLoop();     // Start the main loop.  glutMainLoop never returns.
    return(0);          // Compiler requires this to be here. (Never reached)
}
