/*
 *  Kerr Metric Plotter
 *
 * Example: ./h-kerr-gl  0 10 .01 10000
 *
 * (c) 2018-2022 m4r35n357@gmail.com (Ian Smith), for licencing see the LICENCE file
 */

#include <stdio.h>
#include <time.h>
#include <GL/freeglut.h>
#include "symplectic.h"
#include "h-kerr.h"
#include "opengl.h"

static controls *c;    // integrator controls
static parameters *m;  // the model

static display d;
static char hud[128];
static clock_t since;
static double elapsed, cpu;

static float light_pos[] = { -100.0F, 100.0F, -100.0F, 0.0F };

static _Bool finished = 0, stopped = 0, stepping = 0, running = 1;

void SpecialKeyFunc (int Key, int x, int y) { (void)x; (void)y;
    switch (Key) {
        case    GLUT_KEY_UP: m->view_latitude += 1.0F; break;
        case  GLUT_KEY_DOWN: m->view_latitude -= 1.0F; break;
        case  GLUT_KEY_LEFT: m->view_longitude += 1.0F; break;
        case GLUT_KEY_RIGHT: m->view_longitude -= 1.0F; break;
    }
}

void KeyPressFunc (unsigned char Key, int x, int y) { (void)x; (void)y;
    switch (Key) {
        case 'B': case 'b': m->ball_scale /= 1.1F; break;
        case 'G': case 'g': m->ball_scale *= 1.1F; break;
        case 'A': case 'a': m->view_radius -= 0.1F; break;
        case 'Z': case 'z': m->view_radius += 0.1F; break;
        case 'R': case 'r': running = !running; stopped = 0; break;
        case 'S': case 's': stepping = !stepping; stopped = 0; break;
        case 'F': case 'f': glutFullScreenToggle(); break;
        case  27: exit(1); // Escape key
    }
}

void Animate (void) {
    SetupView(m->view_radius, m->view_latitude, m->view_longitude, light_pos);

    glColor3f(0.0F, 0.0F, 0.5F);
    glutWireSphere(m->horizon, 20, 20);

    if (d == BOTH || d == LINES) {
        glBegin(GL_LINE_STRIP);
        for (int k = m->oldest; k != m->newest; k = (k + 1) % m->max_points) {  // read buffers
            glColor3f(m->colour.a, m->colour.b, m->colour.c);
            glVertex3f(m->track[k].a, m->track[k].b, m->track[k].c);
        }
        glEnd();
    }

    point p = m->track[m->newest];
    if (d == BOTH || d == BALLS) {
        glTranslatef(p.a, p.b, p.c);
        glColor3f(m->colour.a, m->colour.b, m->colour.c);
        glutSolidSphere(m->ball_scale, 10, 10);
    }

    int window_height = glutGet(GLUT_WINDOW_HEIGHT);
    sprintf(hud, "t: %.1Lf  r:% 5.1Lf  theta:% 6.1Lf  phi:% 6.1Lf  ", c->step * c->step_size, r(m), theta(m), phi(m));
    osd(10, window_height - 20, 0.0F, 0.5F, 0.5F, hud);

    pair speed = gamma(m);
    sprintf(hud, "gamma: %.1Lf  v:% .6Lf", speed.a, speed.b);
    osd(10, window_height - 40, 0.0F, 0.5F, 0.5F, hud);

    sprintf(hud, "Elapsed: %.1fs  CPU: %.1fs  %.1f %%",
                  elapsed = finished ? elapsed : ((float)(glutGet(GLUT_ELAPSED_TIME)) / 1000.0F),
                  cpu = finished ? cpu : (double)(clock() - since) / CLOCKS_PER_SEC,
                  (float)(100 * c->step) / (float)c->steps);
    osd(10, 10, 0.0F, 0.5F, 0.5F, hud);

    if (!finished && !stopped) {
        if (generate(c, m)) {
            if (d == BOTH || d == LINES) {  // write buffers
                buffer_point(m->max_points, &m->oldest, &m->newest, &m->buffers_full);
                m->track[m->newest] = to_xyz(m);
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

int main (int argc, char** argv) {
    d = (display)strtol(argv[1], NULL, BASE);
    c = get_c(argv);
    m = get_p(argc, argv, 5);
    since = clock();

    m->max_points = (int)strtol(argv[5], NULL, BASE);
    m->oldest = m->newest = m->buffers_full = 0;
    m->colour = (rgb){0.0F, 0.5F, 0.0F};
    m->track = calloc((size_t)m->max_points, sizeof (components));
    m->track[m->newest] = to_xyz(m);
    m->ball_scale = 0.1F;
    m->view_radius = 20.0F;
    m->view_longitude = 0.0F;
    m->view_latitude = 90.0F;

    ApplicationInit(argc, argv, "Black Hole Orbit Plotter");
    glutMainLoop();     // Start the main loop.  glutMainLoop never returns.
    return(0);          // Compiler requires this to be here. (Never reached)
}
