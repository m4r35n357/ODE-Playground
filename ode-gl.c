/*
 *  ODE OpenGL display
 *
 * (c) 2018-2022 m4r35n357@gmail.com (Ian Smith), for licencing see the LICENCE file
 */

#include <stdio.h>
#include <time.h>
#include <math.h>
#include <GL/freeglut.h>    // OpenGL Graphics Utility Library
#include "taylor-ode.h"
#include "opengl.h"

static controls *c;
static particle *ball;
static void *p;
static series3 *jets;

static display d;
static char hud[128];
static clock_t since;
static double elapsed, cpu;

static float light_pos[] = { -100.0F, 100.0F, -100.0F, 0.0F };

static _Bool finished = 0, stopped = 0, stepping = 0, running = 1;

void SpecialKeyFunc (int Key, int x, int y) { (void)x; (void)y;
    switch (Key) {
        case    GLUT_KEY_UP: ball->view_latitude += 1.0F; break;
        case  GLUT_KEY_DOWN: ball->view_latitude -= 1.0F; break;
        case  GLUT_KEY_LEFT: ball->view_longitude += 1.0F; break;
        case GLUT_KEY_RIGHT: ball->view_longitude -= 1.0F; break;
    }
}

void KeyPressFunc (unsigned char Key, int x, int y) { (void)x; (void)y;
    switch (Key) {
        case 'B': case 'b': ball->ball_size /= 1.1F; break;
        case 'G': case 'g': ball->ball_size *= 1.1F; break;
        case 'A': case 'a': ball->view_radius -= 0.1F; break;
        case 'Z': case 'z': ball->view_radius += 0.1F; break;
        case 'R': case 'r': running = !running; stopped = 0; break;
        case 'S': case 's': stepping = !stepping; stopped = 0; break;
        case 'F': case 'f': glutFullScreenToggle(); break;
        case  27: exit(1); // Escape key
    }
}

void Animate (void) {
    SetupView(ball->view_radius, ball->view_latitude, ball->view_longitude, light_pos);

    if (d == BOTH || d == LINES) {
        glBegin(GL_LINE_STRIP);
        for (int k = ball->oldest; k != ball->newest; k = (k + 1) % ball->max_points) {  // read buffers
            glColor3f(ball->colour.a, ball->colour.b, ball->colour.c);
            glVertex3f(ball->track[k].a, ball->track[k].b, ball->track[k].c);
        }
        glEnd();
    }

    if (d == BOTH || d == BALLS) {
        glTranslatef((float)jets->x[0], (float)jets->y[0], (float)jets->z[0]);
        glColor3f(ball->colour.a, ball->colour.b, ball->colour.c);
        glutSolidSphere(ball->ball_size, 10, 10);
    }

    sprintf(hud, "t: %.1Lf  x: % .1Lf  y: % .1Lf  z: % .1Lf  ",
                  c->step * c->step_size, jets->x[0], jets->y[0], jets->z[0]);
    osd(10, glutGet(GLUT_WINDOW_HEIGHT) - 20, 0.0F, 0.5F, 0.5F, hud);

    sprintf(hud, "Elapsed: %.1fs  CPU: %.1fs  %.1f %%",
                  elapsed = finished ? elapsed : ((float)(glutGet(GLUT_ELAPSED_TIME)) / 1000.0F),
                  cpu = finished ? cpu : (double)(clock() - since) / CLOCKS_PER_SEC,
                  (float)(100 * c->step) / (float)c->steps);
    osd(10, 10, 0.0F, 0.5F, 0.5F, hud);

    if (!finished && !stopped) {
        if (tsm_gen(c, jets, p)) {
            if (d == BOTH || d == LINES) {  // write buffers
                buffer_point(ball->max_points, &ball->oldest, &ball->newest, &ball->buffers_full);
                ball->track[ball->newest] = (point){(float)jets->x[0], (float)jets->y[0], (float)jets->z[0]};
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
    c = get_c_tsm(argv);
    p = get_p(argc, argv, c->order);
    since = clock();

    jets = malloc(sizeof (series3));
    jets->x = t_jet(c->order + 1); jets->x[0] = strtold(argv[5], NULL);
    jets->y = t_jet(c->order + 1); jets->y[0] = strtold(argv[6], NULL);
    jets->z = t_jet(c->order + 1); jets->z[0] = strtold(argv[7], NULL);

    ball = malloc(sizeof (particle));
    ball->max_points = c->steps / 2;
    ball->oldest = ball->newest = ball->buffers_full = 0;
    ball->colour = (rgb){0.0F, 0.5F, 0.0F };
    ball->track = calloc((size_t)ball->max_points, sizeof (components));
    ball->track[ball->newest] = (point){(float)jets->x[0], (float)jets->y[0], (float)jets->z[0]};
    ball->ball_size = 0.1F;
    ball->view_radius = 20.0F;
    ball->view_longitude = 0.0F;
    ball->view_latitude = 90.0F;

    ApplicationInit(argc, argv, "ODE Plottter");
    glutMainLoop();     // Start the main loop.  glutMainLoop never returns.
    return(0);          // Compiler requires this to be here. (Never reached)
}
