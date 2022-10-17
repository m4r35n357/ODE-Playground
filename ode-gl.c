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

static particle *ball;
static void *m;

static int colour_index = 13;

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
        case 'A': case 'a': ball->view_radius -= 0.1F; break;
        case 'Z': case 'z': ball->view_radius += 0.1F; break;
        case 'D': case 'd': colour_index += 1; break;
        case 'C': case 'c': colour_index -= 1; break;
        case 'G': case 'g': ball->ball_size *= 1.1F; break;
        case 'B': case 'b': ball->ball_size /= 1.1F; break;
        case 'R': case 'r': running = !running; stopped = 0; break;
        case 'S': case 's': stepping = !stepping; stopped = 0; break;
        case 'F': case 'f': glutFullScreenToggle(); break;
        case  27: exit(1); // Escape key
    }
}

static point point_from_model (series3 *j) {
    return (point){(float)j->x[0], (float)j->y[0], (float)j->z[0]};
}

void Animate (void) {
    SetupView(ball->view_radius, ball->view_latitude, ball->view_longitude, light_position);

    point p = ball->track[ball->newest];
    glBegin(GL_LINES);
    glColor3f(0.3F, 0.3F, 0.3F);
    glVertex3f(0.0F, 0.0F, 0.0F);
    glVertex3f(p.a, p.b, p.c);
    glEnd();

    ball->colour = get_colour(colour_index);
    glColor3f(ball->colour.a, ball->colour.b, ball->colour.c);

    if (mode == BOTH || mode == LINES) {
        glBegin(GL_LINE_STRIP);
        for (int i = ball->oldest; i != ball->newest; i = (i + 1) % ball->max_points) {  // read buffers
            glVertex3f(ball->track[i].a, ball->track[i].b, ball->track[i].c);
        }
        glEnd();
    }

    if (mode == BOTH || mode == BALLS) {
        glTranslatef(p.a, p.b, p.c);
        glutSolidSphere(ball->ball_size, 10, 10);
    }

    sprintf(hud, "t: %.1Lf  x: % .1lf  y: % .1lf  z: % .1lf  ", c->step * c->step_size, p.a, p.b, p.c);
    osd(10, glutGet(GLUT_WINDOW_HEIGHT) - 20, 0.0F, 0.5F, 0.5F, hud);

    sprintf(hud, "Elapsed: %.1fs  CPU: %.1fs  %.1f %%",
                  elapsed = finished ? elapsed : ((float)(glutGet(GLUT_ELAPSED_TIME)) / 1000.0F),
                  cpu = finished ? cpu : (double)(clock() - since) / CLOCKS_PER_SEC,
                  (float)(100 * c->step) / (float)c->steps);
    osd(10, 10, 0.0F, 0.5F, 0.5F, hud);

    if (!finished && !stopped) {
        if (tsm_gen(c, ball->jets, m)) {
            buffer_point(ball->max_points, &ball->oldest, &ball->newest, &ball->buffers_full);
            ball->track[ball->newest] = point_from_model(ball->jets);
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

    ball = malloc(sizeof (particle));
    ball->jets = malloc(sizeof (series3));
    ball->jets->x = t_jet(c->order + 1); ball->jets->x[0] = strtold(argv[5], NULL);
    ball->jets->y = t_jet(c->order + 1); ball->jets->y[0] = strtold(argv[6], NULL);
    ball->jets->z = t_jet(c->order + 1); ball->jets->z[0] = strtold(argv[7], NULL);
    ball->max_points = c->steps / 2;
    ball->oldest = ball->newest = ball->buffers_full = 0;
    ball->track = calloc((size_t)ball->max_points, sizeof (components));
    ball->track[ball->newest] = point_from_model(ball->jets);
    ball->ball_size = 0.1F;
    ball->view_radius = 20.0F;
    ball->view_longitude = 0.0F;
    ball->view_latitude = 90.0F;

    ApplicationInit(argc, argv, "ODE Plottter");
    glutMainLoop();     // Start the main loop.  glutMainLoop never returns.
    return(0);          // Compiler requires this to be here. (Never reached)
}
