/*
 *  Kerr Metric OpenGL display
 *
 * (c) 2018-2022 m4r35n357@gmail.com (Ian Smith), for licencing see the LICENCE file
 */

#include <stdio.h>
#include <math.h>
#include <GL/freeglut.h>
#include "symplectic.h"
#include "h-kerr.h"
#include "opengl.h"

static kerr *k;  // the model

void SpecialKeyFunc (int Key, int x, int y) { (void)x; (void)y;
    switch (Key) {
        case    GLUT_KEY_UP: k->view_latitude += 1.0F; break;
        case  GLUT_KEY_DOWN: k->view_latitude -= 1.0F; break;
        case  GLUT_KEY_LEFT: k->view_longitude += 1.0F; break;
        case GLUT_KEY_RIGHT: k->view_longitude -= 1.0F; break;
    }
}

void KeyPressFunc (unsigned char Key, int x, int y) { (void)x; (void)y;
    switch (Key) {
        case 'B': case 'b': k->ball_scale /= 1.1F; break;
        case 'G': case 'g': k->ball_scale *= 1.1F; break;
        case 'A': case 'a': k->view_radius -= 0.1F; break;
        case 'Z': case 'z': k->view_radius += 0.1F; break;
        case 'R': case 'r': running = !running; stopped = 0; break;
        case 'S': case 's': stepping = !stepping; stopped = 0; break;
        case 'F': case 'f': glutFullScreenToggle(); break;
        case  27: exit(1); // Escape key
    }
}

static point point_from_model (kerr *k) {
    real ra_sth = sqrtl(k->ra2.val) * sinl(k->q_theta);
    return (point){(float)(ra_sth * cosl(k->q_phi)), (float)(ra_sth * sinl(k->q_phi)), (float)(k->q_r * cosl(k->q_theta))};
}

void Animate (void) {
    SetupView(k->view_radius, k->view_latitude, k->view_longitude, light_position);

    glColor3f(0.0F, 0.0F, 0.5F);
    glutWireSphere(k->horizon, 20, 20);
    glColor3f(k->colour.a, k->colour.b, k->colour.c);

    if (mode == BOTH || mode == LINES) {
        glBegin(GL_LINE_STRIP);
        for (int i = k->oldest; i != k->newest; i = (i + 1) % k->max_points) {  // read buffers
            glVertex3f(k->track[i].a, k->track[i].b, k->track[i].c);
        }
        glEnd();
    }

    point p = k->track[k->newest];
    if (mode == BOTH || mode == BALLS) {
        glTranslatef(p.a, p.b, p.c);
        glutSolidSphere(k->ball_scale, 10, 10);
    }

    int window_height = glutGet(GLUT_WINDOW_HEIGHT);
    real S = sigma(k);
    k->tau += k->step * S;
    sprintf(hud, "tau: %.1Lf  t: %.1Lf  r:% 5.1Lf  theta:% 6.1Lf  phi:% 6.1Lf  ",
                  k->tau, k->q_t, k->q_r, k->q_theta * RAD_TO_DEG - 90.0L, fmodl(k->q_phi * RAD_TO_DEG + 180.0L, 360.0L));
    osd(10, window_height - 20, 0.0F, 0.5F, 0.5F, hud);

    pair speed = gamma_v(k, S);
    sprintf(hud, "gamma: %.1Lf  v:% .6Lf", speed.a, speed.b);
    osd(10, window_height - 40, 0.0F, 0.5F, 0.5F, hud);

    sprintf(hud, "Elapsed: %.1fs  CPU: %.1fs  %.1f %%",
                  elapsed = finished ? elapsed : ((float)(glutGet(GLUT_ELAPSED_TIME)) / 1000.0F),
                  cpu = finished ? cpu : (double)(clock() - since) / CLOCKS_PER_SEC,
                  (float)(100 * c->step) / (float)c->steps);
    osd(10, 10, 0.0F, 0.5F, 0.5F, hud);

    if (!finished && !stopped) {
        if (generate(c, k)) {
            buffer_point(k->max_points, &k->oldest, &k->newest, &k->buffers_full);
            k->track[k->newest] = point_from_model(k);
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
    c = get_c_symp(argv);
    k = get_p_kerr(argc, argv);
    since = clock();

    k->max_points = (int)strtol(argv[5], NULL, BASE);
    k->oldest = k->newest = k->buffers_full = 0;
    k->colour = (rgb){0.0F, 0.5F, 0.0F};
    k->track = calloc((size_t)k->max_points, sizeof (components));
    k->track[k->newest] = point_from_model(k);
    k->ball_scale = 0.1F;
    k->view_radius = 20.0F;
    k->view_longitude = 0.0F;
    k->view_latitude = 90.0F;

    ApplicationInit(argc, argv, "Black Hole Orbit Plotter");
    glutMainLoop();     // Start the main loop.  glutMainLoop never returns.
    return(0);          // Compiler requires this to be here. (Never reached)
}
