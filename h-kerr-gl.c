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

static kerr *m;  // the model

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

static point point_from_model (kerr *p) {
    real ra_sth = sqrtl(p->ra2.val) * sinl(p->q_theta);
    return (point){(float)(ra_sth * cosl(p->q_phi)), (float)(ra_sth * sinl(p->q_phi)), (float)(p->q_r * cosl(p->q_theta))};
}

void Animate (void) {
    SetupView(m->view_radius, m->view_latitude, m->view_longitude, light_position);

    glColor3f(0.0F, 0.0F, 0.5F);
    glutWireSphere(m->horizon, 20, 20);
    glColor3f(m->colour.a, m->colour.b, m->colour.c);

    if (mode == BOTH || mode == LINES) {
        glBegin(GL_LINE_STRIP);
        for (int k = m->oldest; k != m->newest; k = (k + 1) % m->max_points) {  // read buffers
            glVertex3f(m->track[k].a, m->track[k].b, m->track[k].c);
        }
        glEnd();
    }

    point p = m->track[m->newest];
    if (mode == BOTH || mode == BALLS) {
        glTranslatef(p.a, p.b, p.c);
        glutSolidSphere(m->ball_scale, 10, 10);
    }

    int window_height = glutGet(GLUT_WINDOW_HEIGHT);
    real S = sigma(m);
    m->tau += m->step * S;
    sprintf(hud, "tau: %.1Lf  t: %.1Lf  r:% 5.1Lf  theta:% 6.1Lf  phi:% 6.1Lf  ",
                  m->tau, m->q_t, m->q_r, m->q_theta * RAD_TO_DEG - 90.0L, fmodl(m->q_phi * RAD_TO_DEG + 180.0L, 360.0L));
    osd(10, window_height - 20, 0.0F, 0.5F, 0.5F, hud);

    pair speed = gamma_v(m, S);
    sprintf(hud, "gamma: %.1Lf  v:% .6Lf", speed.a, speed.b);
    osd(10, window_height - 40, 0.0F, 0.5F, 0.5F, hud);

    sprintf(hud, "Elapsed: %.1fs  CPU: %.1fs  %.1f %%",
                  elapsed = finished ? elapsed : ((float)(glutGet(GLUT_ELAPSED_TIME)) / 1000.0F),
                  cpu = finished ? cpu : (double)(clock() - since) / CLOCKS_PER_SEC,
                  (float)(100 * c->step) / (float)c->steps);
    osd(10, 10, 0.0F, 0.5F, 0.5F, hud);

    if (!finished && !stopped) {
        if (generate(c, m)) {
            buffer_point(m->max_points, &m->oldest, &m->newest, &m->buffers_full);
            m->track[m->newest] = point_from_model(m);
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
    m = get_p_kerr(argc, argv);
    since = clock();

    m->max_points = (int)strtol(argv[5], NULL, BASE);
    m->oldest = m->newest = m->buffers_full = 0;
    m->colour = (rgb){0.0F, 0.5F, 0.0F};
    m->track = calloc((size_t)m->max_points, sizeof (components));
    m->track[m->newest] = point_from_model(m);
    m->ball_scale = 0.1F;
    m->view_radius = 20.0F;
    m->view_longitude = 0.0F;
    m->view_latitude = 90.0F;

    ApplicationInit(argc, argv, "Black Hole Orbit Plotter");
    glutMainLoop();     // Start the main loop.  glutMainLoop never returns.
    return(0);          // Compiler requires this to be here. (Never reached)
}
