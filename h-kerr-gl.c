/*
 *  Kerr Metric OpenGL display
 *
 * (c) 2018-2022 m4r35n357@gmail.com (Ian Smith), for licencing see the LICENCE file
 */

#include <stdio.h>
#include <math.h>
#include <GL/freeglut.h>
#include "symplectic.h"
#include "opengl.h"
#include "h-kerr.h"

static kerr *k;  // the model
static track *t;

static real RAD_TO_DEG;

point point_from_model (void *model) {
    kerr *k = (kerr *)model;
    real ra_sth = sqrtl(k->ra2.val) * sinl(k->q_theta);
    return (point){(float)(ra_sth * cosl(k->q_phi)), (float)(ra_sth * sinl(k->q_phi)), (float)(k->q_r * cosl(k->q_theta))};
}

void Animate (void) {
    SetupView(view_radius, view_latitude, view_longitude, light_position);

    point p = t->points[newest];
    glBegin(GL_LINES);
    glColor3f(0.3F, 0.3F, 0.3F);
    glVertex3f(0.0F, 0.0F, 0.0F);
    glVertex3f(p.a, p.b, p.c);
    glEnd();

    glColor3f(0.0F, 0.0F, 0.5F);
    glutWireSphere(k->horizon, 20, 20);

    t->colour = get_colour(colour_index);
    glColor3f(t->colour.a, t->colour.b, t->colour.c);

    if (mode == BOTH || mode == LINES) {
        glBegin(GL_LINE_STRIP);
        for (int i = oldest; i != newest; i = (i + 1) % max_points) {  // read buffers
            glVertex3f(t->points[i].a, t->points[i].b, t->points[i].c);
        }
        glEnd();
    }

    if (mode == BOTH || mode == BALLS) {
        glTranslatef(p.a, p.b, p.c);
        glutSolidSphere(ball_scale, 10, 10);
    }

    int window_height = glutGet(GLUT_WINDOW_HEIGHT);
    real S = sigma(k);
    k->tau += k->step * S;
    glColor3f(0.0F, 0.5F, 0.5F);
    sprintf(hud, "tau: %.0Lf  t: %.0Lf  r:% 5.1Lf  theta:% 4.0Lf  phi:% 4.0Lf  ",
                  k->tau, k->q_t, k->q_r, k->q_theta * RAD_TO_DEG - 90.0L, fmodl(k->q_phi * RAD_TO_DEG + 180.0L, 360.0L));
    osd(10, window_height - 20, hud);

    pair speed = gamma_v(k, S);
    sprintf(hud, "gamma: %.1Lf  v:% .6Lf", speed.a, speed.b);
    osd(10, window_height - 40, hud);

    sprintf(hud, "Elapsed: %.1fs  CPU: %.1fs  %.0f %%",
                  elapsed = finished ? elapsed : ((float)(glutGet(GLUT_ELAPSED_TIME)) / 1000.0F),
                  cpu = finished ? cpu : (double)(clock() - since) / CLOCKS_PER_SEC,
                  (float)(100 * c->step / c->steps));
    osd(10, 10, hud);

    if (!finished && !stopped) {
        if (generate(c, k)) {
            buffer_point(max_points, &oldest, &newest, &buffers_full);
            t->points[newest] = point_from_model(k);
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
    colour_index = DARK_GREEN;
    RAD_TO_DEG = 180.0L / acosl(-1.0L);

    max_points = (int)strtol(argv[5], NULL, BASE);
    oldest = newest = buffers_full = 0;
    t = malloc(sizeof (track));
    t->points = calloc((size_t)max_points, sizeof (point));
    t->points[newest] = point_from_model(k);

    ApplicationInit(argc, argv, "Black Hole Orbit Plotter");
    glutMainLoop();     // Start the main loop.  glutMainLoop never returns.
    return(0);          // Compiler requires this to be here. (Never reached)
}
