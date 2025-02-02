/*
 *  Kerr Metric OpenGL display
 *
 * (c) 2018-2025 m4r35n357@gmail.com (Ian Smith), for licencing see the LICENCE file
 */
#include <stdio.h>
#include <math.h>
#include <GL/freeglut.h>
#include "symplectic.h"
#include "h-kerr.h"
#include "opengl.h"

static model *k;  // the model

static real RAD_TO_DEG;

point get_current_point (void *data) {
    model *_ = (model *)data;
    real ra_sth = sqrtl(_->ra2.val) * sinl(_->q_th);
    return (point){(float)(ra_sth * cosl(_->q_ph)), (float)(ra_sth * sinl(_->q_ph)), (float)(_->q_r * cosl(_->q_th))};
}

void Animate () {
    SetupView();

    glColor3f(0.0F, 0.0F, 0.5F);
    solid ? glutSolidSphere(k->horizon, 2 * mesh, 2 * mesh) : glutWireSphere(k->horizon, 2 * mesh, 2 * mesh);

    t->colour = get_colour(colour_index);

    if (mode == BOTH || mode == TRAIL) line_trail(t);

    point p = t->points[newest];
    if (mode == BOTH || mode == POSITION) line_position(p, t->colour, 1.0F);

    if (osd_active) {
        int window_height = glutGet(GLUT_WINDOW_HEIGHT);
        real S = sigma(k);
        k->tau += k->step_size * S;
        glColor3f(0.0F, 0.5F, 0.5F);
        sprintf(hud, "tau: %.0Lf  t: %.0Lf  r:% 5.1Lf  theta:% 4.0Lf  phi:% 4.0Lf  ",
                      k->tau, k->q_t, k->q_r, k->q_th * RAD_TO_DEG - 90.0L, fmodl(k->q_ph * RAD_TO_DEG + 180.0L, 360.0L));
        osd(10, window_height - 20, hud);
        pair speed = gamma_v(k, S);
        sprintf(hud, "gamma: %.1Lf  v:% .6Lf", speed.a, speed.b);
        osd(10, window_height - 40, hud);
        osd_summary();
    }

    if (!finished && !paused) {
        if (generate(c, k)) {
            buffer_point();
            t->points[newest] = get_current_point(k);
        } else finished = true;
        if (stepping) paused = true;
    }

    ReDraw();
}

int main (int argc, char **argv) {
    since = clock();
    c = symp_get_c(argc, argv);
    k = kerr_get_p(argc, argv, c->h);
    RAD_TO_DEG = 180.0L / acosl(-1.0L);

    length = (int)strtol(argv[1], NULL, BASE); CHECK(length >= 0 && length <= c->steps);
    t = malloc(sizeof (trail)); CHECK(t);
    t->points = malloc((size_t)length * sizeof (point)); CHECK(t->points);
    t->points[newest] = get_current_point(k);

    ApplicationInit(argc, argv, "Black Hole Orbit Plotter");
    glutMainLoop();     // Start the main loop.  glutMainLoop never returns.
    return 0 ;          // Compiler requires this to be here. (Never reached)
}
