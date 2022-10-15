/*
 *  N-Body OpenGL display
 *
 * (c) 2018-2022 m4r35n357@gmail.com (Ian Smith), for licencing see the LICENCE file
 */

#include <stdio.h>
#include <GL/freeglut.h>
#include "symplectic.h"
#include "h-nbody.h"
#include "opengl.h"

static nbody *m;     // the model

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

static point point_from_model (body *b) {
    return (point){(float)b->x, (float)b->y, (float)b->z};
}

void Animate (void) {
    SetupView(m->view_radius, m->view_latitude, m->view_longitude, light_position);

    body *b = m->bodies;
    if (mode == BOTH || mode == LINES) {
        for (int i = 0; i < m->n; i += 1) {
            glBegin(GL_LINE_STRIP);
            for (int k = m->oldest; k != m->newest; k = (k + 1) % m->max_points) {  // read buffers
                glColor3f(b[i].colour.a, b[i].colour.b, b[i].colour.c);
                glVertex3f(b[i].track[k].a, b[i].track[k].b, b[i].track[k].c);
            }
            glEnd();
        }
    }

    if (mode == BOTH || mode == BALLS) {
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
            buffer_point(m->max_points, &m->oldest, &m->newest, &m->buffers_full);
            for (int i = 0; i < m->n; i += 1) {
                b[i].track[m->newest] = point_from_model(&b[i]);
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

    rgb colours[] = {
        (rgb){1.0F, 1.0F, 0.0F}, (rgb){0.0F, 1.0F, 1.0F}, (rgb){1.0F, 0.0F, 1.0F},
        (rgb){1.0F, 0.0F, 0.0F}, (rgb){0.0F, 1.0F, 0.0F}, (rgb){0.0F, 0.0F, 1.0F},
        (rgb){0.2F, 0.2F, 0.2F}, (rgb){0.8F, 0.8F, 0.8F}, (rgb){0.5F, 0.5F, 0.5F}
    };
    m->max_points = (int)strtol(argv[5], NULL, BASE);
    m->oldest = m->newest = m->buffers_full = 0;
    for (int i = 0; i < m->n; i += 1) {
        m->bodies[i].colour = colours[i % 9];
        m->bodies[i].track = calloc((size_t)m->max_points, sizeof (components));
        m->bodies[i].track[0] = point_from_model(&m->bodies[i]);
    }
    m->ball_scale = 0.1F;
    m->view_radius = 20.0F;
    m->view_longitude = 0.0F;
    m->view_latitude = 90.0F;

    ApplicationInit(argc, argv, "N-Body Plotter");
    glutCloseFunc(CloseWindow);
    glutMainLoop();     // Start the main loop.  glutMainLoop never returns.
    return(0);          // Compiler requires this to be here. (Never reached)
}
