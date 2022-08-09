/*
 *  Kerr Metric simulator
 *
 * Example: ./h-kerr-gl  0 10 .01 10000
 *
 * (c) 2018-2022 m4r35n357@gmail.com (Ian Smith), for licencing see the LICENCE file
 */

#include <stdio.h>
#include <assert.h>
#include <string.h>
#include <time.h>
#include <GL/glew.h>
#include <GL/freeglut.h>    // OpenGL Graphics Utility Library
#include "symplectic.h"
#include "h-kerr.h"
#include "opengl.h"

static controls *c;
static parameters *bh;

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
        case GLUT_KEY_UP: bh->view_latitude += 1.0F; break;
        case GLUT_KEY_DOWN: bh->view_latitude -= 1.0F; break;
        case GLUT_KEY_LEFT: bh->view_longitude += 1.0F; break;
        case GLUT_KEY_RIGHT: bh->view_longitude -= 1.0F; break;
    }
}

void KeyPressFunc (unsigned char Key, int x, int y) { (void)x; (void)y;
    switch (Key) {
        case 'B': case 'b': bh->ball_size /= 1.1F; break;
        case 'G': case 'g': bh->ball_size *= 1.1F; break;
        case 'A': case 'a': bh->view_radius -= 0.1F; break;
        case 'Z': case 'z': bh->view_radius += 0.1F; break;
        case 'R': case 'r': running = !running; stopped = GL_FALSE; break;
        case 'S': case 's': stepping = !stepping; stopped = GL_FALSE; break;
        case 'F': case 'f': glutFullScreenToggle(); break;
        case 27: exit(1); // Escape key
    }
}

void osd (int x, int y, float r, float g, float b, char *string) {
    glColor3f(r, g, b);
    glWindowPos2i(x, y);
    int len = (int)strlen(string);
    for (int i = 0; i < len; i++) {
        glutBitmapCharacter(GLUT_BITMAP_9_BY_15, string[i]);
    }
}

void Animate (void) {
    // Clear the rendering window
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    glLoadIdentity();
    glTranslatef(0.0F, 0.0F, - bh->view_radius);
    glRotatef(bh->view_latitude, 1.0F, 0.0F, 0.0F);
    glRotatef(bh->view_longitude, 0.0F, 0.0F, 1.0F);

    glLightfv(GL_LIGHT0, GL_POSITION, &light_pos[0]);

    glColor3f(0.0F, 0.0F, 0.5F);
    glutWireSphere(bh->horizon, 20, 20);

    if (d == BOTH || d == LINES) {
        glBegin(GL_LINE_STRIP);
        for (int k = bh->oldest; k != bh->newest; k = (k + 1) % bh->max_points) {  // read buffers
            glColor3f(bh->colour.a, bh->colour.b, bh->colour.c);
            glVertex3f(bh->track[k].a, bh->track[k].b, bh->track[k].c);
        }
        glEnd();
    }

    point p = bh->track[bh->newest];
    if (d == BOTH || d == BALLS) {
        glTranslatef(p.a, p.b, p.c);
        glColor3f(bh->colour.a, bh->colour.b, bh->colour.c);
        glutSolidSphere(bh->ball_size, 10, 10);
    }

    int window_height = glutGet(GLUT_WINDOW_HEIGHT);
    sprintf(hud, "t: %.1Lf  r:% 5.1Lf  theta:% 6.1Lf  phi:% 6.1Lf  ", c->step * c->step_size, r(bh), theta(bh), phi(bh));
    osd(10, window_height - 20, 0.0F, 0.5F, 0.5F, hud);

    pair speed = gamma(bh);
    sprintf(hud, "gamma: %.1Lf  v:% .6Lf", speed.a, speed.b);
    osd(10, window_height - 40, 0.0F, 0.5F, 0.5F, hud);

    sprintf(hud, "Elapsed: %.1fs  CPU: %.1fs  %.1f %%",
                  elapsed = finished ? elapsed : ((float)(glutGet(GLUT_ELAPSED_TIME)) / 1000.0F),
                  cpu = finished ? cpu : (double)(clock() - since) / CLOCKS_PER_SEC,
                  (float)(100 * c->step) / (float)c->steps);
    osd(10, 10, 0.0F, 0.5F, 0.5F, hud);

    if (!finished && !stopped) {
        if (generate(c, bh)) {
            if (d == BOTH || d == LINES) {  // write buffers
                bh->newest += 1;
                if (!bh->buffers_full && (bh->newest == bh->max_points)) {
                    bh->buffers_full = 1;
                }
                if (bh->buffers_full) {
                    bh->oldest = (bh->newest + 1) % bh->max_points;
                    bh->newest %= bh->max_points;
                }
                bh->track[bh->newest] = to_xyz(bh);
            }
        } else {
            finished = GL_TRUE;
        }
        if (stepping) {
            stopped = GL_TRUE;
        }
    }

    // Flush the pipeline, and swap the buffers
    glFlush();
    glutSwapBuffers();
    glutPostRedisplay();        // Request a re-draw for animation purposes
}

void OpenGLInit (void) {
    glShadeModel(GL_FLAT);
    glClearColor(0.0F, 0.0F, 0.0F, 0.0F);
    glClearDepth(1.0F);
    glEnable(GL_DEPTH_TEST);
    glEnable(GL_LIGHTING);
    glEnable(GL_LIGHT0);
    glEnable(GL_COLOR_MATERIAL);
}

void ResizeWindow (int w, int h) {
    float aspectRatio;
    h = (h == 0) ? 1 : h;
    w = (w == 0) ? 1 : w;
    glViewport(0, 0, w, h);   // View port uses whole window
    aspectRatio = (float)w / (float)h;
    // Set up the projection view matrix (not very well!)
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    gluPerspective(60.0F, aspectRatio, 1.0F, 100.0F);
    // Select the Modelview matrix
    glMatrixMode(GL_MODELVIEW);
}

int main (int argc, char** argv) {
    d = (display)strtol(argv[1], NULL, BASE);
    c = get_c(argv);
    bh = get_p(argc, argv, 5);
    since = clock();

    bh->max_points = (int)strtol(argv[5], NULL, BASE);
    bh->oldest = bh->newest = bh->buffers_full = 0;
    bh->colour = (rgb){0.0F, 0.5F, 0.0F};
    bh->track = calloc((size_t)bh->max_points, sizeof (components));
    bh->track[bh->newest] = to_xyz(bh);
    bh->ball_size = 0.1F;
    bh->view_radius = 20.0F;
    bh->view_longitude = 0.0F;
    bh->view_latitude = 90.0F;

    // Initialize GLUT
    glutInit(&argc, argv);
    glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB | GLUT_DEPTH);  // Need to double buffer for animation
    // Create and position the graphics window
    glutInitWindowPosition(0, 0);
    glutInitWindowSize(640, 480);
    glutCreateWindow("Black Hole Demo");
    // Initialize GLEW
    glewInit();
    // Initialize OpenGL
    OpenGLInit();
    // Set up callback functions
    glutKeyboardFunc(KeyPressFunc);
    glutSpecialFunc(SpecialKeyFunc);
    glutReshapeFunc(ResizeWindow);
    glutDisplayFunc(Animate);
    // Start the main loop.  glutMainLoop never returns.
    glutMainLoop();
    return(0);          // Compiler requires this to be here. (Never reached)
}
