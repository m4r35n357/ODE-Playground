/*
 *  N-Body simulator
 *
 * Example: ./h-nbody-gl  0 10 .01 10000
 *
 * (c) 2018-2022 m4r35n357@gmail.com (Ian Smith), for licencing see the LICENCE file
 */

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <string.h>
#include <time.h>
#include <math.h>
#include <GL/freeglut.h>    // OpenGL Graphics Utility Library
#include "taylor-ode.h"
#include "h-nbody-gl.h"

static controls *c;
static body *nb;
static void *p;

static display d;
static char hud[128];
static clock_t since;
static double elapsed, cpu;

static GLenum finished = GL_FALSE;
static GLenum stopped = GL_FALSE;
static GLenum running = GL_TRUE;
static GLenum stepping = GL_FALSE;

static void Key_up (void) {
    nb->view_latitude += 1.0F;
}

static void Key_down (void) {
    nb->view_latitude -= 1.0F;
}

static void Key_left (void) {
    nb->view_longitude += 1.0F;
}

static void Key_right (void) {
    nb->view_longitude -= 1.0F;
}

// glutKeyboardFunc is called below to set this function to handle all normal key presses.
static void KeyPressFunc (unsigned char Key, int x, int y) { (void)x; (void)y;
    switch (Key) {
        case 'A':
        case 'a':
            nb->radius -= 0.1F;
            break;
        case 'Z':
        case 'z':
            nb->radius += 0.1F;
            break;
        case 'R':
        case 'r':
            running = !running;
            stopped = GL_FALSE;
            break;
        case 's':
        case 'S':
            stepping = !stepping;
            stopped = GL_FALSE;
            break;
        case 'F':
        case 'f':
            glutFullScreenToggle();
            break;
        case 27:    // Escape key
            exit(1);
    }
}

// glutSpecialFunc is called below to set this function to handle
//      all special key presses.  See glut.h for the names of special keys.
static void SpecialKeyFunc (int Key, int x, int y) { (void)x; (void)y;
    switch (Key) {
        case GLUT_KEY_UP:
            Key_up();
            break;
        case GLUT_KEY_DOWN:
            Key_down();
            break;
        case GLUT_KEY_LEFT:
            Key_left();
            break;
        case GLUT_KEY_RIGHT:
            Key_right();
            break;
    }
}

static void output(int x, int y, float r, float g, float b, char *string) {
    glColor3f( r, g, b );
    glWindowPos2i(x, y);
    int len = (int)strlen(string);
    for (int i = 0; i < len; i++) {
        glutBitmapCharacter(GLUT_BITMAP_9_BY_15, string[i]);
    }
}

/*
 * Animate() handles the animation and the redrawing of the graphics window contents.
 */
static void Animate (void) {
    // Clear the redering window
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    glLoadIdentity();

    glTranslatef(0.0F, 0.0F, - nb->radius);
    glRotatef(nb->view_latitude, 1.0F, 0.0F, 0.0F);
    glRotatef(nb->view_longitude, 0.0F, 0.0F, 1.0F);

    if (d == BOTH || d == LINES) {
        glBegin( GL_LINE_STRIP );
        for (int k = 0; k < nb->track->newest; k += 1) {
            components point = nb->track->buffer[k];
            glColor3f(nb->colour.x, nb->colour.y, nb->colour.z);
            glVertex3f(point.x, point.y, point.z);
        }
        glEnd();
    }

    if (d == BOTH || d == BALLS) {
        glTranslatef((float)nb->coordinates->x, (float)nb->coordinates->y, (float)nb->coordinates->z);
        glColor3f(nb->colour.x, nb->colour.y, nb->colour.z);
        glutWireSphere((float)0.1F, 10, 10);
    }

    sprintf(hud, "t:%5.1Lf  x:%5.1Lf  y:%5.1Lf  z:%5.1Lf  ",
                  c->step * c->step_size, nb->coordinates->x, nb->coordinates->y, nb->coordinates->z);
    output(10, glutGet(GLUT_WINDOW_HEIGHT) - 20, 0.0F, 0.5F, 0.5F, hud);

    sprintf(hud, "Elapsed: %.1fs  CPU: %.1fs  %.1f %%",
                  elapsed = finished ? elapsed : glutGet(GLUT_ELAPSED_TIME) / 1000.0F,
                  cpu = finished ? cpu : (double)(clock() - since) / CLOCKS_PER_SEC,
                  (100.0F * c->step) / c->steps);
    output(10, 10, 0.0F, 0.5F, 0.5F, hud);

    if (! finished && !stopped) {
        if (tsm_gen(c, nb->coordinates, p)) {
            if (d == BOTH || d == LINES) {
                nb->track->buffer[nb->track->newest++] = *nb->coordinates;
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

// Initialize OpenGL's rendering modes
void OpenGLInit (void) {
    glShadeModel(GL_FLAT);
    glClearColor(0.0F, 0.0F, 0.0F, 0.0F);
    glClearDepth(1.0F);
    glEnable(GL_DEPTH_TEST);
}

// ResizeWindow is called when the window is resized
static void ResizeWindow (int w, int h) {
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

// Set up OpenGL, hook up callbacks, and start the main loop
int main (int argc, char** argv) {
    d = (int)strtol(argv[1], NULL, 10);
    c = get_c(argv);
    p = get_p(argc, argv, c->order);
    since = clock();

    nb = malloc(sizeof (body));
    nb->coordinates = malloc(sizeof (components));
    nb->coordinates->x = strtold(argv[5], NULL);
    nb->coordinates->y = strtold(argv[6], NULL);
    nb->coordinates->z = strtold(argv[7], NULL);
    nb->track = malloc(sizeof (line));
    nb->track->buffer = calloc((size_t)c->steps, sizeof (components));
    nb->track->buffer[nb->track->newest++] = *nb->coordinates;
    nb->track->newest = 0;
    nb->colour = (components) { .x = 0.0F, .y = 0.5F, .z = 0.0F };
    nb->radius = 20.0L;
    nb->view_longitude = 0.0L;
    nb->view_latitude = 0.0L;

    // Need to double buffer for animation
    glutInit(&argc,argv);
    glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB | GLUT_DEPTH );

    // Create and position the graphics window
    glutInitWindowPosition(0, 0);
    glutInitWindowSize(640, 480);
    glutCreateWindow("TSM ODE Demo");

    // Initialize OpenGL.
    OpenGLInit();

    // Set up callback functions for key presses
    glutKeyboardFunc(KeyPressFunc);
    glutSpecialFunc(SpecialKeyFunc);

    // Set up the callback function for resizing windows
    glutReshapeFunc(ResizeWindow);

    // Callback for graphics image redrawing
    glutDisplayFunc(Animate);

    // Start the main loop.  glutMainLoop never returns.
    glutMainLoop();

    return(0);          // Compiler requires this to be here. (Never reached)
}
