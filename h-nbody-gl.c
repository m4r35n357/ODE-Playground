/*
 *  N-Body simulator
 *
 * Example:  ./h-nbody-gl  0 10 .01 10000  .05  100.0 0.0 0.0 0.0 0.0 0.0 0.0 2.0 0.0 4.5 0.4 -0.2 0.0 1.8 3.0 -6.0 0.0 -0.4 0.0 -2.0 1.0 5.0 3.0 0.0 -0.2 0.0 5.8 -0.2 4.0 0.0 -4.0 0.1 -3.6 0.0 0.2 3.0 -4.0 0.0 -0.1 0.0 -0.2 -2.6 3.0 8.0 0.0 -0.3 0.0 2.0 -0.2 4.0 0.0 4.0 -0.2 -4.8 0.0 -0.2
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
#include "symplectic.h"
#include "h-nbody.h"
#include "h-nbody-gl.h"

static controls *c;
static nbody *nb;

static display d;
static char hud[128];
static clock_t since;
static double elapsed, cpu;

static GLenum finished = GL_FALSE;
static GLenum stopped = GL_FALSE;
static GLenum running = GL_TRUE;
static GLenum stepping = GL_FALSE;

// glutSpecialFunc is called below to set this function to handle special key presses - see glut.h for the names of special keys
static void SpecialKeyFunc (int Key, int x, int y) { (void)x; (void)y;
    switch (Key) {
        case GLUT_KEY_UP: nb->view_latitude += 1.0F; break;
        case GLUT_KEY_DOWN: nb->view_latitude -= 1.0F; break;
        case GLUT_KEY_LEFT: nb->view_longitude += 1.0F; break;
        case GLUT_KEY_RIGHT: nb->view_longitude -= 1.0F; break;
    }
}

// glutKeyboardFunc is called below to set this function to handle all normal key presses.
static void KeyPressFunc (unsigned char Key, int x, int y) { (void)x; (void)y;
    switch (Key) {
        case 'B': case 'b': nb->ball_scale /= 1.1F; break;
        case 'G': case 'g': nb->ball_scale *= 1.1F; break;
        case 'A': case 'a': nb->view_radius -= 0.1F; break;
        case 'Z': case 'z': nb->view_radius += 0.1F; break;
        case 'R': case 'r': running = !running; stopped = GL_FALSE; break;
        case 'S': case 's': stepping = !stepping; stopped = GL_FALSE; break;
        case 'F': case 'f': glutFullScreenToggle(); break;
        case 27: exit(1); // Escape key
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

    glTranslatef(0.0F, 0.0F, - nb->view_radius);
    glRotatef(nb->view_latitude, 1.0F, 0.0F, 0.0F);
    glRotatef(nb->view_longitude, 0.0F, 0.0F, 1.0F);

    if (d == BOTH || d == LINES) {
        for (int i = 0; i < nb->n; i += 1) {
            body *b = &nb->bodies[i];
            glBegin( GL_LINE_STRIP );
            for (int k = 0; k < b->track->newest; k += 1) {
                components point = b->track->buffer[k];
                glColor3f((float)b->colour.r, (float)b->colour.g, (float)b->colour.b);
                glVertex3f((float)point.x, (float)point.y, (float)point.z);
            }
            glEnd();
        }
    }

    if (d == BOTH || d == BALLS) {
        body *b = nb->bodies;

        glTranslatef((float)(b[0].q_x - nb->centre.x),
                     (float)(b[0].q_y - nb->centre.y),
                     (float)(b[0].q_z - nb->centre.z));
        glColor3f(b[0].colour.r, b[0].colour.g, b[0].colour.b);
        glutWireSphere((float)powl(nb->ball_scale * b[0].m, 1.0L/3.0L), 10, 10);

        glTranslatef((float)(b[1].q_x - b[0].q_x),
                     (float)(b[1].q_y - b[0].q_y),
                     (float)(b[1].q_z - b[0].q_z));
        glColor3f(b[1].colour.r, b[1].colour.g, b[1].colour.b);
        glutWireSphere((float)powl(nb->ball_scale * b[1].m, 1.0L/3.0L), 10, 10);

        glTranslatef((float)(b[2].q_x - b[1].q_x),
                     (float)(b[2].q_y - b[1].q_y),
                     (float)(b[2].q_z - b[1].q_z));
        glColor3f(b[2].colour.r, b[2].colour.g, b[2].colour.b);
        glutWireSphere((float)powl(nb->ball_scale * b[2].m, 1.0L/3.0L), 10, 10);

        glTranslatef((float)(b[3].q_x - b[2].q_x),
                     (float)(b[3].q_y - b[2].q_y),
                     (float)(b[3].q_z - b[2].q_z));
        glColor3f(b[3].colour.r, b[3].colour.g, b[3].colour.b);
        glutWireSphere((float)powl(nb->ball_scale * b[3].m, 1.0L/3.0L), 10, 10);

        glTranslatef((float)(b[4].q_x - b[3].q_x),
                     (float)(b[4].q_y - b[3].q_y),
                     (float)(b[4].q_z - b[3].q_z));
        glColor3f(b[4].colour.r, b[4].colour.g, b[4].colour.b);
        glutWireSphere((float)powl(nb->ball_scale * b[4].m, 1.0L/3.0L), 10, 10);

        glTranslatef((float)(b[5].q_x - b[4].q_x),
                     (float)(b[5].q_y - b[4].q_y),
                     (float)(b[5].q_z - b[4].q_z));
        glColor3f(b[5].colour.r, b[5].colour.g, b[5].colour.b);
        glutWireSphere((float)powl(nb->ball_scale * b[5].m, 1.0L/3.0L), 10, 10);

        glTranslatef((float)(b[6].q_x - b[5].q_x),
                     (float)(b[6].q_y - b[5].q_y),
                     (float)(b[6].q_z - b[5].q_z));
        glColor3f(b[6].colour.r, b[6].colour.g, b[6].colour.b);
        glutWireSphere((float)powl(nb->ball_scale * b[6].m, 1.0L/3.0L), 10, 10);

        glTranslatef((float)(b[7].q_x - b[6].q_x),
                     (float)(b[7].q_y - b[6].q_y),
                     (float)(b[7].q_z - b[6].q_z));
        glColor3f(b[7].colour.r, b[7].colour.g, b[7].colour.b);
        glutWireSphere((float)powl(nb->ball_scale * b[7].m, 1.0L/3.0L), 10, 10);
    }

    sprintf(hud, "h0: %.9Le  h: %.9Le", nb->h0, h(nb));
    output(10, glutGet(GLUT_WINDOW_HEIGHT) - 20, 0.0F, 0.5F, 0.5F, hud);

    sprintf(hud, "Elapsed: %.1fs  CPU: %.1fs  %.1f %%",
                  elapsed = finished ? elapsed : ((float)(glutGet(GLUT_ELAPSED_TIME)) / 1000.0F),
                  cpu = finished ? cpu : (double)(clock() - since) / CLOCKS_PER_SEC,
                  (float)(100 * c->step) / (float)c->steps);
    output(10, 10, 0.0F, 0.5F, 0.5F, hud);

    if (!finished && !stopped) {
        if (generate(c, nb)) {
            cog(nb);
            if (d == BOTH || d == LINES) {
                for (int i = 0; i < nb->n; i += 1) {
                    body *b = &nb->bodies[i];
                    b->track->buffer[b->track->newest++] = (components) { b->q_x, b->q_y, b->q_z };
                }
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
    gluPerspective(60.0F, aspectRatio, 1.0F, 30.0F);

    // Select the Modelview matrix
    glMatrixMode(GL_MODELVIEW);
}

// Set up OpenGL, hook up callbacks, and start the main loop
int main (int argc, char** argv) {
    c = get_c(argv);
    d = (display)strtol(argv[1], NULL, 10);
    nb = (nbody *)get_p(argc, argv, (argc - 6) / 7);
    since = clock();

    for (int i = 0; i < nb->n; i += 1) {
        body *b = &nb->bodies[i];
        b->track = malloc(sizeof (line));
        b->track->newest = 0;
        b->track->buffer = calloc((size_t)c->steps, sizeof (components));
        b->track->buffer[b->track->newest] = (components) { b->q_x, b->q_y, b->q_z };
    }

    // Need to double buffer for animation
    glutInit(&argc,argv);
    glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB | GLUT_DEPTH);

    // Create and position the graphics window
    glutInitWindowPosition(0, 0);
    glutInitWindowSize(640, 480);
    glutCreateWindow("N-Body Demo");

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
