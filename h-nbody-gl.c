/*
 *  gcc -pedantic -Wall -Wshadow -Wpointer-arith -Wcast-qual -Wstrict-prototypes -Wmissing-prototypes -Wextra -Wconversion -Wredundant-decls -Wmissing-field-initializers -Wmissing-declarations -Wuninitialized -Wunsuffixed-float-constants -frounding-math -fsignaling-nans symplectic.c h-nbody.c h-nbody-gl.c -lm -lglut -lGLU -lGL

 *
 * Solar.c
 *
 * Program to demonstrate how to use a local
 * coordinate method to position parts of a
 * model in relation to other model parts.
 *
 * Draws a simple solar system, with a sun, planet and moon.
 * Based on sample code from the OpenGL programming guide
 *      by Woo, Neider, Davis.  Addison-Wesley.
 *
 * Author: Samuel R. Buss
 *
 * Software accompanying the book
 *      3D Computer Graphics: A Mathematical Introduction with OpenGL,
 *      by S. Buss, Cambridge University Press, 2003.
 *
 * Software is "as-is" and carries no warranty.  It may be used without
 *   restriction, but if you modify it, please change the filenames to
 *   prevent confusion between different versions.
 * Bug reports: Sam Buss, sbuss@ucsd.edu.
 * Web page: http://math.ucsd.edu/~sbuss/MathCG
 *
 * USAGE:
 *    Press "r" key to toggle (off and on) running the animation
 *    Press "s" key to single-step animation
 *    The up arrow key and down array key control the
 *          time step used in the animation rate.  Each key
 *          press multiplies or divides the times by a factor
 *          of two (2).
 *    Press ESCAPE to exit.
 *
 */

#include <stdlib.h>
#include <GL/freeglut.h>    // OpenGL Graphics Utility Library
#include "math.h"
#include "symplectic.h"
#include "h-nbody.h"
#include "h-nbody-gl.h"

static controls *c;
static nbody *p;

static GLenum stopped = GL_FALSE;
static GLenum running = GL_TRUE;
static GLenum stepping = GL_FALSE;

// These three variables control the animation's state and speed.
static float AnimateIncrement = 24.0F;  // Time step for animation (hours)

static void Key_up (void) {
    AnimateIncrement *= 2.0F;           // Double the animation time step
}

static void Key_down (void) {
    AnimateIncrement /= 2.0F;           // Halve the animation time step
}

// glutKeyboardFunc is called below to set this function to handle all normal key presses.
static void KeyPressFunc (unsigned char Key, int x, int y) { (void)x; (void)y;
    switch (Key) {
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
    }
}

/*
 * Animate() handles the animation and the redrawing of the graphics window contents.
 */
static void Animate (void) {
    // Clear the redering window
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    nbody *nb = (nbody *)p;

    // Clear the current matrix (Modelview)
    glLoadIdentity();
    glTranslatef((float)nb->bodies[0].q_x, (float)nb->bodies[0].q_y, (float)nb->bodies[0].q_z -20.0F);
    glColor3f(1.0F, 1.0F, 0.0F);
    glutWireSphere(0.4F, 10, 10);

    glLoadIdentity();
    glTranslatef((float)nb->bodies[0].q_x, (float)nb->bodies[0].q_y, (float)nb->bodies[0].q_z -20.0F);
    glColor3f(0.0F, 1.0F, 1.0F);
    glutWireSphere(0.4F, 10, 10);

    if (!stopped) {
        if (!(nb = (nbody *)generate(c, nb))) glutLeaveMainLoop();
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
    p = (nbody *)get_p(argc, argv, 5);

    // Need to double buffer for animation
    glutInit(&argc,argv);
    glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB | GLUT_DEPTH );

    // Create and position the graphics window
    glutInitWindowPosition(0, 0);
    glutInitWindowSize(640, 480);
    glutCreateWindow("Simulator Demo");

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
