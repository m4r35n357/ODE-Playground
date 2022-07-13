/*
 gcc -pedantic -Wall -Wshadow -Wpointer-arith -Wcast-qual -Wstrict-prototypes -Wmissing-prototypes -Wextra -Wconversion -Wredundant-decls -Wmissing-field-initializers -Wmissing-declarations -Wuninitialized -Wunsuffixed-float-constants -frounding-math -fsignaling-nans symplectic.c dual.c h-newton.c h-newton-gl.c -lm -lglut -lGLU -lGL

 * Example: ./h-newton-gl  6 8 1 1000  1 12 .6
 *
 * (c) 2018-2022 m4r35n357@gmail.com (Ian Smith), for licencing see the LICENCE file
 */

#include <stdlib.h>
#include <GL/freeglut.h>    // OpenGL Graphics Utility Library
#include "math.h"
#include "symplectic.h"
#include "h-newton.h"
#include "h-newton-gl.h"

static controls *c;
static parameters *p;

static GLenum stopped = GL_FALSE;
static GLenum running = GL_TRUE;
static GLenum stepping = GL_FALSE;

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

/*
 * Animate() handles the animation and the redrawing of the graphics window contents.
 */
static void Animate (void) {
    // Clear the redering window
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    // Clear the current matrix (Modelview)
    glLoadIdentity();
    glTranslatef(0.0F, 0.0F, -20.0F);
    glColor3f(1.0F, 0.0F, 0.0F);
    glutWireSphere(1.0F, 10, 10);

    glLoadIdentity();
    glTranslatef((float)(p->q_r * cosl(p->q_phi)), (float)(p->q_r * sinl(p->q_phi)), -20.0F);
    glColor3f(0.0F, 1.0F, 0.0F);
    glutWireSphere(0.4F, 10, 10);

    if (!stopped) {
        if (!(p = (parameters *)generate(c, p))) glutLeaveMainLoop();
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
    p = (parameters *)get_p(argc, argv, 5);

    // Need to double buffer for animation
    glutInit(&argc,argv);
    glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB | GLUT_DEPTH );

    // Create and position the graphics window
    glutInitWindowPosition(0, 0);
    glutInitWindowSize(640, 480);
    glutCreateWindow("Newtonian 2-Body Demo");

    // Initialize OpenGL.
    OpenGLInit();

    // Set up callback functions for key presses
    glutKeyboardFunc(KeyPressFunc);

    // Set up the callback function for resizing windows
    glutReshapeFunc(ResizeWindow);

    // Callback for graphics image redrawing
    glutDisplayFunc(Animate);

    // Start the main loop.  glutMainLoop never returns.
    glutMainLoop();

    return(0);          // Compiler requires this to be here. (Never reached)
}
