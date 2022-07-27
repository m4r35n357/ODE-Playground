/*
 *  N-Body simulator
 *
 * (c) 2018-2022 m4r35n357@gmail.com (Ian Smith), for licencing see the LICENCE file
 */

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <string.h>
#include <time.h>
#include <math.h>
#include <GL/glew.h>
#include <GL/freeglut.h>    // OpenGL Graphics Utility Library
#include "symplectic.h"
#include "h-nbody.h"

static controls *c;
static nbody *nb;

static display d;
static char hud[128];
static clock_t since;
static double elapsed, cpu;

static float light_pos[] = { -100.0f, 100.0f, -100.0f, 0.0f };

static GLenum finished = GL_FALSE;
static GLenum stopped = GL_FALSE;
static GLenum running = GL_TRUE;
static GLenum stepping = GL_FALSE;

// glutSpecialFunc is called below to set this function to handle special key presses - see glut.h for the names of special keys
void SpecialKeyFunc (int Key, int x, int y) { (void)x; (void)y;
    switch (Key) {
        case GLUT_KEY_UP: nb->view_latitude += 1.0F; break;
        case GLUT_KEY_DOWN: nb->view_latitude -= 1.0F; break;
        case GLUT_KEY_LEFT: nb->view_longitude += 1.0F; break;
        case GLUT_KEY_RIGHT: nb->view_longitude -= 1.0F; break;
    }
}

// glutKeyboardFunc is called below to set this function to handle all normal key presses.
void KeyPressFunc (unsigned char Key, int x, int y) { (void)x; (void)y;
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

void output(int x, int y, float r, float g, float b, char *string) {
    glColor3f(r, g, b);
    glWindowPos2i(x, y);
    int len = (int)strlen(string);
    for (int i = 0; i < len; i++) {
        glutBitmapCharacter(GLUT_BITMAP_9_BY_15, string[i]);
    }
}

/*
 * Animate() handles the animation and the redrawing of the graphics window contents.
 */
void Animate (void) {
    // Clear the rendering window
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    glLoadIdentity();
    glTranslatef(0.0F, 0.0F, - nb->view_radius);
    glRotatef(nb->view_latitude, 1.0F, 0.0F, 0.0F);
    glRotatef(nb->view_longitude, 0.0F, 0.0F, 1.0F);

    glLightfv(GL_LIGHT0, GL_POSITION, &light_pos[0]);

    if (d == BOTH || d == LINES) {
        body *b = nb->bodies;
        for (int i = 0; i < nb->n; i += 1) {
            glBegin( GL_LINE_STRIP );
            for (int k = nb->oldest; k != nb->newest; k = (k + 1) % nb->max_points) {
                components point = b[i].track->buffer[k];
                glColor3f(b[i].colour.r, b[i].colour.g, b[i].colour.b);
                glVertex3f((float)point.x, (float)point.y, (float)point.z);
            }
            glEnd();
        }
    }

    if (d == BOTH || d == BALLS) {
        body *b = nb->bodies;
        glTranslatef((float)(b[0].x - nb->centre.x), (float)(b[0].y - nb->centre.y), (float)(b[0].z - nb->centre.z));
        glColor3f(b[0].colour.r, b[0].colour.g, b[0].colour.b);
        glutSolidSphere((float)powl(nb->ball_scale * b[0].m, 1.0L/3.0L), 10, 10);
        for (int i = 1; i < nb->n; i += 1) {
            glTranslatef((float)(b[i].x - b[i - 1].x), (float)(b[i].y - b[i - 1].y), (float)(b[i].z - b[i - 1].z));
            glColor3f(b[i].colour.r, b[i].colour.g, b[i].colour.b);
            glutSolidSphere((float)powl(nb->ball_scale * b[i].m, 1.0L/3.0L), 10, 10);
        }
    }

    sprintf(hud, "t: %.1Lf  h: %.6Le  ~sf: %.1Lf",
                  c->step * c->step_size, nb->h, error(nb->h - nb->h0));
    output(10, glutGet(GLUT_WINDOW_HEIGHT) - 20, 0.0F, 0.5F, 0.5F, hud);

    sprintf(hud, "Elapsed: %.1fs  CPU: %.1fs  %.1f %%",
                  elapsed = finished ? elapsed : ((float)(glutGet(GLUT_ELAPSED_TIME)) / 1000.0F),
                  cpu = finished ? cpu : (double)(clock() - since) / CLOCKS_PER_SEC,
                  (float)(100 * c->step) / (float)c->steps);
    output(10, 10, 0.0F, 0.5F, 0.5F, hud);

    if (!finished && !stopped) {
        if (generate(c, nb)) {
            cog(nb);
            nb->h = h(nb) > nb->h ? h(nb) : nb->h;
            if (d == BOTH || d == LINES) {
                nb->newest += 1;
                if (!nb->full && (nb->newest == nb->max_points)) {
                    nb->full = 1;
                }
                if (nb->full) {
                    nb->oldest = (nb->newest + 1) % nb->max_points;
                    nb->newest %= nb->max_points;
                }
                body *b = nb->bodies;
                for (int i = 0; i < nb->n; i += 1) {
                    b[i].track->buffer[nb->newest] = (components) { b[i].x, b[i].y, b[i].z };
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

void CloseWindow (void) {
    fprintf(stderr, "H : %+.18Le\n", nb->h);
}

// Initialize OpenGL's rendering modes
void OpenGLInit (void) {
    glShadeModel(GL_FLAT);
    glClearColor(0.0F, 0.0F, 0.0F, 0.0F);
    glClearDepth(1.0F);
    glEnable(GL_DEPTH_TEST);
    glEnable(GL_LIGHTING);
    glEnable(GL_LIGHT0);
    glEnable(GL_COLOR_MATERIAL);
}

// ResizeWindow is called when the window is resized
void ResizeWindow (int w, int h) {
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
    d = (display)strtol(argv[1], NULL, BASE);
    c = get_c(argv);
    nb = (nbody *)get_p(argc, argv, (argc - 6) / 7);
    body *b = nb->bodies;
    for (int i = 0; i < nb->n; i += 1) {
        b[i].track = malloc(sizeof (line));
        b[i].track->buffer = calloc((size_t)c->steps, sizeof (components));
        b[i].track->buffer[nb->newest] = (components) { b[i].x, b[i].y, b[i].z };
    }

    fprintf(stderr, "\n");
    fprintf(stderr, "H0: %+.18Le\n", nb->h);
    since = clock();

    // Need to double buffer for animation
    glutInit(&argc,argv);
    glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB | GLUT_DEPTH);

    // Create and position the graphics window
    glutInitWindowPosition(0, 0);
    glutInitWindowSize(640, 480);
    glutCreateWindow("N-Body Demo");

    glewInit();

    // Initialize OpenGL.
    OpenGLInit();

    // Set up callback functions
    glutKeyboardFunc(KeyPressFunc);
    glutSpecialFunc(SpecialKeyFunc);
    glutReshapeFunc(ResizeWindow);
    glutDisplayFunc(Animate);
    glutCloseFunc(CloseWindow);

    // Start the main loop.  glutMainLoop never returns.
    glutMainLoop();

    return(0);          // Compiler requires this to be here. (Never reached)
}
