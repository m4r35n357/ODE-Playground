/*
 *  N-Body simulator
 *
 * (c) 2018-2022 m4r35n357@gmail.com (Ian Smith), for licencing see the LICENCE file
 */

#include <stdio.h>
#include <string.h>
#include <time.h>
#include <math.h>
#include <GL/glew.h>
#include <GL/freeglut.h>    // OpenGL Graphics Utility Library
#include "taylor-ode.h"
#include "opengl.h"

static controls *c;
static particle *ball;
static void *p;

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
        case GLUT_KEY_UP: ball->view_latitude += 1.0F; break;
        case GLUT_KEY_DOWN: ball->view_latitude -= 1.0F; break;
        case GLUT_KEY_LEFT: ball->view_longitude += 1.0F; break;
        case GLUT_KEY_RIGHT: ball->view_longitude -= 1.0F; break;
    }
}

void KeyPressFunc (unsigned char Key, int x, int y) { (void)x; (void)y;
    switch (Key) {
        case 'B': case 'b': ball->size /= 1.1F; break;
        case 'G': case 'g': ball->size *= 1.1F; break;
        case 'A': case 'a': ball->view_radius -= 0.1F; break;
        case 'Z': case 'z': ball->view_radius += 0.1F; break;
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
    glTranslatef(0.0F, 0.0F, - ball->view_radius);
    glRotatef(ball->view_latitude, 1.0F, 0.0F, 0.0F);
    glRotatef(ball->view_longitude, 0.0F, 0.0F, 1.0F);

    glLightfv(GL_LIGHT0, GL_POSITION, &light_pos[0]);

    if (d == BOTH || d == LINES) {
        glBegin(GL_LINE_STRIP);
        for (int k = 0; k < ball->current; k += 1) {
            glColor3f(ball->colour.a, ball->colour.b, ball->colour.c);
            glVertex3f(ball->track[k].a, ball->track[k].b, ball->track[k].c);
        }
        glEnd();
    }

    if (d == BOTH || d == BALLS) {
        glTranslatef((float)ball->coordinates->x, (float)ball->coordinates->y, (float)ball->coordinates->z);
        glColor3f(ball->colour.a, ball->colour.b, ball->colour.c);
        glutSolidSphere(ball->size, 10, 10);
    }

    sprintf(hud, "t: %.1Lf  x: % .1Lf  y: % .1Lf  z: % .1Lf  ",
                  c->step * c->step_size, ball->coordinates->x, ball->coordinates->y, ball->coordinates->z);
    osd(10, glutGet(GLUT_WINDOW_HEIGHT) - 20, 0.0F, 0.5F, 0.5F, hud);

    sprintf(hud, "Elapsed: %.1fs  CPU: %.1fs  %.1f %%",
                  elapsed = finished ? elapsed : ((float)(glutGet(GLUT_ELAPSED_TIME)) / 1000.0F),
                  cpu = finished ? cpu : (double)(clock() - since) / CLOCKS_PER_SEC,
                  (float)(100 * c->step) / (float)c->steps);
    osd(10, 10, 0.0F, 0.5F, 0.5F, hud);

    if (!finished && !stopped) {
        if (tsm_gen(c, ball->coordinates, p)) {
            if (d == BOTH || d == LINES) {
                ball->track[ball->current++] = (rgb){(float)ball->coordinates->x, (float)ball->coordinates->y, (float)ball->coordinates->z};
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
    p = get_p(argc, argv, c->order);
    since = clock();

    ball = malloc(sizeof (particle));
    ball->coordinates = malloc(sizeof (components));
    ball->coordinates->x = strtold(argv[5], NULL);
    ball->coordinates->y = strtold(argv[6], NULL);
    ball->coordinates->z = strtold(argv[7], NULL);
    ball->current = 0;
    ball->track = calloc((size_t)c->steps, sizeof (components));
    ball->track[ball->current] = (rgb){(float)ball->coordinates->x, (float)ball->coordinates->y, (float)ball->coordinates->z};
    ball->size = 0.1F;
    ball->colour = (rgb){0.0F, 0.5F, 0.0F };
    ball->view_radius = 20.0F;
    ball->view_longitude = 0.0F;
    ball->view_latitude = 90.0F;

    // Initialize GLUT
    glutInit(&argc, argv);
    glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB | GLUT_DEPTH);  // Need to double buffer for animation
    // Create and position the graphics window
    glutInitWindowPosition(0, 0);
    glutInitWindowSize(640, 480);
    glutCreateWindow("TSM ODE Demo");
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
