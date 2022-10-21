/*
 *  OpenGL common code
 *
 * (c) 2018-2022 m4r35n357@gmail.com (Ian Smith), for licencing see the LICENCE file
 */

#include <string.h>
#include <GL/glew.h>
#include <GL/freeglut.h>
#include "opengl.h"

/*
 * Global variables
 */
controls *c;

display mode = BOTH;

char hud[128];

clock_t since;

double elapsed, cpu;

static float light_position[] = { -100.0F, 100.0F, -100.0F, 0.0F };

_Bool finished = 0, stopped = 0, stepping = 0, running = 1;

float ball_scale = 0.1F;

static float view_radius = 20.0F, view_latitude = 90.0F, view_longitude = 0.0F;

int max_points, oldest = 0, newest = 0, colour_index;

void SpecialKeyFunc (int Key, int x, int y) { (void)x; (void)y;
    switch (Key) {
        case    GLUT_KEY_UP: view_latitude  += 1.0F; break;
        case  GLUT_KEY_DOWN: view_latitude  -= 1.0F; break;
        case  GLUT_KEY_LEFT: view_longitude += 1.0F; break;
        case GLUT_KEY_RIGHT: view_longitude -= 1.0F; break;
        case  GLUT_KEY_HOME: view_radius    -= 1.0F; break;
        case   GLUT_KEY_END: view_radius    += 1.0F; break;
    }
}

void KeyPressFunc (unsigned char Key, int x, int y) { (void)x; (void)y;
    switch (Key) {
        case 'D': case 'd': colour_index += 1; break;
        case 'C': case 'c': colour_index -= 1; break;
        case 'G': case 'g': ball_scale *= 1.1F; break;
        case 'B': case 'b': ball_scale /= 1.1F; break;
        case 'R': case 'r': running = !running; stopped = 0; break;
        case 'S': case 's': stepping = !stepping; stopped = 0; break;
        case 'F': case 'f': glutFullScreenToggle(); break;
        case 'V': case 'v': mode = (mode + 1) % 3; break;
        case  27: exit(1); // Escape key
    }
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

void ApplicationInit (int argc, char** argv, char *title) {
    // Initialize GLUT
    glutInit(&argc, argv);
    glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB | GLUT_DEPTH);  // Need to double buffer for animation
    // Create and position the graphics window
    glutInitWindowPosition(0, 0);
    glutInitWindowSize(640, 480);
    glutCreateWindow(title);
    // Initialize GLEW
    glewInit();
    // Initialize OpenGL
    OpenGLInit();
    // Set up callback functions
    glutKeyboardFunc(KeyPressFunc);
    glutSpecialFunc(SpecialKeyFunc);
    glutReshapeFunc(ResizeWindow);
    glutDisplayFunc(Animate);
}

void SetupView () {
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT); // Clear the rendering window
    glLoadIdentity();
    glTranslatef(0.0F, 0.0F, - view_radius);
    glRotatef(view_latitude, 1.0F, 0.0F, 0.0F);
    glRotatef(view_longitude, 0.0F, 0.0F, 1.0F);
    glLightfv(GL_LIGHT0, GL_POSITION, light_position);
}

rgb get_colour (int index) {
    rgb colours[] = {
        (rgb){1.0F, 1.0F, 0.0F}, (rgb){0.0F, 1.0F, 1.0F}, (rgb){1.0F, 0.0F, 1.0F},
        (rgb){1.0F, 0.0F, 0.0F}, (rgb){0.0F, 1.0F, 0.0F}, (rgb){0.0F, 0.0F, 1.0F},
        (rgb){0.2F, 0.2F, 0.2F}, (rgb){0.8F, 0.8F, 0.8F}, (rgb){0.5F, 0.5F, 0.5F},
        (rgb){0.5F, 0.5F, 0.0F}, (rgb){0.0F, 0.5F, 0.5F}, (rgb){0.5F, 0.0F, 0.5F},
        (rgb){0.5F, 0.0F, 0.0F}, (rgb){0.0F, 0.5F, 0.0F}, (rgb){0.0F, 0.0F, 0.5F}
    };
    return colours[index % 15];
}

void osd (int x, int y, char *string) {
    glWindowPos2i(x, y);
    int len = (int)strlen(string);
    for (int i = 0; i < len; i++) {
        glutBitmapCharacter(GLUT_BITMAP_9_BY_15, string[i]);
    }
}

void buffer_point () {
    static _Bool buffers_full = 0;
    newest += 1;
    if (! buffers_full && (newest == max_points)) {
        buffers_full = 1;
    }
    if (buffers_full) {
        oldest = (newest + 1) % max_points;
        newest %= max_points;
    }
}

void ReDraw (void) {
    glFlush();                  // Flush the pipeline, and swap the buffers
    glutSwapBuffers();
    glutPostRedisplay();        // Request a re-draw for animation purposes
}
