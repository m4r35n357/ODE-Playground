/*
 *  OpenGL common code
 *
 * (c) 2018-2022 m4r35n357@gmail.com (Ian Smith), for licencing see the LICENCE file
 */

#include <string.h>
#include <GL/glew.h>
#include <GL/freeglut.h>
#include "opengl.h"

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

void SetupView (float view_radius, float view_latitude, float view_longitude, float *light_pos) {
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT); // Clear the rendering window
    glLoadIdentity();
    glTranslatef(0.0F, 0.0F, - view_radius);
    glRotatef(view_latitude, 1.0F, 0.0F, 0.0F);
    glRotatef(view_longitude, 0.0F, 0.0F, 1.0F);
    glLightfv(GL_LIGHT0, GL_POSITION, light_pos);
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

void osd (int x, int y, float r, float g, float b, char *string) {
    glColor3f(r, g, b);
    glWindowPos2i(x, y);
    int len = (int)strlen(string);
    for (int i = 0; i < len; i++) {
        glutBitmapCharacter(GLUT_BITMAP_9_BY_15, string[i]);
    }
}

void buffer_point (int max_points, int *oldest, int *newest, _Bool *full) {
    *newest += 1;
    if (! *full && (*newest == max_points)) {
        *full = 1;
    }
    if (*full) {
        *oldest = (*newest + 1) % max_points;
        *newest %= max_points;
    }
}

void ReDraw (void) {
    glFlush();                  // Flush the pipeline, and swap the buffers
    glutSwapBuffers();
    glutPostRedisplay();        // Request a re-draw for animation purposes
}
