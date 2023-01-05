/*
 *  OpenGL common code
 *
 * (c) 2018-2023 m4r35n357@gmail.com (Ian Smith), for licencing see the LICENCE file
 */

#include <stdio.h>
#include <GL/glew.h>
#include <GL/freeglut.h>
#include "opengl.h"

controls *c;

trail *t;

display mode = BOTH;

char hud[128];

clock_t since;

_Bool finished = 0, paused = 0, stepping = 0, running = 1, osd_active = 1, solid = 1;

int length, oldest = 0, newest = 0, colour_index = DARK_GREEN, mesh = 10;

static float elapsed, cpu, radius = 20.0F, latitude = 90.0F, longitude = 0.0F, ball_size = 0.1F;

void SpecialKeyFunc (int Key, int x, int y) { (void)x; (void)y;
    switch (Key) {
        case        GLUT_KEY_UP: latitude += 1.0F; break;
        case      GLUT_KEY_DOWN: latitude -= 1.0F; break;
        case      GLUT_KEY_LEFT: longitude += 1.0F; break;
        case     GLUT_KEY_RIGHT: longitude -= 1.0F; break;
        case      GLUT_KEY_HOME: radius -= 0.2F; break;
        case       GLUT_KEY_END: radius += 0.2F; break;
        case    GLUT_KEY_INSERT: solid = !solid; break;
        case   GLUT_KEY_PAGE_UP: mesh++; break;
        case GLUT_KEY_PAGE_DOWN: mesh = mesh > 2 ? mesh - 1 : mesh; break;
        default: break;
    }
}

void KeyPressFunc (unsigned char Key, int x, int y) { (void)x; (void)y;
    switch (Key) {
        case 'D': case 'd': colour_index++; break;
        case 'C': case 'c': colour_index--; break;
        case 'G': case 'g': ball_size *= 1.1F; break;
        case 'B': case 'b': ball_size /= 1.1F; break;
        case 'S': case 's': running = !running; paused = 0; break;
        case 'P': case 'p': stepping = !stepping; paused = 0; break;
        case 'F': case 'f': glutFullScreenToggle(); break;
        case 'V': case 'v': mode = (mode + 1) % 3; break;
        case 'H': case 'h': osd_active = !osd_active; break;
        case 'Q': case 'q': case 27: exit(1);  // Code 27 is the Escape key
        default: break;
    }
}

void OpenGLInit () {
    glShadeModel(GL_FLAT);
    glClearColor(0.0F, 0.0F, 0.0F, 0.0F);
    glClearDepth(1.0F);
    glEnable(GL_DEPTH_TEST);
    glEnable(GL_LIGHTING);
    glEnable(GL_LIGHT0);
    glEnable(GL_COLOR_MATERIAL);
}

void ResizeWindow (int w, int h) {
    h = (h == 0) ? 1 : h;
    w = (w == 0) ? 1 : w;
    glViewport(0, 0, w, h);  // View port uses whole window
    glMatrixMode(GL_PROJECTION);  // Set up the projection view matrix (not very well!)
    glLoadIdentity();
    gluPerspective(60.0F, (float)w / (float)h, 1.0F, 100.0F);
    glMatrixMode(GL_MODELVIEW);  // Select the Modelview matrix
}

void ApplicationInit (int argc, char** argv, char *title) {
    glutInit(&argc, argv);  // Initialize GLUT
    glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB | GLUT_DEPTH);  // Need to double buffer for animation
    glutInitWindowPosition(0, 0);  // Create and position the graphics window
    glutInitWindowSize(640, 480);
    glutCreateWindow(title);
    glewInit();  // Initialize GLEW
    OpenGLInit();  // Initialize OpenGL
    glutKeyboardFunc(KeyPressFunc);  // Set up callback functions
    glutSpecialFunc(SpecialKeyFunc);
    glutReshapeFunc(ResizeWindow);
    glutDisplayFunc(Animate);
    fprintf(stderr, "\n  OpenGL: %s\n", glGetString(GL_VERSION));
    fprintf(stderr, "FreeGLUT: %d\n", glutGet(GLUT_VERSION));
    fprintf(stderr, "    GLEW: %s\n\n", glewGetString(GLEW_VERSION));
}

void SetupView () {
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);  // Clear the rendering window
    glLoadIdentity();
    glTranslatef(0.0F, 0.0F, - radius);
    glRotatef(latitude, 1.0F, 0.0F, 0.0F);
    glRotatef(longitude, 0.0F, 0.0F, 1.0F);
    glLightfv(GL_LIGHT0, GL_POSITION, (float []){-100.0F, 100.0F, -100.0F, 0.0F});
}

void ReDraw () {
    glFlush();  // Flush the pipeline, and swap the buffers
    glutSwapBuffers();
    glutPostRedisplay();  // Request a re-draw for animation purposes
}

rgb get_colour (int index) {
    return (rgb []){
        (rgb){1.0F, 1.0F, 0.0F}, (rgb){0.0F, 1.0F, 1.0F}, (rgb){1.0F, 0.0F, 1.0F},
        (rgb){1.0F, 0.0F, 0.0F}, (rgb){0.0F, 1.0F, 0.0F}, (rgb){0.0F, 0.0F, 1.0F},
        (rgb){0.2F, 0.2F, 0.2F}, (rgb){0.8F, 0.8F, 0.8F}, (rgb){0.5F, 0.5F, 0.5F},
        (rgb){0.5F, 0.5F, 0.0F}, (rgb){0.0F, 0.5F, 0.5F}, (rgb){0.5F, 0.0F, 0.5F},
        (rgb){0.5F, 0.0F, 0.0F}, (rgb){0.0F, 0.5F, 0.0F}, (rgb){0.0F, 0.0F, 0.5F}
    }[index % 15];
}

void line_trail (trail *track) {
    glColor3f(track->colour.a, track->colour.b, track->colour.c);
    glBegin(GL_LINE_STRIP);
    for (int i = oldest; i != newest; i = (i + 1) % length) {  // read buffers
        glVertex3f(track->points[i].a, track->points[i].b, track->points[i].c);
    }
    glVertex3f(track->points[newest].a, track->points[newest].b, track->points[newest].c);
    glEnd();
}

void line_position (point p, rgb colour, float scale) {
    glColor3f(0.3F, 0.3F, 0.3F);
    glBegin(GL_LINES);
    glVertex3f(0.0F, 0.0F, 0.0F);
    glVertex3f(p.a, p.b, p.c);
    glEnd();
    glColor3f(colour.a, colour.b, colour.c);
    glPushMatrix();
    glTranslatef(p.a, p.b, p.c);
    solid ? glutSolidSphere(ball_size * scale, mesh, mesh) : glutWireSphere(ball_size * scale, mesh, mesh);
    glPopMatrix();
}

void osd (int x, int y, char *string) {
    glWindowPos2i(x, y);
    glutBitmapString(GLUT_BITMAP_9_BY_15, (const unsigned char *)string);
}

void osd_summary () {
    sprintf(hud, "Elapsed: %.1fs  CPU: %.1fs  %.0f%%",
                  elapsed = finished ? elapsed : 0.001F * (float)glutGet(GLUT_ELAPSED_TIME),
                  cpu = finished ? cpu : (float)(clock() - since) / CLOCKS_PER_SEC,
                  (float)(100.0L * c->step / c->steps));
    osd(10, 10, hud);
}

void buffer_point () {
    static _Bool full = 0;
    newest++;
    if (!full && newest == length) full = 1;
    if (full) {
        oldest = (newest + 1) % length;
        newest %= length;
    }
}
