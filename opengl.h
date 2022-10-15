#include <time.h>
#include "real.h"

/*
 * Global variables
 */
static controls *c;

typedef enum Display {BOTH=0, BALLS=1, LINES=2} display;

static display mode;
static char hud[128];
static clock_t since;
static double elapsed, cpu;

static float light_position[] = { -100.0F, 100.0F, -100.0F, 0.0F };

static _Bool finished = 0, stopped = 0, stepping = 0, running = 1;

/*
 * OpenGL set-up functions
 */
void OpenGLInit (void);

void ApplicationInit (int argc, char** argv, char *title);

/*
 * Callbacks
 */
void SpecialKeyFunc (int Key, int x, int y);

void KeyPressFunc (unsigned char Key, int x, int y);

void ResizeWindow (int w, int h);

void Animate (void);

void CloseWindow (void);

/*
 * View/Camera boilerplate (no logic)
 */
void SetupView (float view_radius, float view_latitude, float view_longitude, float *light_pos);

/*
 * OSD/HUD
 */
void osd (int x, int y, float r, float g, float b, char *string);

/*
 * Push latest point to the track buffer
 */
void buffer_point (int max_points, int *oldest, int *newest, _Bool *full);

/*
 * Re-draw boilerplate (no logic)
 */
void ReDraw (void);
