/*
 * (c) 2018-2025 m4r35n357@gmail.com (Ian Smith), for licencing see the LICENCE file
 */
#pragma once
#include <time.h>
#include "real.h"

/*
 * Triple of floats
 */
typedef struct triple_f {
    float a, b, c;
} rgb, point;

/*
 * Particle/Body tracks
 */
typedef struct Trail {
    struct triple_f colour, *points;
} trail;

typedef enum Display {BOTH=0, POSITION=1, TRAIL=2} display;

/*
 * Global variables
 */
extern controls *c;

extern trail *t;

extern display mode;

extern char hud[];

extern clock_t since;

extern bool finished, paused, stepping, running, osd_active, solid;

extern int length, oldest, newest, colour_index, mesh;

/*
 * OpenGL set-up functions
 */
void OpenGLInit (void);

void ApplicationInit (int argc, char **argv, char *title);

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
void SetupView (void);

/*
 * Colours
 */
rgb get_colour (int index);

/*
 * Lines & balls
 */
void line_trail (trail *track);

void line_position (point p, rgb colour, float scale);

/*
 * OSD/HUD
 */
void osd (int x, int y, char *string);

void osd_summary (void);

/*
 * Extract current coordinates from data
 */
point get_current_point (void *data);

/*
 * Push latest point to the track buffer
 */
void buffer_point (void);

/*
 * Re-draw boilerplate (no logic)
 */
void ReDraw (void);
