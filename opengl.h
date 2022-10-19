/*
 * (c) 2018-2022 m4r35n357@gmail.com (Ian Smith), for licencing see the LICENCE file
 */

#pragma once

#include <time.h>
#include "real.h"

/*
 * Particle/Body tracks
 */
typedef struct Track {
    struct triple_f colour, *points;
} track;

/*
 * Global variables
 */
extern controls *c;

typedef enum Display {BOTH=0, BALLS=1, LINES=2} display;

typedef enum CName {YELLOW=0, CYAN=1, MAGENTA=2, RED=3, GREEN=4, BLUE=5, DARK_GREY=6, LIGHT_GREY=7, GREY=8,
                    DARK_YELLOW=9, DARK_CYAN=10, DARK_MAGENTA=11, DARK_RED=12, DARK_GREEN=13, DARK_BLUE=14} colour_name;

extern display mode;

extern char hud[];

extern clock_t since;

extern double elapsed, cpu;

extern float light_position[];

extern _Bool finished, stopped, stepping, running, buffers_full;

extern float ball_scale, view_radius, view_latitude, view_longitude;

extern int max_points, oldest, newest, colour_index;

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
 * Colours
 */
rgb get_colour (int index);

/*
 * OSD/HUD
 */
void osd (int x, int y, char *string);

/*
 * Extract current coordinates from model
 */
point point_from_model (void *model);

/*
 * Push latest point to the track buffer
 */
void buffer_point (int max_points, int *oldest, int *newest, _Bool *full);

/*
 * Re-draw boilerplate (no logic)
 */
void ReDraw (void);
