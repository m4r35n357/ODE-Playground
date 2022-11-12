/*
 * (c) 2018-2022 m4r35n357@gmail.com (Ian Smith), for licencing see the LICENCE file
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

/*
 * Global variables
 */
extern controls *c;

extern trail *t;

typedef enum Display {BOTH=0, POSITION=1, TRAIL=2} display;

typedef enum CName {YELLOW=0, CYAN=1, MAGENTA=2, RED=3, GREEN=4, BLUE=5, DARK_GREY=6, LIGHT_GREY=7, GREY=8,
                    DARK_YELLOW=9, DARK_CYAN=10, DARK_MAGENTA=11, DARK_RED=12, DARK_GREEN=13, DARK_BLUE=14} colour_name;

extern display mode;

extern char hud[];

extern clock_t since;

extern _Bool finished, paused, stepping, running, osd_active, solid;

extern float elapsed, cpu, ball_scale;

extern int length, oldest, newest, colour_index, mesh;

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
void SetupView (void);

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
void buffer_point (void);

/*
 * Re-draw boilerplate (no logic)
 */
void ReDraw (void);
