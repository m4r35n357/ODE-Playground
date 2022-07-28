#include "real.h"

typedef enum Display {BOTH=0, BALLS=1, LINES=2} display;

/*
 * OpenGL stuff . . .
 */
void OpenGLInit(void);

/*
 * Callbacks
 */
void SpecialKeyFunc(int Key, int x, int y);

void KeyPressFunc(unsigned char Key, int x, int y);

void ResizeWindow(int w, int h);

void Animate(void);

void CloseWindow(void);

/*
 * OSD/HUD
 */
void osd(int x, int y, float r, float g, float b, char *string);
