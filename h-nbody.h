
#include "real.h"

/*
 * Calculate the centre of mass of the system
 */
void cog (nbody *nb);

/*
 * Hamiltonian
 */
real h (nbody *nb);

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
 * OSD
 */
void output(int x, int y, float r, float g, float b, char *string);
