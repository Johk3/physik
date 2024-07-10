#ifndef CONTROLPANEL_H
#define CONTROLPANEL_H

#include <GLFW/glfw3.h>

// Control panel related functions and variables
extern double g_trailSpacing;
extern double g_maxTrailLength;
extern double g_bounceFactor;
extern double g_scaleFactor;
extern double g_G;
extern double g_MAX_FORCE;
extern int g_circleSegments;
extern int g_refreshRate;
extern int g_enable_trail;
extern int g_enable_lod;
extern int g_simulate;

void updateFPS();
void renderControlPanel(GLFWwindow* controlWindow);

#endif // CONTROL_PANEL_H