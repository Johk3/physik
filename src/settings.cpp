// src/Settings.cpp
#include "../include/settings.h"
#include "../include/constants.h"
#include <algorithm>

namespace Settings {
    int g_circleSegments;
    double g_trailSpacing;
    double g_maxTrailLength;
    double g_bounceFactor;
    double g_scaleFactor;
    double g_G;
    double g_MAX_FORCE;
    double g_arrowProportionality;
    int g_refreshRate;
    bool g_enableTrail;
    bool g_drawArrow;
    bool g_drawShadow;
    bool g_enableLOD;
    bool g_simulate = true;
    bool g_enableRotation = true;
    const double EPSILON = 1e-6;

    void initialize() {
        g_circleSegments = CIRCLE_SEGMENTS;
        g_trailSpacing = SPACING;
        g_maxTrailLength = TRAIL_LENGTH;
        g_bounceFactor = BOUNCE_FACTOR;
        g_scaleFactor = SCALE_FACTOR;
        g_G = GRAVITY;
        g_MAX_FORCE = MAX_FORCE_INIT;
        g_refreshRate = REFRESH_RATE;
        g_enableTrail = ENABLE_TRAIL;
        g_drawArrow = DRAW_ARROW;
        g_drawShadow = DRAW_SHADOW;
        g_enableLOD = ENABLE_LOD;
        g_simulate = SIMULATE;
        g_enableRotation = true;
    }

    void update() {
        // This function can be used to update settings based on user input
        // It will be called from the control panel
        g_circleSegments = std::clamp(g_circleSegments, 3, 100);
        g_trailSpacing = std::clamp(g_trailSpacing, 0.001, 0.1);
        g_maxTrailLength = std::clamp(g_maxTrailLength, 5.0, 1000.0);
        g_bounceFactor = std::clamp(g_bounceFactor, 0.0, 10.0);
        g_scaleFactor = std::clamp(g_scaleFactor, 0.1, 100.0);
        g_G = std::clamp(g_G, 1e-20, 1e10);
        g_MAX_FORCE = std::clamp(g_MAX_FORCE, 0.01, 1e10);
        g_refreshRate = std::clamp(g_refreshRate, 1, 100000);
        g_arrowProportionality = std::clamp(g_arrowProportionality, 1.0, 25.0);
    }
}