// include/Settings.h
#ifndef SETTINGS_H
#define SETTINGS_H

namespace Settings {;
    extern double g_trailSpacing;
    extern double g_maxTrailLength;
    extern double g_bounceFactor;
    extern double g_scaleFactor;
    extern double g_G;
    extern double g_MAX_FORCE;
    extern const double EPSILON;
    extern int g_circleSegments;
    extern int g_refreshRate;
    extern bool g_enableTrail;
    extern bool g_enableLOD;
    extern bool g_simulate;

    void initialize();
    void update();
}

#endif // SETTINGS_H