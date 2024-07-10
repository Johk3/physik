//
// Created by Hugo on 08/07/2024.
//

#ifndef SETTINGS_H
#define SETTINGS_H

// PHYSICS SETTINGS
#define BOUNCE_FACTOR 1.0           // Restitution. Collision stuff
#define CONST_RADIUS 10.0f          // Default radius of ball (og code)
#define SCALE_FACTOR 10             // used to scale all objects
#define MIN_GRAVITY_MASS 1e7        // Minimum mass for an object to have a gravitational field
#define GRAVITY 6.67430e-11
#define MAX_FORCE 10000.0;

// TRAIL SETTINGS
#define ENABLE_TRAIL true           // Enable or disable trail
#define SPACING 0.05f               // for trail spacing
#define TRAIL_LENGTH 10             // How many trail "dots" per object
#define TRAIL_SCALE 0.75            // Trail size as a ratio to object

// GRAPHICS SETTINGS
#define CIRCLE_SEGMENTS 32          // Number of segments in circle rendering
#define WINDOW_HEIGHT 1440
#define WINDOW_WIDTH 1440

// PERFORMANCE SETTINGS
#define REFRESH_RATE 5000           // Hz. Improves physics accuracy but will slow time down
#define ENABLE_LOD true             // Enable or disable LOD
#define MIN_LOD 5                   // Minimum number of circle segments
#define PREVENT_DRAW_DISTANCE 1.0   // Distance for object to not be drawn

#endif //SETTINGS_H
