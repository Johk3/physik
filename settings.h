//
// Created by Hugo on 08/07/2024.
//

#ifndef SETTINGS_H
#define SETTINGS_H

// PHYSICS SETTINGS
#define BOUNCE_FACTOR 1.0       // Restitution. Collision stuff
#define CONST_RADIUS 10.0f      // Default radius of ball (og code)
#define SCALE_FACTOR 10         // used to scale all objects
#define MIN_GRAVITY_MASS 1e7    // Minimum mass for an object to have a gravitational field

// TRAIL SETTINGS
#define ENABLE_TRAIL true       // Enable or disable trail
#define SPACING 0.05f          // for trail spacing
#define TRAIL_LENGTH 10         // How many trail "dots" per object

// GRAPHICS SETTINGS
#define CIRCLE_SEGMENTS 32      // Number of segments in circle rendering
#define WINDOW_HEIGHT 1440
#define WINDOW_WIDTH 1440

// PERFORMANCE SETTINGS
#define ENABLE_THREADING true
#define REFRESH_RATE 5000       // Hz. Improves physics accuracy but will slow time down

#endif //SETTINGS_H
