//
// Created by Hugo on 08/07/2024.
//

#ifndef SETTINGS_H
#define SETTINGS_H

// PHYSICS SETTINGS
#define BOUNCE_FACTOR 0.8       // Restitution. Collision stuff

// TRAIL SETTINGS
#define ENABLE_TRAIL false      // Enable or disable trail
#define SPACING 0.05f           // for trail spacing
#define TRAIL_LENGTH 30         // How many trail "dots" per object

// OTHER SETTINGS
#define CONST_RADIUS 10.0f      // Default radius of ball (og code)
#define SCALE_FACTOR 6          // used to scale all objects
#define REFRESH_RATE 2500       // Hz. Improves physics accuracy but will slow time down if computer is too slow
#define CIRCLE_SEGMENTS 32      // Number of segments in circle rendering


#endif //SETTINGS_H
