//
// Created by Hugo on 08/07/2024.
//

#ifndef SETTINGS_H
#define SETTINGS_H

// PHYSICS SETTINGS
#define BOUNCE_FACTOR 1.0       // Restitution. Collision stuff
#define SIMULATION_ACCURACY 0.3     // Low value -> High accuracy and bad fps

// TRAIL SETTINGS
#define ENABLE_TRAIL true       // Enable or disable trail
#define SPACING 0.003f          // for trail spacing
#define TRAIL_LENGTH 50         // How many trail "dots" per object
#define MAX_TRAIL_SIZE 1.0      // The max size of the trail compared to the main object

// OTHER SETTINGS
#define CONST_RADIUS 10.0f      // Default radius of ball (og code)
#define SCALE_FACTOR 6          // used to scale all objects
#define REFRESH_RATE 5000       // Hz. Improves physics accuracy but will slow time down
#define CIRCLE_SEGMENTS 32      // Number of segments in circle rendering
#define WINDOW_HEIGHT 1440
#define WINDOW_WIDTH 1440
#define ENABLE_THREADING true


#endif //SETTINGS_H
