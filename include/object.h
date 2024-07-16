#ifndef OBJECT_H
#define OBJECT_H

#include "vector2.h"
#include <vector>
#include <CL/cl_platform.h>

// Represents a physical object in the simulation
// object.h
class Object {
public:
    Object(const Vector2& pos = {0, 0},
           const Vector2& vel = {0, 0},
           double m = 1.0,
           double d = 1.0,
           double red = 1.0,
           double green = 1.0,
           double blue = 1.0);
    Vector2 position;
    Vector2 velocity;
    Vector2 acceleration;
    double mass;
    double density;
    float radius;
    double r, g, b;
    std::vector<Vector2> trail;
    float TRAIL_SPACING;
    size_t MAX_TRAIL_LENGTH;

    // GPU data struct and methods (as previously defined)
    struct GPUData {
        cl_float2 position;
        cl_float2 velocity;
        cl_float2 acceleration;
        cl_float mass;    // Changed from cl_double to cl_float
        cl_float radius;  // Changed from cl_double to cl_float
    };

    GPUData getGPUData() const;
    void updateFromGPUData(const GPUData& data);
};
#endif // OBJECT_H