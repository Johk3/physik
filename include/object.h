#ifndef OBJECT_H
#define OBJECT_H

#include "Vector2.h"
#include <vector>

// Represents a physical object in the simulation
class Object {
public:
    Object(const Vector2 pos, const Vector2 vel, const double m, const double d, const double red, const double green, const double blue);

    Vector2 position;
    Vector2 velocity;
    Vector2 acceleration;

    double mass;
    double density;
    double r, g, b;
    double radius;

    std::vector<Vector2> trail;
    float TRAIL_SPACING;
    size_t MAX_TRAIL_LENGTH;
};

#endif // OBJECT_H