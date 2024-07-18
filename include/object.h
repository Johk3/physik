#ifndef OBJECT_H
#define OBJECT_H

#include "vector2.h"
#include <vector>
enum class ObjectShape {
    SPHERE,
    TRIANGLE,
    FLAT_SURFACE,
    CONTAINER,
    COW
};
// Represents a physical object in the simulation
class Object {
public:
    Object(const Vector3 pos, const Vector3 vel, const double m, const double d, const double red, const double green, const double blue, ObjectShape shape);

    Vector3 position;
    Vector3 velocity;
    Vector3 acceleration;

    double mass;
    double density;
    double r, g, b;
    double radius;

    std::vector<Vector3> trail;
    double TRAIL_SPACING;
    size_t MAX_TRAIL_LENGTH;

    ObjectShape shape;  // New member to store the object's shape
};

#endif // OBJECT_H