#include "../include/object.h"
#include "../include/constants.h"
#include <cmath>

Object::Object(const Vector3 pos, const Vector3 vel, const double m, const double d, const double red, const double green, const double blue, ObjectShape shape)
: position(pos), velocity(vel), mass(m), density(d), r(red), g(green), b(blue), shape(shape) {
    radius = (std::cbrt((3 * mass) / (4 * 3.141592 * density)) / (1000 * SCALE_FACTOR));
    angularVelocity = Vector3{0, 0, 0};
    rotation = Vector3{0, 0, 0};

    // Calculate moment of inertia based on shape
    switch(shape) {
        case ObjectShape::SPHERE:
            momentOfInertia = 2.0 / 5.0 * mass * radius * radius;
            rotationalFriction = 0.01;  // Low friction for spheres
            break;
        case ObjectShape::TRIANGLE:
            momentOfInertia = mass * (radius * radius) / 6;
            rotationalFriction = 0.03;  // Medium friction for triangles
            break;
        case ObjectShape::FLAT_SURFACE:
        case ObjectShape::CONTAINER:
            momentOfInertia = mass * (radius * radius) / 6;
            rotationalFriction = 0.05;  // Higher friction for flat surfaces and containers
            break;
        case ObjectShape::COW:
            momentOfInertia = mass * radius * radius;
            rotationalFriction = 0.04;  // Medium-high friction for cows
            break;
    }
}

void Object::updateRotation(double deltaTime) {
    // Apply rotational friction
    Vector3 frictionTorque = -angularVelocity * rotationalFriction;
    Vector3 angularAcceleration = frictionTorque / momentOfInertia;

    angularVelocity += angularAcceleration * deltaTime;

    // Update rotation
    rotation += angularVelocity * deltaTime;

    // Normalize rotation to keep it between 0 and 2π
    rotation.x = fmod(rotation.x, 2 * M_PI);
    rotation.y = fmod(rotation.y, 2 * M_PI);
    rotation.z = fmod(rotation.z, 2 * M_PI);

    // Stop rotation if angular velocity becomes very small
    const double minAngularVelocity = 0.01;
    if (angularVelocity.length() < minAngularVelocity) {
        angularVelocity = Vector3{0, 0, 0};
    }
}