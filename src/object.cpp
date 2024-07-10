#include "../include/object.h"
#include "../include/constants.h"
#include <cmath>

Object::Object(const Vector2 pos, const Vector2 vel, const double m, const double d, const double red, const double green, const double blue)
: position(pos), velocity(vel), mass(m), density(d), r(red), g(green), b(blue) {
    radius = (std::cbrt((3 * mass) / (4 * 3.141592 * density)) / (1000 * SCALE_FACTOR));
}