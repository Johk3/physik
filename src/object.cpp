#include "../include/object.h"
#include "../include/constants.h"
#include <cmath>

Object::Object(const Vector3 pos, const Vector3 vel, const double m, const double d, const double red, const double green, const double blue, ObjectShape shape)
: position(pos), velocity(vel), mass(m), density(d), r(red), g(green), b(blue), shape(shape) {
    radius = (std::cbrt((3 * mass) / (4 * 3.141592 * density)) / (1000 * SCALE_FACTOR));
}
