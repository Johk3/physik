#include "../include/object.h"
#include "../include/constants.h"
#include "../include/settings.h"
#include <cmath>

Object::Object(const Vector2& pos, const Vector2& vel, double m, double d, double red, double green, double blue)
    : position(pos),
      velocity(vel),
      acceleration{0.0, 0.0},  // Initialize acceleration to zero
      mass(m),
      density(d),
      r(red),
      g(green),
      b(blue),
      TRAIL_SPACING(Settings::g_trailSpacing),
      MAX_TRAIL_LENGTH(static_cast<size_t>(Settings::g_maxTrailLength))
{
    radius = static_cast<float>(std::cbrt((3 * mass) / (4 * M_PI * density)) / (1000 * Settings::g_scaleFactor));

    // Initialize an empty trail
    trail.reserve(MAX_TRAIL_LENGTH);
}

Object::GPUData Object::getGPUData() const {
    return {
            {static_cast<cl_float>(position.x), static_cast<cl_float>(position.y)},
            {static_cast<cl_float>(velocity.x), static_cast<cl_float>(velocity.y)},
            {static_cast<cl_float>(acceleration.x), static_cast<cl_float>(acceleration.y)},
            static_cast<cl_float>(mass),
            static_cast<cl_float>(radius)
        };
}

void Object::updateFromGPUData(const GPUData& data) {
    position = {data.position.x, data.position.y};
    velocity = {data.velocity.x, data.velocity.y};
    acceleration = {data.acceleration.x, data.acceleration.y};
    // Note: mass and radius are not updated as they should remain constant
}