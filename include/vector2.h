// include/Vector2.h
#ifndef VECTOR2_H
#define VECTOR2_H

#include <cmath>

struct Vector2 {
    double x, y;

    // Vector operations
    constexpr Vector2 operator-(const Vector2& other) const {
        return {x - other.x, y - other.y};
    }

    constexpr Vector2 operator+(const Vector2& other) const {
        return {x + other.x, y + other.y};
    }

    constexpr Vector2 operator*(const float scalar) const {
        return {x * scalar, y * scalar};
    }

    constexpr Vector2 operator/(const float scalar) const {
        return {x / scalar, y / scalar};
    }

    // Length of vector, modulus
    [[nodiscard]] constexpr double length() const {
        return std::sqrt(x*x + y*y);
    }

    // Normal
    [[nodiscard]] constexpr Vector2 normal() const {
        double len = length();
        return {x / len, y / len};
    }

    // Dot product
    [[nodiscard]] constexpr double dot(const Vector2& other) const {
        return (x * other.x + y * other.y);
    }
};

struct Vector3 {
    double x, y, z;

    Vector3 operator+(const Vector3& other) const {
        return {x + other.x, y + other.y, z + other.z};
    }

    Vector3 operator-(const Vector3& other) const {
        return {x - other.x, y - other.y, z - other.z};
    }

    Vector3 operator*(double scalar) const {
        return {x * scalar, y * scalar, z * scalar};
    }

    [[nodiscard]] double length() const {
        return std::sqrt(x*x + y*y + z*z);
    }

    [[nodiscard]] Vector3 normal() const {
        double len = length();
        return {x / len, y / len, z / len};
    }

    [[nodiscard]] double dot(const Vector3& other) const {
        return x * other.x + y * other.y + z * other.z;
    }
};


#endif // VECTOR2_H