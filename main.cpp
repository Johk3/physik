#include <GLFW/glfw3.h>
#include "linmath.h"
#include <cmath>
#include <chrono>
#include <iterator>
#include <thread>
#include <vector>
#include <iostream>
#include <memory>
//#include <bits/stl_algo.h>

#include "settings.h"

struct Vector2 {
    double x, y;
    // Vector operations which can be performed locally
    // constexpr is for compile-time constant use
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
        return {x / length(), y / length()};
    }

    // Dot product
    [[nodiscard]] constexpr double dot(const Vector2& other) const {
        return (x * other.x + y * other.y);
    }

};

//Physical Object
class Object {
public:
    Object(const Vector2 pos, const Vector2 vel, const double m, const double d, const double red, const double green, const double blue)
        : position(pos), velocity(vel), mass(m), density(d), r(red), g(green), b(blue) {
        radius = (std::cbrt((3 * mass) / (4 * 3.141592 * density)) / (1000 * SCALE_FACTOR));
    }
    // Physical location and movement. Default center, no movement
    Vector2 position = {0.0f, 0.0f};
    Vector2 velocity = {0.0f, 0.0f};
    Vector2 acceleration = {0.0f, 0.0f};

    // Properties. Default 100kg, white dot
    double mass = 100.0;
    double density = 1.0;
    double r = 1.0;
    double g = 1.0;
    double b = 1.0;
    double radius = CONST_RADIUS / (1000 * SCALE_FACTOR);

    // Trail
    std::vector<Vector2> trail;
    static constexpr double TRAIL_SPACING = SPACING;  // Fixed distance between trail dots
    static constexpr size_t MAX_TRAIL_LENGTH = TRAIL_LENGTH;  // Maximum number of trail dots
};


void updateObjectTrail(Object& obj) {
    if (obj.trail.empty()) {
        obj.trail.push_back(obj.position);
        return;
    }

    const Vector2 lastPosition = obj.trail.back();
    const Vector2 diff = obj.position - lastPosition;
    double distance = diff.length();

    if (distance >= Object::TRAIL_SPACING) {
        // Calculate how many new dots we need to add
        const int numNewDots = static_cast<int>(distance / Object::TRAIL_SPACING);

        for (int i = 0; i < numNewDots; ++i) {
            const double t = (i + 1) * Object::TRAIL_SPACING / distance;
            Vector2 newDotPosition = {
                lastPosition.x + diff.x * t,
                lastPosition.y + diff.y * t
            };
            obj.trail.push_back(newDotPosition);

            // Remove oldest dot if we exceed the maximum trail length
            if (obj.trail.size() > Object::MAX_TRAIL_LENGTH) {
                obj.trail.erase(obj.trail.begin());
            }
        }
    }
}

// Check if two circular objects are overlapping
bool checkCollision(const Object& obj1, const Object& obj2) {
    const double distance = (obj1.position - obj2.position).length();
    return distance <= (obj1.radius + obj2.radius);
}

void handleCollision(Object& obj1, Object& obj2) {
    // Calculate vector from obj1 to obj2
    const Vector2 delta = obj2.position - obj1.position;
    const float distance = delta.length();

    // Normalize the delta vector
    const Vector2 normal = delta.normal();

    // Calculate relative velocity
    const Vector2 relativeVelocity = obj2.velocity - obj1.velocity;

    // Calculate relative velocity along the normal using dot product
    float velocityAlongNormal = relativeVelocity.dot(normal);

    // Do not resolve if velocities are separating
    if (velocityAlongNormal > 0)
        return;

    // Calculate impulse scalar
    float impulseScalar = - (1 + BOUNCE_FACTOR) * velocityAlongNormal;
    impulseScalar /= 1 / obj1.mass  +  1 / obj2.mass;

    // Apply impulse
    Vector2 impulse = normal * impulseScalar;

    // Update velocities
    obj1.velocity.x -= impulse.x / obj1.mass;
    obj1.velocity.y -= impulse.y / obj1.mass;
    obj2.velocity.x += impulse.x / obj2.mass;
    obj2.velocity.y += impulse.y / obj2.mass;

    // Separate the objects to prevent overlapping
    const float overlap = (obj1.radius + obj2.radius) - distance;
    if (overlap > 0) {
        Vector2 separationVector = normal * overlap;

        obj1.position.x -= separationVector.x * obj2.mass/(obj1.mass + obj2.mass);
        obj1.position.y -= separationVector.y * obj2.mass/(obj1.mass + obj2.mass);
        obj2.position.x += separationVector.x * obj1.mass/(obj1.mass + obj2.mass);
        obj2.position.y += separationVector.y * obj1.mass/(obj1.mass + obj2.mass);
    }
}

// Calculate acceleration due to gravity between two objects, acceleration for first object is returned
Vector2 gravity(const Object& object1, const std::vector<Object>& objects) {
    Vector2 totalAcceleration = {0.0, 0.0};

    for (const auto& object2 : objects) {
        if (&object1 == &object2) continue;

        Vector2 diff = object2.position - object1.position;
        double r2 = diff.x * diff.x + diff.y * diff.y;

        if (r2 < 1e-12) continue;  // Avoid division by zero

        double r = std::sqrt(r2);

        if (r <= object1.radius + object2.radius) {
            handleCollision(const_cast<Object&>(object1), const_cast<Object&>(object2));
            continue;
        }

        double force = std::min(6.67430e-11 * object2.mass / (r2 * r), 10000.0);
        totalAcceleration = totalAcceleration + diff * force;
    }

    return totalAcceleration;
}

// Updates the physical state of an object over a given time step
void updateObject(Object& obj, std::vector<Object>& all_objects, const double delta_time) {
    obj.acceleration = gravity(obj, all_objects);
    obj.velocity = obj.velocity + obj.acceleration * delta_time;
    obj.position = obj.position + obj.velocity * delta_time;

    // Update the object's trail for rendering motion trail
    updateObjectTrail(obj);
}

// -- Objects which can be drawn
void drawSquare(const double x, const double y, const double r, const double g, const double b, const double alpha, const double size) {
    mat4x4 mvp;
    mat4x4_identity(mvp);
    mat4x4_translate(mvp, x, y, 0.0f);

    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    glMultMatrixf((const GLfloat*)mvp);

    glColor4f(r, g, b, alpha);
    glPointSize(size);
    glBegin(GL_POINTS);
    glVertex2f(0.0f, 0.0f);
    glEnd();
}

void drawCircle(const double x, const double y, const double radius, const double r, const double g, const double b, const double alpha) {
    glColor4d(r, g, b, alpha);
    glBegin(GL_TRIANGLE_FAN);
    glVertex2d(x, y);
    for (int i = 0; i <= CIRCLE_SEGMENTS; i++) {
        const double theta = 2.0 * M_PI * static_cast<double>(i) / static_cast<double>(CIRCLE_SEGMENTS);
        const double dx = radius * std::cos(theta);
        const double dy = radius * std::sin(theta);
        glVertex2d(x + dx, y + dy);
    }
    glEnd();
}
// -- Objects which can be drawn

// Draws an object and its trail on the screen
void drawObject(const Object& obj) {

    // Draw trail
    if constexpr (ENABLE_TRAIL) {
        for (size_t i = 0; i < obj.trail.size(); ++i) {
            double alpha = static_cast<double>(i) / obj.trail.size();
            double trailRadius = obj.radius * 0.5 * alpha;
            drawCircle(obj.trail[i].x, obj.trail[i].y, trailRadius, obj.r, obj.g, obj.b, alpha);
        }
    };

    // Draw main object
    drawCircle(obj.position.x, obj.position.y, obj.radius, obj.r, obj.g, obj.b, 1.0);
}

std::vector<Object> get_objects() {
    std::vector<Object> all_objects;

    // Create and add objects
    for (int i=0; i < 25; i++) {
        for (int j=0; j < 25; j++) {
            const float r = ((rand() % 51) / 10) + 0.5;
            const float g = ((rand() % 51) / 10) + 0.5;
            const float b = ((rand() % 51) / 10) + 0.5;

            all_objects.push_back(Object({-0.5 + i * 0.04f,  j * 0.04f}, {(rand()%3-2)*r*b, (rand()%3-2) * g*r}, 1e5, 0.1, r, g, b));  // White object
        }
    }

    all_objects.push_back(Object({0.02,  -0.5f}, {0.0f, 0.0f}, 5e12, 5e3, 0.0, 0.0, 1.0));

    return all_objects;

};

// Render physics aka calculate physics
void render_physics(std::vector<Object>& all_objects, const double delta_time) {

    // param all_objects: list of objects, all objects to be calculated
    // param delta_time, double. Time taken from last update. Smaller values gives more accuracy in physics


    // Update all objects
    for (auto& obj : all_objects) {
        updateObject(obj, all_objects, delta_time);
    }

};


// Renders the screen window based on objects on it
void render_screen(const std::vector<Object>& all_objects, GLFWwindow* window) {

    // Clear screen?
    glClearColor(0.0f, 0.0f, 0.0f, 1.0f);
    glClear(GL_COLOR_BUFFER_BIT);

    // Draw all objects
    for (const auto& obj : all_objects) {
        drawObject(obj);
    }

    // Swap the back buffer with the front
    glfwSwapBuffers(window);
    // Listen for any events
    glfwPollEvents();

};

int main() {
    // Initialize GLFW
    if (!glfwInit()) {
        return -1;
    }

    // Create a windowed mode window and its OpenGL context
    auto window_deleter = [](GLFWwindow* w) { glfwDestroyWindow(w); glfwTerminate(); };
    std::unique_ptr<GLFWwindow, decltype(window_deleter)> window(
        glfwCreateWindow(WINDOW_WIDTH, WINDOW_HEIGHT, "GRAVITY", nullptr, nullptr),
        window_deleter
    );
    if (!window) {
        glfwTerminate();
        return -1;
    }

    // Make the window's context current
    glfwMakeContextCurrent(window.get());

    // Initialize objects, time, and grid
    std::vector<Object> all_objects = get_objects();
    const double delta_time = 1.0f / REFRESH_RATE;

    // Window loop. Loop until the user closes the window
    while (!glfwWindowShouldClose(window.get())) {

        render_physics(all_objects, delta_time);
        render_screen(all_objects, window.get());

    }

    glfwTerminate();
    return 0;
}

/*
 * arity: 3
 * C(I(Z,Z), C(P1,P3,P2), Z)
 *
 * C(I(Z,Z), C(P1,P3,P2), Z)(x,y,z)
 * I(Z,Z)(C(P1,P3,P2), Z)(x,y,z)
 * I(Z,Z)(P1(P3,P2), Z)(x,y,z)
 * I(Z,Z)(P1(P3,P2)(x,y,z), Z(x,y,z))
 * BASE CASE: I(Z,Z)(0, Z(x,y,z)) = Z(Z(x,y,z)) = 0
 * RECURSIVE CASE: I(Z,Z)(P1(P3,P2)(x,y,z), Z(x,y,z))
 * = Z(I(Z,Z)(P1(P3,P2)(x,y,z), Z(x,y,z), P1(P3,P2)(x,y,z), Z(x,y,z)) = 0,
*/
