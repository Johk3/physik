#include <GLFW/glfw3.h>
#include "linmath.h"
#include <math.h>
#include <cmath>
#include <chrono>
#include <iterator>
#include <thread>
#include <vector>
#include <iostream>
#include <memory>
#include "settings.h"


struct Vector2 {
    double x, y;
    // Vector operations which can be performed locally
    Vector2 operator-(const Vector2& other) const {
        return {x - other.x, y - other.y};
    }


    Vector2 operator+(const Vector2& other) const {
        return {x + other.x, y + other.y};
    }

    Vector2 operator*(float scalar) const {
        return {x * scalar, y * scalar};
    }
    Vector2 operator/(float scalar) const {
        return {x / scalar, y / scalar};
    }

    double length() const {
        return std::sqrt(x*x + y*y);
    }

    Vector2 normal() const {
        return {x / length(), y / length()};
    }

    double dot(Vector2& other) const {

        return (x * other.x + y * other.y);

    }

};

//Physical Object
class Object {
public:
    Object(Vector2 pos, Vector2 vel, double m, double d, double red, double green, double blue)
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

Object createObject(Vector2 pos, Vector2 vel, double m, double d, double red, double green, double blue) {
    return Object(pos, vel, m, d, red, green, blue);
}

void updateObjectTrail(Object& obj) {
    if (obj.trail.empty()) {
        obj.trail.push_back(obj.position);
        return;
    }

    Vector2 lastPosition = obj.trail.back();
    Vector2 diff = obj.position - lastPosition;
    double distance = diff.length();

    if (distance >= Object::TRAIL_SPACING) {
        // Calculate how many new dots we need to add
        int numNewDots = static_cast<int>(distance / Object::TRAIL_SPACING);

        for (int i = 0; i < numNewDots; ++i) {
            double t = (i + 1) * Object::TRAIL_SPACING / distance;
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

bool checkCollision(const Object& obj1, const Object& obj2) {
    double dx = obj1.position.x - obj2.position.x;
    double dy = obj1.position.y - obj2.position.y;
    double distance = std::sqrt(dx*dx + dy*dy);
    return distance < (obj1.radius + obj2.radius);
}

void handleCollision(Object& obj1, Object& obj2) {
    // Calculate vector from obj1 to obj2
    Vector2 delta = obj2.position - obj1.position;
    float distance = delta.length();

    // Normalize the delta vector
    Vector2 normal = delta.normal();

    // Calculate relative velocity
    Vector2 relativeVelocity = obj2.velocity - obj1.velocity;

    // Calculate relative velocity along the normal using dot product
    float velocityAlongNormal = relativeVelocity.dot(normal);

    // Do not resolve if velocities are separating
    if (velocityAlongNormal > 0)
        return;

    // Calculate restitution (bounce factor)
    float restitution = BOUNCE_FACTOR;

    // Calculate impulse scalar
    float impulseScalar = -(1 + restitution) * velocityAlongNormal;
    impulseScalar /= 1/obj1.mass + 1/obj2.mass;

    // Apply impulse
    Vector2 impulse = normal * impulseScalar;

    // Update velocities
    obj1.velocity.x -= impulse.x / obj1.mass;
    obj1.velocity.y -= impulse.y / obj1.mass;
    obj2.velocity.x += impulse.x / obj2.mass;
    obj2.velocity.y += impulse.y / obj2.mass;

    // Separate the objects to prevent overlapping
    float overlap = (obj1.radius + obj2.radius) - distance;
    if (overlap > 0) {
        Vector2 separationVector = normal * (overlap * 0.5);

        obj1.position.x -= separationVector.x * obj2.mass/(obj1.mass + obj2.mass);
        obj1.position.y -= separationVector.y * obj2.mass/(obj1.mass + obj2.mass);
        obj2.position.x += separationVector.x * obj1.mass/(obj1.mass + obj2.mass);
        obj2.position.y += separationVector.y * obj1.mass/(obj1.mass + obj2.mass);
    }
}

// Calculate acceleration due to gravity between two objects, acceleration for first object is returned
Vector2 gravity(Object object1, const std::vector<Object>& objects) {
    const double G = 6.6743e-11f; // Gravitational constant
    Vector2 totalAcceleration = {0.0f, 0.0f};

    for (const auto& object2 : objects) {
        // Skip if it's the same object
        if (&object1 == &object2) continue;

        Vector2 diff = object2.position - object1.position;
        double r2 = diff.x * diff.x + diff.y * diff.y;
        double r = std::sqrt(r2);

        // Avoid division by zero
        if (r < 1e-6f) continue;

        double force = G * object2.mass / (r2 * r);

        // Avoid teleporatation due to obscene force
        double MAX_FORCE = 10000.0f;
        force = std::min(force, MAX_FORCE);

        totalAcceleration = totalAcceleration + diff * force;
    }

    return totalAcceleration;
}

// Updates the physical state of an object over a given time step
void updateObject(Object& obj, const std::vector<Object>& allObjects, double delta_time) {
    obj.acceleration = gravity(obj, allObjects);
    obj.velocity = obj.velocity + obj.acceleration * delta_time;
    obj.position = obj.position + obj.velocity * delta_time;

    // Update the object's trail for rendering motion trail
    updateObjectTrail(obj);
}

// -- Objects which can be drawn
void drawSquare(double x, double y, double r, double g, double b, double alpha, double size) {
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

void drawCircle(double x, double y, double radius, double r, double g, double b, double alpha) {
    int num_segments = CIRCLE_SEGMENTS;
    glColor4d(r, g, b, alpha);
    glBegin(GL_TRIANGLE_FAN);
    glVertex2d(x, y);
    for (int i = 0; i <= num_segments; i++) {
        double theta = 2.0 * M_PI * static_cast<double>(i) / static_cast<double>(num_segments);
        double dx = radius * std::cos(theta);
        double dy = radius * std::sin(theta);
        glVertex2d(x + dx, y + dy);
    }
    glEnd();
}
// -- Objects which can be drawn

// Draws an object and its trail on the screen
void drawObject(const Object& obj) {

    // Draw trail
    if (ENABLE_TRAIL) {
        for (size_t i = 0; i < obj.trail.size(); ++i) {
            double alpha = static_cast<double>(i) / obj.trail.size();
            double trailRadius = obj.radius * 0.5 * alpha;
            drawCircle(obj.trail[i].x, obj.trail[i].y, trailRadius, obj.r, obj.g, obj.b, alpha);
        }
    };

    // Draw main object
    drawCircle(obj.position.x, obj.position.y, obj.radius, obj.r, obj.g, obj.b, 1.0);
}

//Used to set the color of an object
void set_color(Object& obj, double r, double g, double b) {

    obj.r = r;
    obj.g = g;
    obj.b = b;

}

int main() {
    // Initialize GLFW
    if (!glfwInit()) {
        return -1;
    }

    // Create a windowed mode window and its OpenGL context
    GLFWwindow* window = glfwCreateWindow(1000, 1000, "GRAVITY", NULL, NULL);
    if (!window) {
        glfwTerminate();
        return -1;
    }

    // Make the window's context current
    glfwMakeContextCurrent(window);
    // Initialize objects

    std::vector<Object> allObjects;

    // Create and add objects
    for (int i=0; i < 25; i++) {
        for (int j=0; j < 25; j++) {
            allObjects.push_back(createObject({-0.5 + i * 0.04f,  j * 0.04f}, {0.0f, 0.0f}, 1e5, 0.1, 1.0, 1.0, 1.0));  // White object
        }
    }

    allObjects.push_back(createObject({0.02,  -0.5f}, {0.0f, 0.0f}, 5e11, 5e3, 0.0, 0.0, 1.0));

    // Loop until the user closes the window
    while (!glfwWindowShouldClose(window)) {
        glClearColor(0.0f, 0.0f, 0.0f, 1.0f);
        glClear(GL_COLOR_BUFFER_BIT);

        double delta_time = 1.0f / REFRESH_RATE;

        // Update all objects
        for (auto& obj : allObjects) {
            updateObject(obj, allObjects, delta_time);
        }

        // Draw all objects
        for (const auto& obj : allObjects) {
            drawObject(obj);
        }

        // Check for and handle collisions
        for (size_t i = 0; i < allObjects.size(); i++) {
            for(size_t j = i +1; j < allObjects.size(); j++) {
                if (checkCollision(allObjects[i], allObjects[j])) {
                    handleCollision(allObjects[i], allObjects[j]);
                }
            }
        }

        // Swap the back buffer with the front
        glfwSwapBuffers(window);
        // Listen for any events
        glfwPollEvents();

        // Set the refresh rate
        std::this_thread::sleep_for(std::chrono::milliseconds(1000/REFRESH_RATE));
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
