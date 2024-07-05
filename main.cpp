#include <GLFW/glfw3.h>
#include "linmath.h"
#include <math.h>
#include <cmath>
#include <chrono>
#include <iterator>
#include <thread>
#include <vector>
#include <iostream>
#define TRAIL_RADIUS 2.0f
#define CONST_RADIUS 10.0f
#define SCALE_FACTOR 5 // used to scale all objects
#define REFRESH_RATE 1000 // Hz. Improves physics accuracy or might fuck up if computer too slow


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

    double length() const {
        return std::sqrt(x*x + y*y);
    }
};

//Physical Object
struct Object {
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
    static constexpr double TRAIL_SPACING = 0.02f;  // Fixed distance between trail dots
    static constexpr size_t MAX_TRAIL_LENGTH = 20;  // Maximum number of trail dots
};

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
    Vector2 delta = {obj2.position.x - obj1.position.x, obj2.position.y - obj1.position.y};
    float distance = std::sqrt(delta.x * delta.x + delta.y * delta.y);

    // Normalize the delta vector
    Vector2 normal = {delta.x / distance, delta.y / distance};

    // Calculate relative velocity
    Vector2 relativeVelocity = {
        obj2.velocity.x - obj1.velocity.x,
        obj2.velocity.y - obj1.velocity.y
    };

    // Calculate relative velocity along the normal
    float velocityAlongNormal = relativeVelocity.x * normal.x + relativeVelocity.y * normal.y;

    // Do not resolve if velocities are separating
    if (velocityAlongNormal > 0)
        return;

    // Calculate restitution (bounce factor)
    float restitution = 0.8f;

    // Calculate impulse scalar
    float impulseScalar = -(1 + restitution) * velocityAlongNormal;
    impulseScalar /= 1/obj1.mass + 1/obj2.mass;

    // Apply impulse
    Vector2 impulse = {impulseScalar * normal.x, impulseScalar * normal.y};

    // Update velocities
    obj1.velocity.x -= impulse.x / obj1.mass;
    obj1.velocity.y -= impulse.y / obj1.mass;
    obj2.velocity.x += impulse.x / obj2.mass;
    obj2.velocity.y += impulse.y / obj2.mass;

    // Separate the objects to prevent overlapping
    float overlap = (obj1.radius + obj2.radius) - distance;
    if (overlap > 0) {
        Vector2 separationVector = {
            normal.x * overlap * 0.5f,
            normal.y * overlap * 0.5f
        };
        obj1.position.x -= separationVector.x;
        obj1.position.y -= separationVector.y;
        obj2.position.x += separationVector.x;
        obj2.position.y += separationVector.y;
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

void drawDot(double x, double y, double r, double g, double b, double alpha, double size) {
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

// Draws an object and its trail on the screen
void drawObject(const Object& obj) {
    // Draw trail
    for (size_t i = 0; i < obj.trail.size(); ++i) {
        double alpha = static_cast<double>(i) / obj.trail.size();
        drawDot(obj.trail[i].x, obj.trail[i].y, obj.r, obj.g, obj.b, alpha, TRAIL_RADIUS);
    }

    // Draw main object
    drawDot(obj.position.x, obj.position.y, obj.r, obj.g, obj.b, 1.0f, obj.radius * 1000);
}

// Used to update object radius after mass or density change
void update_radius(Object& obj) {

    obj.radius = (std::cbrt((3 * obj.mass) / (4 * 3.141592 * obj.density)) / (1000 * SCALE_FACTOR));

};

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
    // Make two objects
    Object obj1, obj2, obj3, obj4, obj5, obj6;

    obj1.position.x = 0.6f;
    obj1.velocity.x = -2.0f;
    obj1.mass = 1000.0f;
    obj1.r = 0.0f;
    obj1.b = 0.0f;
    update_radius(obj1);

    obj2.position.x = 0.0f;
    obj2.velocity.x = 2.0f;
    obj2.mass = 1000.0f;
    obj2.density = 0.01;
    obj2.r = 0.3f;
    obj2.g = 0.8f;
    obj2.b = 1.0f;
    update_radius(obj2);


    obj3.position.x = 0.3f;
    obj3.position.y = 0.2f;
    obj3.mass = 10000000000.0f;
    obj3.density = 1000.0;
    obj3.g = 0.0f;
    obj3.b = 0.0f;
    update_radius(obj3);

    std::cout << obj3.radius;

    obj4.position.x = 0.2f;
    obj4.position.y = 0.4f;
    obj4.velocity.x = -1.0f;
    obj4.velocity.y = 0.0f;
    obj4.acceleration.x = 0.0f;
    obj4.acceleration.y = 0.0f;
    obj4.mass = 10.0f;
    obj4.r = 0.0f;
    obj4.g = 1.0f;
    obj4.b = 0.0f;
    obj4.radius = CONST_RADIUS/1000;

    obj5.position.x = 0.1f;
    obj5.position.y = -0.5f;
    obj5.velocity.x = 0.0f;
    obj5.velocity.y = 0.0f;
    obj5.acceleration.x = 0.0f;
    obj5.acceleration.y = 0.0f;
    obj5.mass = 100.0f;
    obj5.r = 0.3f;
    obj5.g = 0.8f;
    obj5.b = 1.0f;
    obj5.radius = CONST_RADIUS/1000;

    obj6.position.x = -0.3f;
    obj6.position.y = 0.0f;
    obj6.velocity.x = 0.0f;
    obj6.velocity.y = 0.0f;
    obj6.acceleration.x = 0.0f;
    obj6.acceleration.y = 0.0f;
    obj6.mass = 10000.0f;
    obj6.r = 1.0f;
    obj6.g = 0.0f;
    obj6.b = 0.0f;
    obj6.radius = CONST_RADIUS/1000;
    std::vector<Object> allObjects = {obj1, obj2, obj3};

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
