#include <GLFW/glfw3.h>
#include "linmath.h"
#include <math.h>
#include <cmath>
#include <chrono>
#include <iterator>
#include <thread>
#include <vector>
#define TRAIL_LENGTH 20
#define CONST_RADIUS 10.0f
#define REFRESH_RATE 1000 // Hz. Improves physics accuracy or might fuck up if computer too slow



// Position for object
struct Position {
    double x;
    double y;
};

// Velocity for object
struct Velocity {
    double x;
    double y;
};

// Acceleration for object
struct Acceleration {
    double x;
    double y;
};

//Physical Object
struct Object {
    // Physical location and movement
    Position position;
    Velocity velocity;
    Acceleration acceleration;
    // Properties
    double mass;
    double r;
    double g;
    double b;
    double radius;
    // Trail
    Position trail[TRAIL_LENGTH];
    int trailIndex;
};

bool checkCollision(const Object& obj1, const Object& obj2) {
    double dx = obj1.position.x - obj2.position.x;
    double dy = obj1.position.y - obj2.position.y;
    double distance = std::sqrt(dx*dx + dy*dy);
    return distance < (obj1.radius + obj2.radius);
}

void handleCollision(Object& obj1, Object& obj2) {
    // Calculate normal vector
    double nx = obj2.position.x - obj1.position.x;
    double ny = obj2.position.y - obj1.position.y;
    double dist = std::sqrt(nx*nx + ny*ny);
    nx /= dist;
    ny /= dist;

    // Calculate relative velocity
    double vx = obj1.velocity.x - obj2.velocity.x;
    double vy = obj1.velocity.y - obj2.velocity.y;

    // Calculate relative velocity in terms of the normal direction
    double velocityAlongNormal = vx * nx + vy * ny;

    // Do not resolve if velocities are separating
    if (velocityAlongNormal > 0)
        return;

    // Calculate restitution (bounce factor)
    double restitution = 0.8f;

    // Calculate impulse scalar
    double impulseScalar = -(1 + restitution) * velocityAlongNormal;
    impulseScalar /= 1/obj1.mass + 1/obj2.mass;

    // Apply impulse
    double impulseX = impulseScalar * nx;
    double impulseY = impulseScalar * ny;

    obj1.velocity.x += impulseX / obj1.mass;
    obj1.velocity.y += impulseY / obj1.mass;
    obj2.velocity.x -= impulseX / obj2.mass;
    obj2.velocity.y -= impulseY / obj2.mass;

    // Separate the objects to prevent sticking
    double overlap = (obj1.radius + obj2.radius) - dist;
    double separationX = overlap * nx * 0.5f;
    double separationY = overlap * ny * 0.5f;
    obj1.position.x -= separationX;
    obj1.position.y -= separationY;
    obj2.position.x += separationX;
    obj2.position.y += separationY;

    // Add a small random perturbation to prevent objects from getting stuck in a perfect bounce loop
    double perturbation = 0.01f;
    obj1.velocity.x += perturbation * (rand() / (float)RAND_MAX - 0.5f);
    obj1.velocity.y += perturbation * (rand() / (float)RAND_MAX - 0.5f);
    obj2.velocity.x += perturbation * (rand() / (float)RAND_MAX - 0.5f);
    obj2.velocity.y += perturbation * (rand() / (float)RAND_MAX - 0.5f);
}

// Draw a cube
void drawDude(GLfloat r, GLfloat g, GLfloat b, GLfloat rad, GLfloat posx, GLfloat posy) {
    glColor3f(r,g,b);
    glPointSize(rad);
    glBegin(GL_POINTS);
    glVertex2f(posx, posy);
    glEnd();
}



// Calculate acceleration due to gravity between two objects, acceleration for first object is returned
Acceleration gravity(Object object1, const std::vector<Object>& objects) {
    const double G = 6.6743e-11f; // Gravitational constant
    Acceleration totalAcceleration = {0.0f, 0.0f};

    for (const auto& object2 : objects) {
        // Skip if it's the same object
        if (&object1 == &object2) continue;

        double dx = object2.position.x - object1.position.x;
        double dy = object2.position.y - object1.position.y;
        double r2 = dx*dx + dy*dy;
        double r = std::sqrt(r2);

        // Avoid division by zero
        if (r < 1e-6f) continue;

        double force = G * object2.mass / (r2 * r);

        // Avoid teleporatation due to obscene force
        double MAX_FORCE = 1000.0f;
        force = std::min(force, MAX_FORCE);

        totalAcceleration.x += force * dx;
        totalAcceleration.y += force * dy;
    }

    return totalAcceleration;
}

// Updates the physical state of an object over a given time step
void updateObject(Object& obj, const std::vector<Object>& allObjects, float delta_time) {
    obj.acceleration = gravity(obj, allObjects);
    obj.velocity.x += obj.acceleration.x * delta_time;
    obj.velocity.y += obj.acceleration.y * delta_time;
    obj.position.x += obj.velocity.x * delta_time;
    obj.position.y += obj.velocity.y * delta_time;

    // Update the object's trail for rendering motion trail
    obj.trail[obj.trailIndex] = obj.position;
    obj.trailIndex = (obj.trailIndex + 1) % TRAIL_LENGTH;
}
// Draw trail
void drawTrail(Object obj) {
    for (int i = 0; i < TRAIL_LENGTH; i++) {
        int index = (obj.trailIndex - i - 1 + TRAIL_LENGTH) % TRAIL_LENGTH;
        double alpha = 1.0f - (double)i / TRAIL_LENGTH;

        mat4x4 mvp;
        mat4x4_identity(mvp);
        mat4x4_translate(mvp, obj.trail[index].x, obj.trail[index].y, 0.0f);

        glMatrixMode(GL_PROJECTION);
        glLoadIdentity();
        glMultMatrixf((const GLfloat*)mvp);

        drawDude(obj.r, obj.g, obj.b, 5.0f * alpha, 0.0f, 0.0f);
    }
}

// Draws an object and its trail on the screen
void drawObject(const Object& obj) {
    drawTrail(obj);

    // Set up the model-view-projection matrix for the object
    mat4x4 mvp;
    mat4x4_identity(mvp);
    mat4x4_translate(mvp, obj.position.x, obj.position.y, 0.0f);

    // Apply the transformation and draw the object
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    glMultMatrixf((const GLfloat*)mvp);
    drawDude(obj.r, obj.g, obj.b, CONST_RADIUS, 0.0f, 0.0f);
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

    // Make two objects
    Object obj1;
    Object obj2;
    Object obj3;

    obj1.position.x = 0.5f;
    obj1.position.y = 0.6f;
    obj1.velocity.x = -2.0f;
    obj1.velocity.y = 0.0f;
    obj1.acceleration.x = 0.0f;
    obj1.acceleration.y = 0.0f;
    obj1.mass = 10.0f;
    obj1.r = 0.0f;
    obj1.g = 0.0f;
    obj1.b = 1.0f;
    obj1.radius = CONST_RADIUS/1000;
    obj1.trailIndex = 0;
    obj1.trail[obj1.trailIndex] = obj1.position;

    obj2.position.x = 0.5f;
    obj2.position.y = 0.5f;
    obj2.velocity.x = 0.0f;
    obj2.velocity.y = 0.0f;
    obj2.acceleration.x = 0.0f;
    obj2.acceleration.y = 0.0f;
    obj2.mass = 10000000000.0f;
    obj2.r = 1.0f;
    obj2.g = 0.0f;
    obj2.b = 0.0f;
    obj2.radius = CONST_RADIUS/1000;
    obj2.trailIndex = 0;
    obj2.trail[obj2.trailIndex] = obj2.position;

    obj3.position.x = 0.3f;
    obj3.position.y = 0.2f;
    obj3.velocity.x = 0.0f;
    obj3.velocity.y = 0.0f;
    obj3.acceleration.x = 0.0f;
    obj3.acceleration.y = 0.0f;
    obj3.mass = 10000.0f;
    obj3.r = 1.0f;
    obj3.g = 0.0f;
    obj3.b = 0.0f;
    obj3.radius = CONST_RADIUS/1000;
    obj3.trailIndex = 0;
    obj3.trail[obj2.trailIndex] = obj2.position;
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
