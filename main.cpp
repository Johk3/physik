#include <GLFW/glfw3.h>
#include "linmath.h"
#include <math.h>
#include <cmath>
#include <chrono>
#include <iterator>
#include <thread>
#include <vector>
#define TRAIL_LENGTH 20
#define REFRESH_RATE 1000 // Hz. Improves physics accuracy or might fuck up if computer too slow



// Position for object
struct Position {
    float x;
    float y;
};

// Velocity for object
struct Velocity {
    float x;
    float y;
};

// Acceleration for object
struct Acceleration {
    float x;
    float y;
};

//Physical Object
struct Object {
    // Physical location and movement
    Position position;
    Velocity velocity;
    Acceleration acceleration;
    // Properties
    float mass;
    float r;
    float g;
    float b;
    // Trail
    Position trail[TRAIL_LENGTH];
    int trailIndex;
};

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
    const float G = 6.6743e-11f; // Gravitational constant
    Acceleration totalAcceleration = {0.0f, 0.0f};

    for (const auto& object2 : objects) {
        // Skip if it's the same object
        if (&object1 == &object2) continue;

        float dx = object2.position.x - object1.position.x;
        float dy = object2.position.y - object1.position.y;
        float r2 = dx*dx + dy*dy;
        float r = std::sqrt(r2);

        // Avoid division by zero
        if (r < 1e-6f) continue;

        float force = G * object2.mass / (r2 * r);

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
        float alpha = 1.0f - (float)i / TRAIL_LENGTH;

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
    drawDude(obj.r, obj.g, obj.b, 10.0f, 0.0f, 0.0f);
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
    obj2.trailIndex = 0;
    obj2.trail[obj2.trailIndex] = obj2.position;

    obj3.position.x = 0.0f;
    obj3.position.y = -0.3f;
    obj3.velocity.x = 0.0f;
    obj3.velocity.y = 0.0f;
    obj3.acceleration.x = 0.0f;
    obj3.acceleration.y = 0.0f;
    obj3.mass = 10000000000.0f;
    obj3.r = 1.0f;
    obj3.g = 0.0f;
    obj3.b = 0.0f;
    obj3.trailIndex = 0;
    obj3.trail[obj2.trailIndex] = obj2.position;
    std::vector<Object> allObjects = {obj1, obj2, obj3};

    // Loop until the user closes the window
    while (!glfwWindowShouldClose(window)) {
        glClearColor(0.0f, 0.0f, 0.0f, 1.0f);
        glClear(GL_COLOR_BUFFER_BIT);

        float delta_time = 1.0f / REFRESH_RATE;

        // Update all objects
        for (auto& obj : allObjects) {
            updateObject(obj, allObjects, delta_time);
        }

        // Draw all objects
        for (const auto& obj : allObjects) {
            drawObject(obj);
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
