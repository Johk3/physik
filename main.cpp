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
    obj1.mass = 100.0f;
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
    obj3.position.y = 0.5f;
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

    // Loop until the user closes the window
    while (!glfwWindowShouldClose(window)) {



        // Clear the screen to black
        glClearColor(0.0f, 0.0f, 0.0f, 1.0f);
        glClear(GL_COLOR_BUFFER_BIT);

        // Update the gravity acceleration of both objects
        std::vector<Object> allObjects = {obj1, obj2, obj3};
        obj1.acceleration = gravity(obj1, allObjects);
        obj2.acceleration = gravity(obj2, allObjects);

        // Time taken between frames
        float delta_time = 1.0f / REFRESH_RATE;

        // Calculate velocity from acceleration
        obj1.velocity.x = obj1.velocity.x + obj1.acceleration.x * delta_time;
        obj1.velocity.y = obj1.velocity.y + obj1.acceleration.y * delta_time;
        obj2.velocity.x = obj2.velocity.x + obj2.acceleration.x * delta_time;
        obj2.velocity.y = obj2.velocity.y + obj2.acceleration.y * delta_time;

        // Calculate position form velocity
        obj1.position.x = obj1.position.x + obj1.velocity.x * delta_time;
        obj1.position.y = obj1.position.y + obj1.velocity.y * delta_time;
        obj2.position.x = obj2.position.x + obj2.velocity.x * delta_time;
        obj2.position.y = obj2.position.y + obj2.velocity.y * delta_time;


        // Draw both objects and their trails

        // OBJECT 1

        // Store current position in trail
        obj1.trail[obj1.trailIndex] = obj1.position;
        obj1.trailIndex = (obj1.trailIndex + 1) % TRAIL_LENGTH;

        // Draw trail
        drawTrail(obj1);

        // Draw current position (main dot)
        mat4x4 mvp;
        mat4x4_identity(mvp);
        mat4x4_translate(mvp, obj1.position.x, obj1.position.y, 0.0f);

        // apply the mvp matrix
        glMatrixMode(GL_PROJECTION);
        glLoadIdentity();
        glMultMatrixf((const GLfloat*)mvp);
        // Draw a white point
        drawDude(obj1.r,obj1.g,obj1.b, 10.0f, 0.0f, 0.0f);

        // OBJECT 2

        // Store current position in trail
        obj2.trail[obj2.trailIndex] = obj2.position;
        obj2.trailIndex = (obj2.trailIndex + 1) % TRAIL_LENGTH;

        // Draw trail
        drawTrail(obj2);

        // Draw current position (main dot)
        mat4x4_identity(mvp);
        mat4x4_translate(mvp, obj2.position.x, obj2.position.y, 0.0f);

        // apply the mvp matrix
        glMatrixMode(GL_PROJECTION);
        glLoadIdentity();
        glMultMatrixf((const GLfloat*)mvp);
        // Draw a white point
        drawDude(obj2.r,obj2.g,obj2.b, 10.0f, 0.0f, 0.0f);

        // OBJECT 3

        // Store current position in trail
        obj3.trail[obj3.trailIndex] = obj3.position;
        obj3.trailIndex = (obj3.trailIndex + 1) % TRAIL_LENGTH;

        // Draw trail
        drawTrail(obj3);

        // Draw current position (main dot)
        mat4x4_identity(mvp);
        mat4x4_translate(mvp, obj3.position.x, obj3.position.y, 0.0f);

        // apply the mvp matrix
        glMatrixMode(GL_PROJECTION);
        glLoadIdentity();
        glMultMatrixf((const GLfloat*)mvp);
        // Draw a white point
        drawDude(obj3.r,obj3.g,obj3.b, 10.0f, 0.0f, 0.0f);

        // Swap front and back buffers
        glfwSwapBuffers(window);

        // Poll for and process events
        glfwPollEvents();

        // Refresh rate pause
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
