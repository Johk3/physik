#include <GLFW/glfw3.h>
#include "linmath.h"
#include <math.h>
#include <cmath>
#define TRAIL_LENGTH 20

struct Position {
    float x;
    float y;
};

struct Velocity {
    float x;
    float y;
};

struct Acceleration {
    float x;
    float y;
};

struct Object {
    Position position;
    Velocity velocity;
    Acceleration acceleration;
    float mass;
    float r;
    float g;
    float b;
    Position trail[TRAIL_LENGTH];
    int trailIndex;
};

void drawDude(GLfloat r, GLfloat g, GLfloat b, GLfloat rad, GLfloat posx, GLfloat posy) {
    glColor3f(r,g,b);
    glPointSize(rad);
    glBegin(GL_POINTS);
    glVertex2f(posx, posy);
    glEnd();
}

Acceleration gravity(Object object1, Object object2) {
    float G = 6.6743f * pow(10, -11);
    float r = sqrt(pow(object1.position.x - object2.position.x, 2) + pow(object1.position.y - object2.position.y, 2));
    float base = G * object2.mass / pow(r, 3);

    Acceleration acceleration;
    acceleration.x = base * (object2.position.x - object1.position.x);
    acceleration.y = base * (object2.position.y - object1.position.y);

    return acceleration;

};


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

    Object obj1;
    Object obj2;

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

    // Loop until the user closes the window
    while (!glfwWindowShouldClose(window)) {

        float last_time = (float)glfwGetTime();

        // Clear the screen to black
        glClearColor(0.0f, 0.0f, 0.0f, 1.0f);
        glClear(GL_COLOR_BUFFER_BIT);

        // Update the gravity acceleration of both objects
        obj1.acceleration = gravity(obj1, obj2);
        obj2.acceleration = gravity(obj2, obj1);

        // Get the time taken between frames
        float time = (float)glfwGetTime();
        float delta_time = time - last_time;

        obj1.velocity.x = obj1.velocity.x + obj1.acceleration.x * delta_time;
        obj1.velocity.y = obj1.velocity.y + obj1.acceleration.y * delta_time;
        obj2.velocity.x = obj2.velocity.x + obj2.acceleration.x * delta_time;
        obj2.velocity.y = obj2.velocity.y + obj2.acceleration.y * delta_time;

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



        // Swap front and back buffers
        glfwSwapBuffers(window);

        // Poll for and process events
        glfwPollEvents();

    }

    glfwTerminate();
    return 0;
}

