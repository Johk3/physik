#include <GLFW/glfw3.h>
#include "linmath.h"
#include <math.h>
#define TRAIL_LENGTH 20

struct Position {
    float x;
    float y;
};

Position trail[TRAIL_LENGTH];
int trailIndex = 0;

void drawDude(GLfloat r, GLfloat g, GLfloat b, GLfloat rad, GLfloat posx, GLfloat posy) {
    glColor3f(r,g,b);
    glPointSize(rad);
    glBegin(GL_POINTS);
    glVertex2f(posx, posy);
    glEnd();
}

Position circularMotion(float radius) {
    // Calculate dot position
    float time = (float)glfwGetTime();
    Position position;
    position.x = radius * cosf(time);
    position.y = radius * sinf(time*0.5f);

    return position;
}

void drawTrail() {
    for (int i = 0; i < TRAIL_LENGTH; i++) {
        int index = (trailIndex - i - 1 + TRAIL_LENGTH) % TRAIL_LENGTH;
        float alpha = 1.0f - (float)i / TRAIL_LENGTH;

        mat4x4 mvp;
        mat4x4_identity(mvp);
        mat4x4_translate(mvp, trail[index].x, trail[index].y, 0.0f);

        glMatrixMode(GL_PROJECTION);
        glLoadIdentity();
        glMultMatrixf((const GLfloat*)mvp);

        drawDude(1.0f, 1.0f, 1.0f, 5.0f * alpha, 0.0f, 0.0f);
    }
}

int main() {
    // Initialize GLFW
    if (!glfwInit()) {
        return -1;
    }

    // Create a windowed mode window and its OpenGL context
    GLFWwindow* window = glfwCreateWindow(640, 480, "Moving White Dot", NULL, NULL);
    if (!window) {
        glfwTerminate();
        return -1;
    }

    // Make the window's context current
    glfwMakeContextCurrent(window);

    // Loop until the user closes the window
    while (!glfwWindowShouldClose(window)) {
        // Clear the screen to black
        glClearColor(0.0f, 0.0f, 0.0f, 1.0f);
        glClear(GL_COLOR_BUFFER_BIT);

        Position currentPos = circularMotion(0.5f);

        // Store current position in trail
        trail[trailIndex] = currentPos;
        trailIndex = (trailIndex + 1) % TRAIL_LENGTH;

        // Draw trail
        drawTrail();

        // Draw current position (main dot)
        mat4x4 mvp;
        mat4x4_identity(mvp);
        mat4x4_translate(mvp, currentPos.x, currentPos.y, 0.0f);

        // apply the mvp matrix
        glMatrixMode(GL_PROJECTION);
        glLoadIdentity();
        glMultMatrixf((const GLfloat*)mvp);
        // Draw a white point
        drawDude(1.0f,1.0f,1.0f, 10.0f, 0.0f, 0.0f);

        // Swap front and back buffers
        glfwSwapBuffers(window);

        // Poll for and process events
        glfwPollEvents();
    }

    glfwTerminate();
    return 0;
}

