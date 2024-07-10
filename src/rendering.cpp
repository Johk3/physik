#include "../include/rendering.h"
#include "../include/constants.h"
#include "../include/linmath.h"
#include "../include/settings.h"

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

    int circle_segments = Settings::g_circleSegments;

    // LOD optimization
    if (Settings::g_enableLOD) {
        circle_segments = std::min(int(round(150 * radius + 5)), Settings::g_circleSegments);
    };

    glColor4d(r, g, b, alpha);
    glBegin(GL_TRIANGLE_FAN);
    glVertex2d(x, y);
    for (int i = 0; i <= circle_segments; i++) {
        const double theta = 2.0 * M_PI * static_cast<double>(i) / static_cast<double>(circle_segments);
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
    if (Settings::g_enableTrail) {
        for (size_t i = 0; i < obj.trail.size(); ++i) {
            if (std::abs(obj.trail[i].x) > PREVENT_DRAW_DISTANCE or std::abs(obj.trail[i].y) > PREVENT_DRAW_DISTANCE) continue;
            double alpha = static_cast<double>(i) / obj.trail.size();
            double trailRadius = obj.radius * TRAIL_SCALE * alpha;
            drawCircle(obj.trail[i].x, obj.trail[i].y, trailRadius, obj.r, obj.g, obj.b, alpha);
        }
    };

    // Draw main object
    if (std::abs(obj.position.x) > PREVENT_DRAW_DISTANCE or std::abs(obj.position.y) > PREVENT_DRAW_DISTANCE) {return;};
    drawCircle(obj.position.x, obj.position.y, obj.radius, obj.r, obj.g, obj.b, 1.0);
}

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