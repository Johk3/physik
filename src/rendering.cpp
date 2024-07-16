#include "../include/rendering.h"
#include "../include/constants.h"
#include "../include/linmath.h"
#include "../include/settings.h"

#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/type_ptr.hpp>

void drawSphere(const Vector3& position, double radius, const double r, const double g, const double b, const double alpha) {
    int stacks = 20;
    int slices = 20;

    glPushMatrix();
    glTranslatef(position.x, position.y, position.z);
    glColor4f(r, g, b, alpha);

    for (int i = 0; i <= stacks; ++i) {
        double lat0 = M_PI * (-0.5 + (double)(i - 1) / stacks);
        double z0 = sin(lat0);
        double zr0 = cos(lat0);

        double lat1 = M_PI * (-0.5 + (double)i / stacks);
        double z1 = sin(lat1);
        double zr1 = cos(lat1);

        glBegin(GL_QUAD_STRIP);
        for (int j = 0; j <= slices; ++j) {
            double lng = 2 * M_PI * (double)(j - 1) / slices;
            double x = cos(lng);
            double y = sin(lng);

            glNormal3f(x * zr0, y * zr0, z0);
            glVertex3f(x * zr0 * radius, y * zr0 * radius, z0 * radius);
            glNormal3f(x * zr1, y * zr1, z1);
            glVertex3f(x * zr1 * radius, y * zr1 * radius, z1 * radius);
        }
        glEnd();
    }

    glPopMatrix();
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
            // drawCircle(obj.trail[i].x, obj.trail[i].y, trailRadius, obj.r, obj.g, obj.b, alpha);
            drawSphere(obj.trail[i], trailRadius, obj.r, obj.g, obj.b, alpha);
        }
    };

    // Draw main object
    if (std::abs(obj.position.x) > PREVENT_DRAW_DISTANCE or std::abs(obj.position.y) > PREVENT_DRAW_DISTANCE) {return;};
    //drawCircle(obj.position.x, obj.position.y, obj.radius, obj.r, obj.g, obj.b, 1.0);
    drawSphere(obj.position, obj.radius, obj.r, obj.g, obj.b, 1.0);
}

// Renders the screen window based on objects on it
void render_screen(const std::vector<Object>& all_objects, GLFWwindow* window, const glm::mat4& view, const glm::mat4& projection) {
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    glMatrixMode(GL_PROJECTION);
    glLoadMatrixf(glm::value_ptr(projection));

    glMatrixMode(GL_MODELVIEW);
    glLoadMatrixf(glm::value_ptr(view));

    for (const auto& obj : all_objects) {
        drawObject(obj);
    }

    glfwSwapBuffers(window);
    glfwPollEvents();
};