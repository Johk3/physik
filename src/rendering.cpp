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

    // Set material properties
    GLfloat mat_ambient[] = { static_cast<GLfloat>(r) * 0.3f, static_cast<GLfloat>(g) * 0.3f, static_cast<GLfloat>(b) * 0.3f, static_cast<GLfloat>(alpha) };
    GLfloat mat_diffuse[] = { static_cast<GLfloat>(r), static_cast<GLfloat>(g), static_cast<GLfloat>(b), static_cast<GLfloat>(alpha) };
    GLfloat mat_specular[] = { 0.5f, 0.5f, 0.5f, 1.0f };
    GLfloat mat_shininess[] = { 50.0f };

    glMaterialfv(GL_FRONT, GL_AMBIENT, mat_ambient);
    glMaterialfv(GL_FRONT, GL_DIFFUSE, mat_diffuse);
    glMaterialfv(GL_FRONT, GL_SPECULAR, mat_specular);
    glMaterialfv(GL_FRONT, GL_SHININESS, mat_shininess);

    // Set the color
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
void drawSlimArrow(const Vector3& start, const Vector3& direction, float length, float r, float g, float b) {
    Vector3 end;
    end.x = start.x + direction.x * length;
    end.y = start.y + direction.y * length;
    end.z = start.z + direction.z * length;

    // Draw the arrow body
    glColor3f(r, g, b);
    glLineWidth(2.0f);  // Make the line slimmer
    glBegin(GL_LINES);
    glVertex3f(start.x, start.y, start.z);
    glVertex3f(end.x, end.y, end.z);
    glEnd();

    // Draw the arrow head
    float arrowHeadLength = length * 0.2f; // 20% of the arrow length
    float arrowHeadWidth = arrowHeadLength * 0.5f;

    Vector3 up, right;

    // Calculate an arbitrary 'up' vector
    if (std::abs(direction.x) < std::abs(direction.y) && std::abs(direction.x) < std::abs(direction.z)) {
        up.x = 1; up.y = 0; up.z = 0;
    } else if (std::abs(direction.y) < std::abs(direction.z)) {
        up.x = 0; up.y = 1; up.z = 0;
    } else {
        up.x = 0; up.y = 0; up.z = 1;
    }

    // Calculate right vector
    right.x = direction.y * up.z - direction.z * up.y;
    right.y = direction.z * up.x - direction.x * up.z;
    right.z = direction.x * up.y - direction.y * up.x;
    float rightLength = std::sqrt(right.x * right.x + right.y * right.y + right.z * right.z);
    right.x /= rightLength; right.y /= rightLength; right.z /= rightLength;

    // Recalculate up vector
    up.x = right.y * direction.z - right.z * direction.y;
    up.y = right.z * direction.x - right.x * direction.z;
    up.z = right.x * direction.y - right.y * direction.x;

    glBegin(GL_TRIANGLES);
    // Arrow head in XY plane
    glVertex3f(end.x, end.y, end.z);
    glVertex3f(end.x - direction.x * arrowHeadLength + right.x * arrowHeadWidth,
               end.y - direction.y * arrowHeadLength + right.y * arrowHeadWidth,
               end.z - direction.z * arrowHeadLength + right.z * arrowHeadWidth);
    glVertex3f(end.x - direction.x * arrowHeadLength - right.x * arrowHeadWidth,
               end.y - direction.y * arrowHeadLength - right.y * arrowHeadWidth,
               end.z - direction.z * arrowHeadLength - right.z * arrowHeadWidth);

    // Arrow head in XZ plane
    glVertex3f(end.x, end.y, end.z);
    glVertex3f(end.x - direction.x * arrowHeadLength + up.x * arrowHeadWidth,
               end.y - direction.y * arrowHeadLength + up.y * arrowHeadWidth,
               end.z - direction.z * arrowHeadLength + up.z * arrowHeadWidth);
    glVertex3f(end.x - direction.x * arrowHeadLength - up.x * arrowHeadWidth,
               end.y - direction.y * arrowHeadLength - up.y * arrowHeadWidth,
               end.z - direction.z * arrowHeadLength - up.z * arrowHeadWidth);
    glEnd();

    glLineWidth(1.0f);  // Reset line width
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
    // Draw trail (if enabled)
    if (Settings::g_enableTrail) {
        for (size_t i = 0; i < obj.trail.size(); ++i) {
            if (std::abs(obj.trail[i].x) > PREVENT_DRAW_DISTANCE or std::abs(obj.trail[i].y) > PREVENT_DRAW_DISTANCE) continue;
            double alpha = static_cast<double>(i) / obj.trail.size();
            double trailRadius = obj.radius * TRAIL_SCALE * alpha;
            drawSphere(obj.trail[i], trailRadius, obj.r, obj.g, obj.b, alpha);
        }
    }
    // Draw the sphere
    if (std::abs(obj.position.x) > PREVENT_DRAW_DISTANCE or std::abs(obj.position.y) > PREVENT_DRAW_DISTANCE) {return;};
    drawSphere(obj.position, obj.radius, obj.r, obj.g, obj.b, 1.0);
    if (Settings::g_drawArrow){
    // Draw direction arrow
    float velocityLength = std::sqrt(obj.velocity.x * obj.velocity.x + obj.velocity.y * obj.velocity.y + obj.velocity.z * obj.velocity.z);
    float minVelocityThreshold = obj.radius; // Adjust this value as needed

    if (velocityLength > minVelocityThreshold) {
        Vector3 direction;
        direction.x = obj.velocity.x / velocityLength;
        direction.y = obj.velocity.y / velocityLength;
        direction.z = obj.velocity.z / velocityLength;

        float accelerationLength = std::sqrt(obj.acceleration.x * obj.acceleration.x +
                                             obj.acceleration.y * obj.acceleration.y +
                                             obj.acceleration.z * obj.acceleration.z);

        // Adjust these scale factors as needed
        float maxArrowLength = obj.radius*1.5f; // Maximum arrow length in simulation units
        float minArrowLength = obj.radius; // Minimum arrow length relative to object size
        float arrowLength = std::min((float)obj.radius*1.5f*(1-1/accelerationLength), maxArrowLength);
        arrowLength = std::max(arrowLength, minArrowLength);

        // Use a color that contrasts with both light and dark objects
        float arrowR = 1.0f - obj.r;
        float arrowG = 1.0f - obj.g;
        float arrowB = 1.0f - obj.b;

        // Start the arrow from the edge of the sphere in the direction of motion
        Vector3 arrowStart;
        arrowStart.x = obj.position.x + direction.x * obj.radius;
        arrowStart.y = obj.position.y + direction.y * obj.radius;
        arrowStart.z = obj.position.z + direction.z * obj.radius;

        drawSlimArrow(arrowStart, direction, arrowLength, arrowR, arrowG, arrowB);
    }
    }

}
// Renders the screen window based on objects on it
void render_screen(const std::vector<Object>& all_objects, GLFWwindow* window, const glm::mat4& view, const glm::mat4& projection) {
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    glMatrixMode(GL_PROJECTION);
    glLoadMatrixf(glm::value_ptr(projection));

    glMatrixMode(GL_MODELVIEW);
    glLoadMatrixf(glm::value_ptr(view));

    // Set up a basic light
    GLfloat light_position[] = { 5.0f, 5.0f, 5.0f, 1.0f };
    GLfloat light_ambient[] = { 0.2f, 0.2f, 0.2f, 1.0f };
    GLfloat light_diffuse[] = { 0.8f, 0.8f, 0.8f, 1.0f };
    GLfloat light_specular[] = { 0.5f, 0.5f, 0.5f, 1.0f };

    glLightfv(GL_LIGHT0, GL_POSITION, light_position);
    glLightfv(GL_LIGHT0, GL_AMBIENT, light_ambient);
    glLightfv(GL_LIGHT0, GL_DIFFUSE, light_diffuse);
    glLightfv(GL_LIGHT0, GL_SPECULAR, light_specular);

    for (const auto& obj : all_objects) {
        drawObject(obj);
    }

    glfwSwapBuffers(window);
    glfwPollEvents();
}