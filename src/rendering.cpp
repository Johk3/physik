#include "../include/rendering.h"
#include "../include/constants.h"
#include "../include/linmath.h"
#include "../include/settings.h"


// Camera state
struct Camera {
    double x = 0.0;
    double y = 0.0;
    double zoom = 1.0;
} camera;

// Mouse state
struct MouseState {
    bool leftButtonDown = false;
    double lastX = 0.0;
    double lastY = 0.0;
} mouseState;

// Function to transform world coordinates to screen coordinates
void worldToScreen(double worldX, double worldY, double& screenX, double& screenY) {
    screenX = (worldX - camera.x) * camera.zoom;
    screenY = (worldY - camera.y) * camera.zoom;
}

// -- Objects which can be drawn
void drawSquare(const double x, const double y, const double r, const double g, const double b, const double alpha, const double size) {
    double screenX, screenY;
    worldToScreen(x, y, screenX, screenY);

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
    double screenX, screenY;
    worldToScreen(x, y, screenX, screenY);

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
            double screenX, screenY;
            worldToScreen(obj.trail[i].x, obj.trail[i].y, screenX, screenY);
            if (std::abs(screenX) > PREVENT_DRAW_DISTANCE || std::abs(screenY) > PREVENT_DRAW_DISTANCE) continue;
            double alpha = static_cast<double>(i) / obj.trail.size();
            double trailRadius = obj.radius * TRAIL_SCALE * alpha;
            drawCircle(obj.trail[i].x, obj.trail[i].y, trailRadius, obj.r, obj.g, obj.b, alpha);
        }
    }

    // Draw main object
    double screenX, screenY;
    worldToScreen(obj.position.x, obj.position.y, screenX, screenY);
    if (std::abs(screenX) > PREVENT_DRAW_DISTANCE || std::abs(screenY) > PREVENT_DRAW_DISTANCE) return;
    drawCircle(obj.position.x, obj.position.y, obj.radius, obj.r, obj.g, obj.b, 1.0);
}

// Renders the screen window based on objects on it
void render_screen(const std::vector<Object>& all_objects, GLFWwindow* window) {
    // Clear screen
    glClearColor(0.0f, 0.0f, 0.0f, 1.0f);
    glClear(GL_COLOR_BUFFER_BIT);

    // Set up the view
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    glOrtho(-1.0 / camera.zoom, 1.0 / camera.zoom, -1.0 / camera.zoom, 1.0 / camera.zoom, -1.0, 1.0);

    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
    glTranslatef(-camera.x, -camera.y, 0.0f);

    // Draw all objects
    for (const auto& obj : all_objects) {
        drawObject(obj);
    }

    // Swap the back buffer with the front
    glfwSwapBuffers(window);
};

// Mouse button callback
void mouse_button_callback(GLFWwindow* window, int button, int action, int mods) {
    if (button == GLFW_MOUSE_BUTTON_LEFT) {
        if (action == GLFW_PRESS) {
            mouseState.leftButtonDown = true;
            glfwGetCursorPos(window, &mouseState.lastX, &mouseState.lastY);
        } else if (action == GLFW_RELEASE) {
            mouseState.leftButtonDown = false;
        }
    }
}

// Mouse movement callback
void mouse_callback(GLFWwindow* window, double xpos, double ypos) {
    if (mouseState.leftButtonDown) {
        double dx = xpos - mouseState.lastX;
        double dy = ypos - mouseState.lastY;
        camera.x -= dx / (500.0 * camera.zoom);
        camera.y += dy / (500.0 * camera.zoom);
        mouseState.lastX = xpos;
        mouseState.lastY = ypos;
    }
}

// Scroll callback
void scroll_callback(GLFWwindow* window, double xoffset, double yoffset) {
    camera.zoom *= (1.0 + yoffset * 0.1);
    if (camera.zoom < 0.1) camera.zoom = 0.1;
    if (camera.zoom > 10.0) camera.zoom = 10.0;
}

// Function to handle keyboard input
void process_input(GLFWwindow* window) {
    float cameraSpeed = 0.05f / camera.zoom;
    if (glfwGetKey(window, GLFW_KEY_W) == GLFW_PRESS)
        camera.y += cameraSpeed;
    if (glfwGetKey(window, GLFW_KEY_S) == GLFW_PRESS)
        camera.y -= cameraSpeed;
    if (glfwGetKey(window, GLFW_KEY_A) == GLFW_PRESS)
        camera.x -= cameraSpeed;
    if (glfwGetKey(window, GLFW_KEY_D) == GLFW_PRESS)
        camera.x += cameraSpeed;
}

// Function to set up input callbacks
void setup_input_callbacks(GLFWwindow* window) {
    glfwSetMouseButtonCallback(window, mouse_button_callback);
    glfwSetCursorPosCallback(window, mouse_callback);
    glfwSetScrollCallback(window, scroll_callback);
}
