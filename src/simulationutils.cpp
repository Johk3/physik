#include "../include/simulationutils.h"
#include <glm/glm.hpp>
#include <glm/geometric.hpp>
#include <iostream>
#include "../include/settings.h"
#include <fstream>

void calculate_gravity(Object& object1, const std::vector<Object>& objects, size_t start, size_t end) {
    glm::vec2 acc(0.0f);
    const glm::vec2 pos1 = glm::vec2(object1.position.x, object1.position.y);

    for (size_t i = start; i < end; ++i) {
        const auto& object2 = objects[i];
        if (&object1 == &object2 || object2.mass < Settings::EPSILON) continue;

        glm::vec2 pos2(object2.position.x, object2.position.y);
        glm::vec2 diff = pos2 - pos1;
        float r2 = glm::dot(diff, diff);

        if (r2 < 1e-12f) continue;

        float r = glm::sqrt(r2);
        if (r < 1e-6f) continue;

        float force = static_cast<float>(6.67430e-11 * object2.mass / (r2 * r));
        force = glm::min(force, (float)Settings::g_MAX_FORCE);

        acc += diff * force;
    }

    object1.acceleration = Vector2{acc.x, acc.y};
}

void updateObjectTrail(Object& obj) {
    if (obj.trail.empty()) {
        obj.trail.push_back(obj.position);
        return;
    }

    const glm::vec2 currentPos(obj.position.x, obj.position.y);
    const glm::vec2 lastPos(obj.trail.back().x, obj.trail.back().y);
    const glm::vec2 diff = currentPos - lastPos;
    float distance = glm::length(diff);

    if (distance >= obj.TRAIL_SPACING) {
        const int numNewDots = static_cast<int>(distance / obj.TRAIL_SPACING);
        const glm::vec2 direction = glm::normalize(diff);

        for (int i = 0; i < numNewDots; ++i) {
            const float t = (i + 1) * obj.TRAIL_SPACING;
            glm::vec2 newDotPosition = lastPos + direction * t;
            obj.trail.push_back(Vector2{newDotPosition.x, newDotPosition.y});

            if (obj.trail.size() > obj.MAX_TRAIL_LENGTH) {
                obj.trail.erase(obj.trail.begin());
            }
        }
    }
}

bool checkCollision(const Object& obj1, const Object& obj2) {
    const glm::vec2 pos1(obj1.position.x, obj1.position.y);
    const glm::vec2 pos2(obj2.position.x, obj2.position.y);
    return glm::distance(pos1, pos2) <= (obj1.radius + obj2.radius);
}

void handleCollision(Object& obj1, Object& obj2) {
    const glm::vec2 pos1(obj1.position.x, obj1.position.y);
    const glm::vec2 pos2(obj2.position.x, obj2.position.y);
    const glm::vec2 vel1(obj1.velocity.x, obj1.velocity.y);
    const glm::vec2 vel2(obj2.velocity.x, obj2.velocity.y);

    const glm::vec2 delta = pos2 - pos1;
    const float distance = glm::length(delta);

    // Early exit if objects are too far apart
    if (distance > obj1.radius + obj2.radius) return;

    const glm::vec2 normal = glm::normalize(delta);
    const glm::vec2 relativeVelocity = vel2 - vel1;
    const float velocityAlongNormal = glm::dot(relativeVelocity, normal);

    // Do not resolve if velocities are separating
    if (velocityAlongNormal > 0) return;

    const float restitution = Settings::g_bounceFactor;
    float impulseScalar = -(1.0f + restitution) * velocityAlongNormal;
    impulseScalar /= 1.0f / obj1.mass + 1.0f / obj2.mass;

    const glm::vec2 impulse = normal * impulseScalar;

    // Update velocities
    obj1.velocity = Vector2{
        vel1.x - (impulse.x / obj1.mass),
        vel1.y - (impulse.y / obj1.mass)
    };

    obj2.velocity = Vector2{
        vel2.x + (impulse.x / obj2.mass),
        vel2.y + (impulse.y / obj2.mass)
    };

    // Separate the objects
    const float overlap = (obj1.radius + obj2.radius) - distance;
    if (overlap > 0) {
        const float totalMass = obj1.mass + obj2.mass;
        const glm::vec2 separation = normal * overlap;

        const glm::vec2 pos1New = pos1 - separation * ((float)obj2.mass / totalMass);
        const glm::vec2 pos2New = pos2 + separation * ((float)obj1.mass / totalMass);

        obj1.position = Vector2{pos1New.x, pos1New.y};
        obj2.position = Vector2{pos2New.x, pos2New.y};
    }
}

void updateSimulation(std::vector<Object>& allObjects, SpatialGrid& grid, ThreadPool& pool, double delta_time) {
    const size_t numObjects = allObjects.size();
    const size_t numThreads = pool.getThreadCount();
    const size_t batchSize = std::max(size_t(1), numObjects / (numThreads * 4));
    const float dt = static_cast<float>(delta_time);

    // Step 1: Gravity Calculation
    std::vector<std::future<void>> futures;
    for (size_t i = 0; i < numObjects; i += batchSize) {
        size_t end = std::min(i + batchSize, numObjects);
        futures.push_back(pool.enqueue([&allObjects, i, end]() {
            for (size_t j = i; j < end; ++j) {
                calculate_gravity(allObjects[j], allObjects, 0, allObjects.size());
            }
        }));
    }

    for (auto& future : futures) {
        future.get();
    }
    futures.clear();

    // Step 2: Update positions and trails
    for (size_t i = 0; i < numObjects; i += batchSize) {
        size_t end = std::min(i + batchSize, numObjects);
        futures.push_back(pool.enqueue([&allObjects, i, end, dt]() {
            for (size_t j = i; j < end; ++j) {
                Object& obj = allObjects[j];
                glm::vec2 vel(obj.velocity.x, obj.velocity.y);
                glm::vec2 acc(obj.acceleration.x, obj.acceleration.y);
                glm::vec2 pos(obj.position.x, obj.position.y);

                vel += acc * dt;
                pos += vel * dt;

                obj.velocity = Vector2{vel.x, vel.y};
                obj.position = Vector2{pos.x, pos.y};
                updateObjectTrail(obj);
            }
        }));
    }

    for (auto& future : futures) {
        future.get();
    }
    futures.clear();

    // Step 3: Rebuild spatial grid
    grid.clear();
    for (auto& obj : allObjects) {
        grid.insert(&obj);
    }

    // Step 4: Handle collisions
    for (size_t i = 0; i < numObjects; i += batchSize) {
        size_t end = std::min(i + batchSize, numObjects);
        futures.push_back(pool.enqueue([&allObjects, &grid, i, end]() {
            for (size_t j = i; j < end; ++j) {
                auto neighbors = grid.getNeighbors(allObjects[j]);
                for (auto* neighbor : neighbors) {
                    if (checkCollision(allObjects[j], *neighbor)) {
                        handleCollision(allObjects[j], *neighbor);
                    }
                }
            }
        }));
    }

    for (auto& future : futures) {
        future.get();
    }
}

std::vector<Object> get_objects() {
    std::vector<Object> allObjects;
    allObjects.reserve(626); // Pre-allocate space for all objects (25*25 + 1)

    // Create grid of objects
    for (int i = 0; i < 25; i++) {
        for (int j = 0; j < 25; j++) {
            const float r = glm::mix(0.1f, 1.0f, float(i) / 24.0f);
            const float g = glm::mix(0.1f, 1.0f, float(j) / 24.0f);
            const float b = glm::mix(0.1f, 1.0f, 1.0f - float(i + j) / 48.0f);

            allObjects.emplace_back(
                Vector2{-0.5f + i * 0.04f, j * 0.04f},
                Vector2{0.2f, 0.0f},
                1e4,
                0.1,
                r, g, b
            );
        }
    }

    // Add central object
    allObjects.emplace_back(
        Vector2{0.02f, 0.0001f},
        Vector2{0.2f, 0.0f},
        5e14,
        5e6,
        1.0f, 1.0f, 1.0f
    );

    return allObjects;
}