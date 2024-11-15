#ifndef SIMULATIONUTILS_H
#define SIMULATIONUTILS_H

#include "object.h"
#include "spatialgrid.h"
#include "ThreadPool.h"
#include <vector>

void calculate_gravity(Object& object1, const std::vector<Object>& objects, size_t start, size_t end);

void updateObjectTrail(Object& obj);

bool checkCollision(const Object& obj1, const Object& obj2);

void handleCollision(Object& obj1, Object& obj2);

void updateSimulation(std::vector<Object>& allObjects, SpatialGrid& grid, ThreadPool& pool, double delta_time);
std::vector<Object> get_objects();

#endif // SIMULATION_UTILS_H