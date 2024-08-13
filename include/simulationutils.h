#ifndef SIMULATIONUTILS_H
#define SIMULATIONUTILS_H

#include "object.h"
#include "spatialgrid.h"
#include "ThreadPool.h"
#include <vector>
#ifdef USE_GPU
#include <CL/cl2.hpp>
#endif
// Simulation-related functions
void updateObjectTrail(Object& obj);
bool checkCollision(const Object& obj1, const Object& obj2);
void handleCollision(Object& obj1, Object& obj2);

#ifdef USE_AVX2
void calculate_gravity_simd(Object& object1, const std::vector<Object>& objects, size_t start, size_t end);
#endif

void calculate_gravity_normal(Object& object1, const std::vector<Object>& objects, size_t start, size_t end);

// Updated to include ThreadPool
void updateSimulation(std::vector<Object>& allObjects, SpatialGrid& grid, ThreadPool& pool, double delta_time);

std::vector<Object> get_objects();
#ifdef USE_GPU
// OpenCL related functions
void initializeOpenCL();
void cleanupOpenCL();

// Global OpenCL variables (you might want to wrap these in a namespace or class)
extern cl::Context context;
extern cl::CommandQueue queue;
extern cl::Kernel gravityKernel;
extern cl::Kernel updatePositionKernel;
extern cl::Kernel collisionKernel;
extern cl::Program program;
#endif

#endif // SIMULATION_UTILS_H