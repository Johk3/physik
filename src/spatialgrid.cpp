#pragma once
#include "../include/spatialgrid.h"
#include <cmath>

// Constructor
SpatialGrid::SpatialGrid(float cellSize) : cellSize(cellSize) {}

// Hash function for 3D coordinates
size_t SpatialGrid::hash(int x, int y, int z) const {
    return x * 73856093 ^ y * 19349663 ^ z * 83492791;
}

// Clear the grid
void SpatialGrid::clear() {
    grid.clear();
}

// Insert an object into the grid
void SpatialGrid::insert(Object* obj) {
    int x = static_cast<int>(std::floor(obj->position.x / cellSize));
    int y = static_cast<int>(std::floor(obj->position.y / cellSize));
    int z = static_cast<int>(std::floor(obj->position.z / cellSize));
    size_t key = hash(x, y, z);
    grid[key].push_back(obj);
}

// Get neighbors of an object
std::vector<Object*> SpatialGrid::getNeighbors(const Object& obj) {
    std::vector<Object*> neighbors;
    int x = static_cast<int>(std::floor(obj.position.x / cellSize));
    int y = static_cast<int>(std::floor(obj.position.y / cellSize));
    int z = static_cast<int>(std::floor(obj.position.z / cellSize));

    // Check neighboring cells (including the current cell)
    for (int dx = -1; dx <= 1; ++dx) {
        for (int dy = -1; dy <= 1; ++dy) {
            for (int dz = -1; dz <= 1; ++dz) {
                size_t key = hash(x + dx, y + dy, z + dz);
                auto it = grid.find(key);
                if (it != grid.end()) {
                    neighbors.insert(neighbors.end(), it->second.begin(), it->second.end());
                }
            }
        }
    }
    return neighbors;
}