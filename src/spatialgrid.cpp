#include "../include/spatialgrid.h"
#include <cmath>
#include <algorithm>

namespace {
    // Constants for optimization
    constexpr size_t INITIAL_CELL_CAPACITY = 8;
    constexpr size_t NEIGHBOR_RESERVE_SIZE = 32;

    // Helper functions
    inline int fastFloor(double x) {
        int i = static_cast<int>(x);
        return i - (x < i);
    }

    inline int getCellIndex(int x, int y, int gridWidth) {
        return y * gridWidth + x;
    }
}

SpatialGrid::SpatialGrid(double worldWidth, double worldHeight, double cellSize)
    : worldWidth(worldWidth)
    , worldHeight(worldHeight)
    , cellSize(cellSize)
    , gridWidth(static_cast<int>(std::ceil(worldWidth / cellSize)))
    , gridHeight(static_cast<int>(std::ceil(worldHeight / cellSize))) {

    // Pre-allocate grid with estimated size
    grid.resize(gridWidth * gridHeight);

    // Pre-reserve space in each cell to reduce reallocations
    for (auto& cell : grid) {
        cell.objects.reserve(INITIAL_CELL_CAPACITY);
    }
}

void SpatialGrid::clear() {
    // Clear cells while preserving capacity
    for (auto& cell : grid) {
        cell.objects.clear();
    }
}

void SpatialGrid::insert(Object* obj) {
    // Fast coordinate conversion with offset
    const double halfWidth = worldWidth * 0.5;
    const double halfHeight = worldHeight * 0.5;

    // Calculate grid coordinates
    int minX = fastFloor((obj->position.x - obj->radius + halfWidth) / cellSize);
    int maxX = fastFloor((obj->position.x + obj->radius + halfWidth) / cellSize);
    int minY = fastFloor((obj->position.y - obj->radius + halfHeight) / cellSize);
    int maxY = fastFloor((obj->position.y + obj->radius + halfHeight) / cellSize);

    // Clamp to grid boundaries
    minX = std::clamp(minX, 0, gridWidth - 1);
    maxX = std::clamp(maxX, 0, gridWidth - 1);
    minY = std::clamp(minY, 0, gridHeight - 1);
    maxY = std::clamp(maxY, 0, gridHeight - 1);

    // Insert using linear indexing for better cache coherency
    for (int y = minY; y <= maxY; ++y) {
        const int rowOffset = y * gridWidth;
        for (int x = minX; x <= maxX; ++x) {
            grid[rowOffset + x].objects.push_back(obj);
        }
    }
}

std::vector<Object*> SpatialGrid::getNeighbors(const Object& obj) {
    std::vector<Object*> neighbors;
    neighbors.reserve(NEIGHBOR_RESERVE_SIZE);

    // Temporary vector for tracking unique objects (better cache locality than unordered_set)
    std::vector<Object*> uniqueObjects;
    uniqueObjects.reserve(NEIGHBOR_RESERVE_SIZE);

    // Calculate grid range with a small overlap for safety
    const double halfWidth = worldWidth * 0.5;
    const double halfHeight = worldHeight * 0.5;

    int minX = fastFloor((obj.position.x - obj.radius - cellSize + halfWidth) / cellSize);
    int maxX = fastFloor((obj.position.x + obj.radius + cellSize + halfWidth) / cellSize);
    int minY = fastFloor((obj.position.y - obj.radius - cellSize + halfHeight) / cellSize);
    int maxY = fastFloor((obj.position.y + obj.radius + cellSize + halfHeight) / cellSize);

    // Clamp to grid boundaries
    minX = std::clamp(minX, 0, gridWidth - 1);
    maxX = std::clamp(maxX, 0, gridWidth - 1);
    minY = std::clamp(minY, 0, gridHeight - 1);
    maxY = std::clamp(maxY, 0, gridHeight - 1);

    // Gather neighbors using linear indexing
    for (int y = minY; y <= maxY; ++y) {
        const int rowOffset = y * gridWidth;
        for (int x = minX; x <= maxX; ++x) {
            const auto& cellObjects = grid[rowOffset + x].objects;
            for (Object* neighbor : cellObjects) {
                if (neighbor != &obj &&
                    std::find(uniqueObjects.begin(), uniqueObjects.end(), neighbor) == uniqueObjects.end()) {
                    uniqueObjects.push_back(neighbor);
                    neighbors.push_back(neighbor);
                }
            }
        }
    }

    return neighbors;
}