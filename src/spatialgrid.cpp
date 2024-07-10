#include "../include/spatialgrid.h"
#include <cmath>

SpatialGrid::SpatialGrid(double worldWidth, double worldHeight, double cellSize)
        : worldWidth(worldWidth), worldHeight(worldHeight), cellSize(cellSize) {
        gridWidth = static_cast<int>(std::ceil(worldWidth / cellSize));
        gridHeight = static_cast<int>(std::ceil(worldHeight / cellSize));
        grid.resize(gridWidth * gridHeight);
    }

    void SpatialGrid::clear() {
        for (auto& cell : grid) {
            cell.objects.clear();
        }
    }

    void SpatialGrid::insert(Object* obj) {
        // Calculate the range of cells this object could occupy
        int minX = static_cast<int>((obj->position.x - obj->radius + worldWidth/2) / cellSize);
        int maxX = static_cast<int>((obj->position.x + obj->radius + worldWidth/2) / cellSize);
        int minY = static_cast<int>((obj->position.y - obj->radius + worldHeight/2) / cellSize);
        int maxY = static_cast<int>((obj->position.y + obj->radius + worldHeight/2) / cellSize);

        // Clamp to grid boundaries
        minX = std::max(0, std::min(minX, gridWidth - 1));
        maxX = std::max(0, std::min(maxX, gridWidth - 1));
        minY = std::max(0, std::min(minY, gridHeight - 1));
        maxY = std::max(0, std::min(maxY, gridHeight - 1));

        // Insert the object into all cells it occupies
        for (int y = minY; y <= maxY; ++y) {
            for (int x = minX; x <= maxX; ++x) {
                grid[y * gridWidth + x].objects.push_back(obj);
            }
        }
    }

    std::vector<Object*> SpatialGrid::getNeighbors(const Object& obj) {
        std::vector<Object*> neighbors;
        std::unordered_set<Object*> uniqueNeighbors;  // To avoid duplicates

        // Calculate the range of cells to check
        int minX = static_cast<int>((obj.position.x - obj.radius - cellSize + worldWidth/2) / cellSize);
        int maxX = static_cast<int>((obj.position.x + obj.radius + cellSize + worldWidth/2) / cellSize);
        int minY = static_cast<int>((obj.position.y - obj.radius - cellSize + worldHeight/2) / cellSize);
        int maxY = static_cast<int>((obj.position.y + obj.radius + cellSize + worldHeight/2) / cellSize);

        // Clamp to grid boundaries
        minX = std::max(0, std::min(minX, gridWidth - 1));
        maxX = std::max(0, std::min(maxX, gridWidth - 1));
        minY = std::max(0, std::min(minY, gridHeight - 1));
        maxY = std::max(0, std::min(maxY, gridHeight - 1));

        for (int y = minY; y <= maxY; ++y) {
            for (int x = minX; x <= maxX; ++x) {
                for (Object* neighbor : grid[y * gridWidth + x].objects) {
                    if (neighbor != &obj && uniqueNeighbors.insert(neighbor).second) {
                        neighbors.push_back(neighbor);
                    }
                }
            }
        }

        return neighbors;
    }
