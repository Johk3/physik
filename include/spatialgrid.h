#ifndef SPATIALGRID_H
#define SPATIALGRID_H

#include "object.h"
#include <vector>
#include <unordered_set>

// Optimized spatial partitioning grid for efficient neighbor queries
class SpatialGrid {
private:
    struct Cell {
        std::vector<Object*> objects;
    };

    std::vector<Cell> grid;
    int gridWidth, gridHeight;
    double cellSize;
    double worldWidth, worldHeight;

public:
    SpatialGrid(double worldWidth, double worldHeight, double cellSize);

    void clear();
    void insert(Object* obj);
    std::vector<Object*> getNeighbors(const Object& obj);
};

#endif // SPATIAL_GRID_H