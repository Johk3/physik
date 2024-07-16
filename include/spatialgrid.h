#ifndef SPATIALGRID_H
#define SPATIALGRID_H
#include <vector>
#include <unordered_map>
#include "../include/object.h"

class SpatialGrid {
private:
    float cellSize;
    std::unordered_map<size_t, std::vector<Object*>> grid;

    // Hash function for 3D coordinates
    size_t hash(int x, int y, int z) const;

public:
    // Constructor
    SpatialGrid(float cellSize);

    // Clear the grid
    void clear();

    // Insert an object into the grid
    void insert(Object* obj);

    // Get neighbors of an object
    std::vector<Object*> getNeighbors(const Object& obj);
};

#endif // SPATIAL_GRID_H