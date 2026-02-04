/**
 * @file spatial_index.cpp
 * @brief Implementation of spatial indexing and neighbor search
 */

#include "spatial_index.h"
#include <cmath>
#include <iostream>

SpatialIndex::SpatialIndex(double distThresh)
    : distanceThreshold(distThresh)
{

}

double SpatialIndex::euclideanDist(const SpatialInstance& a, const SpatialInstance& b) {
    // Calculate Euclidean distance using the Pythagorean theorem
    double dx = a.x - b.x;
    double dy = a.y - b.y;
    return std::sqrt(dx * dx + dy * dy);
}

std::vector<std::pair<SpatialInstance, SpatialInstance>> SpatialIndex::findNeighborPair(const std::vector<SpatialInstance>& instances) {
    std::vector<std::pair<SpatialInstance, SpatialInstance>> neighborPairs;

    int N = instances.size();

    // Brute-force approach: compare all pairs of instances
    // Time complexity: O(n²) where n is the number of instances
    for (int i = 0; i < N; ++i) {
        for (int j = i + 1; j < N; ++j) {
            if (instances[i].type != instances[j].type) {
                // Calculate distance between instances i and j
                double dist = euclideanDist(instances[i], instances[j]);

                // If distance is within threshold, they are neighbors
                if (dist <= distanceThreshold) {
                    neighborPairs.push_back({ instances[i], instances[j] });
                }
            }
        }
    }
    
    return neighborPairs;
}