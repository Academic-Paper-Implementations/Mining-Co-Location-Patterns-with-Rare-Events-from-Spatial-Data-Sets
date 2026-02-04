/**
 * @file spatial_index.h
 * @brief Spatial indexing and neighbor search functionality
 */

#pragma once
#include "types.h"
#include <vector>

/**
 * @brief SpatialIndex class for managing spatial indexing and neighbor searches
 * 
 * Provides functionality to find neighboring spatial instances within a distance threshold.
 * Currently uses a brute-force O(n²) approach for neighbor pair detection.
 */
class SpatialIndex {
private:
    double distanceThreshold;  ///< Distance threshold for neighbor determination

    /**
     * @brief Calculate Euclidean distance between two spatial instances
     * 
     * @param a First spatial instance
     * @param b Second spatial instance
     * @return double Euclidean distance between a and b
     */
    double euclideanDist(const SpatialInstance& a, const SpatialInstance& b);
    
public:
    /**
     * @brief Constructor to initialize SpatialIndex with a distance threshold
     * 
     * @param distThresh Maximum distance for two instances to be considered neighbors
     */
    SpatialIndex(double distThresh);

    /**
     * @brief Find all neighbor pairs within the distance threshold
     * 
     * Uses brute-force comparison (O(n²)) to find all pairs of instances
     * that are within the specified distance threshold.
     * 
     * @param instances Vector of all spatial instances to search
     * @return std::vector<std::pair<SpatialInstance, SpatialInstance>> Vector of neighbor pairs
     * @note Time complexity: O(n²) where n is the number of instances
     */
    std::vector<std::pair<SpatialInstance, SpatialInstance>> findNeighborPair(const std::vector<SpatialInstance>& instances);
};