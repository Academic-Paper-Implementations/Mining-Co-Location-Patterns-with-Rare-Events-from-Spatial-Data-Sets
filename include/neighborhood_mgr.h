/**
 * @file neighborhood_mgr.h
 * @brief Star neighborhood management for spatial instances
 */

#pragma once
#include "types.h"
#include <unordered_map>
#include <vector>

/**
 * @brief NeighborhoodMgr class for managing star neighborhoods of spatial instances
 * 
 * Organizes spatial instances into star neighborhoods, where each star consists of
 * a center instance and all its neighbors within the distance threshold.
 */
class NeighborhoodMgr {
private:
    /// Map from feature type to all star neighborhoods of that type
    std::unordered_map<FeatureType, std::vector<StarNeighborhood>> starNeighborhoods;

public:
    /**
     * @brief Build star neighborhoods from neighbor pairs
     * 
     * Constructs star neighborhoods by grouping neighbor pairs. For each instance,
     * creates a star with that instance as center and all its neighbors.
     * 
     * @param pairs Vector of neighbor pairs found by spatial indexing
     */
    void buildFromPairs(const std::vector<std::pair<SpatialInstance, SpatialInstance>>& pairs);

    
    /**
     * @brief Get all star neighborhoods organized by feature type
     * 
     * @return const std::unordered_map<FeatureType, std::vector<StarNeighborhood>>& 
     *         Map from feature type to vector of star neighborhoods
     */
    const std::unordered_map<FeatureType, std::vector<StarNeighborhood>>& getAllStarNeighborhoods() const;
};