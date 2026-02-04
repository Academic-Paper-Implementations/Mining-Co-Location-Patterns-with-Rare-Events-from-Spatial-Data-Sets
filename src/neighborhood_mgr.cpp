/**
 * @file neighborhood_mgr.cpp
 * @brief Implementation of star neighborhood management
 */

#include "neighborhood_mgr.h"
#include "utils.h"
#include <algorithm>

void NeighborhoodMgr::buildFromPairs(const std::vector<std::pair<SpatialInstance, SpatialInstance>>& pairs) {
    // Build star neighborhoods from neighbor pairs
    // A star neighborhood has a center instance and all its neighbors
    
    for (const auto& pair : pairs) {
        const SpatialInstance& center = pair.first;
        const SpatialInstance& neighbor = pair.second;
        
        // Get or create the vector of star neighborhoods for this feature type
        auto& vec = starNeighborhoods[pair.first.type];

        // Search for existing star neighborhood with this center
        auto it = std::find_if(vec.begin(), vec.end(), [&](const StarNeighborhood& sn) {
            return sn.center->id == center.id;
        });

        if (it != vec.end()) {
            // Star neighborhood already exists, add this neighbor to it
            it->neighbors.push_back(&neighbor);
        } else {
            // Create new star neighborhood with this center
            StarNeighborhood newStarNeigh;
            newStarNeigh.center = &center;
            newStarNeigh.neighbors.push_back(&neighbor);
            vec.push_back(newStarNeigh);
        }
    }
    return;
}

const std::unordered_map<FeatureType, std::vector<StarNeighborhood>>& NeighborhoodMgr::getAllStarNeighborhoods() const {
    return starNeighborhoods;
}