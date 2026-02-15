/**
 * @file neighborhood_mgr.cpp
 * @brief Implementation of star neighborhood management
 */

#include "neighborhood_mgr.h"
#include "utils.h"
#include <algorithm>

void NeighborhoodMgr::buildFromPairs(const std::vector<std::pair<const SpatialInstance*, const SpatialInstance*>>& pairs) {
    // Build star neighborhoods from neighbor pairs
    // A star neighborhood has a center instance and all its neighbors
    
    for (const auto& pair : pairs) {
        const SpatialInstance* center = pair.first;
        const SpatialInstance* neighbor = pair.second;
        
        // Get or create the vector of star neighborhoods for this feature type
        auto& vec = starNeighborhoods[center->type];
        // Search for existing star neighborhood with this center
        auto it = std::find_if(vec.begin(), vec.end(), [&](const StarNeighborhood& sn) {
            return sn.center == center;
        });

        if (it != vec.end()) {
            // Star neighborhood already exists, add this neighbor to it
            it->neighbors.push_back(neighbor);
        } else {
            // Create new star neighborhood with this center
            StarNeighborhood newStarNeigh;
            newStarNeigh.center = center;
            newStarNeigh.neighbors.push_back(neighbor);
            vec.push_back(newStarNeigh);
        }
    }
    return;
}

bool NeighborhoodMgr::areNeighbors(const SpatialInstance* instA, const SpatialInstance* instB) const {
    if (!instA || !instB || instA == instB) return false;

    // Check if instB is in instA's star neighborhood
    auto itA = starNeighborhoods.find(instA->type);
    if (itA != starNeighborhoods.end()) {
        for (const auto& star : itA->second) {
            if (star.center == instA) {
                for (const auto* neighbor : star.neighbors) {
                    if (neighbor == instB) return true;
                }
                break;
            }
        }
    }

    // Check if instA is in instB's star neighborhood
    auto itB = starNeighborhoods.find(instB->type);
    if (itB != starNeighborhoods.end()) {
        for (const auto& star : itB->second) {
            if (star.center == instB) {
                for (const auto* neighbor : star.neighbors) {
                    if (neighbor == instA) return true;
                }
                break;
            }
        }
    }

    return false;
}

const std::unordered_map<FeatureType, std::vector<StarNeighborhood>>& NeighborhoodMgr::getAllStarNeighborhoods() const {
    return starNeighborhoods;
}