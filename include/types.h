/**
 * @file types.h
 * @brief Core data type definitions for spatial colocation pattern mining
 * 
 * This file defines the fundamental data structures and type aliases used throughout
 * the joinless colocation mining algorithm implementation.
 */

#pragma once
#include <string>
#include <unordered_map>
#include <vector>

// ============================================================================
// Type Aliases
// ============================================================================

/** @brief Type alias for feature types (e.g., "Restaurant", "Hotel") */
using FeatureType = std::string;

/** @brief Type alias for instance identifiers (e.g., "A1", "B2") */
using instanceID = std::string;

/** @brief Type alias for a colocation pattern (set of feature types) */
using Colocation = std::vector<FeatureType>;

/** @brief Type alias for a colocation instance (set of spatial instance pointers) */
using ColocationInstance = std::vector<const struct SpatialInstance*>;

/** @brief Type alias for colocation rules mapping */
using ColocationRule = std::unordered_map<Colocation, Colocation>;

// ============================================================================
// Data Structures
// ============================================================================

/**
 * @brief Structure representing a spatial data instance
 * 
 * Each spatial instance has a feature type, unique identifier, and 2D coordinates.
 */
struct SpatialInstance {
    FeatureType type;  ///< Feature type of this instance (e.g., "A", "B")
    instanceID id;     ///< Unique identifier (e.g., "A1", "B2")
    double x, y;       ///< 2D spatial coordinates
};

/**
 * @brief Structure representing a star neighborhood
 * 
 * A star neighborhood consists of a center instance and all its neighboring instances
 * within the distance threshold. This is a key concept in the joinless algorithm.
 */
struct StarNeighborhood {
    const SpatialInstance* center;                      ///< Center instance of the star
    std::vector<const SpatialInstance*> neighbors;      ///< All neighbors within distance threshold
};