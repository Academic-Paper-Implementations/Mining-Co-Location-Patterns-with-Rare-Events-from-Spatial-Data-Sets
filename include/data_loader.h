/**
 * @file data_loader.h
 * @brief CSV data loading functionality for spatial instances
 */

#pragma once
#include "types.h"
#include "csv.hpp"
#include <string>
#include <vector>

/**
 * @brief DataLoader class for loading spatial instances from CSV files
 * 
 * Provides static methods to parse CSV datasets containing spatial feature instances.
 */
class DataLoader {
public:
    /**
     * @brief Load spatial instances from a CSV file
     * 
     * Expects CSV with columns: Feature, Instance, LocX, LocY
     * - Feature: Feature type (e.g., "A", "B", "Restaurant")
     * - Instance: Instance number (integer)
     * - LocX: X coordinate (double)
     * - LocY: Y coordinate (double)
     * 
     * @param filepath Path to the CSV file
     * @return std::vector<SpatialInstance> Vector of loaded spatial instances
     * @note Instance IDs are generated as: FeatureType + InstanceNumber (e.g., "A1", "B2")
     */
    static std::vector<SpatialInstance> load_csv(const std::string& filepath);
};