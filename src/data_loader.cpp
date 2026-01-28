/**
 * @file data_loader.cpp
 * @brief Implementation of CSV data loading for spatial instances
 */

#include "data_loader.h"
#include <iostream>

using namespace csv;

std::vector<SpatialInstance> DataLoader::load_csv(const std::string& filepath) {
    CSVReader reader(filepath);
    std::vector<SpatialInstance> instances;

    // Parse each row from the CSV file
    for (auto& row : reader) {
        SpatialInstance instance;
        
        // Map CSV columns to SpatialInstance fields
        instance.type = row["Feature"].get<FeatureType>();
        
        // Generate instance ID by concatenating feature type and instance number
        // Example: Feature "A" with Instance 1 becomes "A1"
        instance.id = instanceID(instance.type + std::to_string(row["Instance"].get<int>()));
        
        instance.x = row["LocX"].get<double>();
        instance.y = row["LocY"].get<double>();
        
        instances.push_back(instance);
    }
    
    return instances;
}