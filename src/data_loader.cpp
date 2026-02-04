/**
 * @file data_loader.cpp
 * @brief Implementation of CSV data loading for spatial instances
 */

#include "data_loader.h"
#include <iostream>

using namespace csv;


/**
 * @brief Load spatial instances from a CSV file
 * @param filepath Path to the CSV file
 * @return std::vector<SpatialInstance> Vector of loaded spatial instances
 * 
 * Expects CSV with columns: Feature, Instance, LocX, LocY.
 * Instance IDs are generated as: FeatureType + InstanceNumber (e.g., "A1", "B2").
 */
std::vector<SpatialInstance> DataLoader::load_csv(const std::string& filepath) {
    CSVReader reader(filepath);
    std::vector<SpatialInstance> instances;

    for (auto& row : reader) {
        SpatialInstance instance;
        
        instance.type = row["Feature"].get<FeatureType>();
        instance.id = instanceID(instance.type + std::to_string(row["Instance"].get<int>()));
        instance.x = row["LocX"].get<double>();
        instance.y = row["LocY"].get<double>();
        
        instances.push_back(instance);
    }
    
    return instances;
}