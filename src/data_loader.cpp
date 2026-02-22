/**
 * @file data_loader.cpp
 * @brief Implementation of CSV data loading for spatial instances
 */

#include "data_loader.h"
#include <iostream>
#include <algorithm>
#include <random>
#include <map>

using namespace csv;


/**
 * @brief Load spatial instances from a CSV file
 * @param filepath Path to the CSV file
 * @return std::vector<SpatialInstance> Vector of loaded spatial instances
 * 
 * Expects CSV with columns: Feature, Instance, LocX, LocY.
 * Instance IDs are generated as: FeatureType + InstanceNumber (e.g., "A1", "B2").
 */
std::vector<SpatialInstance> DataLoader::load_csv(const std::string& filepath, double percentage) {
    CSVReader reader(filepath);
    auto colNames = reader.get_col_names();
    std::string xCol = "LocX";
    std::string yCol = "LocY";
    auto hasColumn = [&](const std::string& name) {
        return std::find(colNames.begin(), colNames.end(), name) != colNames.end();
        };

    if (hasColumn("X")) xCol = "X";
    if (hasColumn("Y")) yCol = "Y";

    std::vector<SpatialInstance> allInstances;

    for (auto& row : reader) {
        SpatialInstance instance;

        instance.type = row["Feature"].get<FeatureType>();
        instance.id = instance.type + std::to_string(row["Instance"].get<int>());
        instance.x = row[xCol].get<double>();
        instance.y = row[yCol].get<double>();

        allInstances.push_back(instance);
    }

    if (percentage >= 1.0 || percentage <= 0.0) {
        return allInstances;
    }

    std::map<std::string, std::vector<SpatialInstance>> featureGroups;
    for (const auto& inst : allInstances) {
        featureGroups[inst.type].push_back(inst);
    }

    std::vector<SpatialInstance> sampledInstances;

    std::random_device rd;
    std::mt19937 g(rd());

    std::cout << "Sampling " << (percentage * 100) << "% data per feature...\n";

    for (auto& pair : featureGroups) {
        auto& group = pair.second;

        std::shuffle(group.begin(), group.end(), g);

        size_t keepCount = static_cast<size_t>(group.size() * percentage);

        if (keepCount == 0 && !group.empty() && percentage > 0) {
            keepCount = 1;
        }

        for (size_t i = 0; i < keepCount; ++i) {
            sampledInstances.push_back(group[i]);
        }
    }

    std::cout << "Reduced dataset from " << allInstances.size()
        << " to " << sampledInstances.size() << " instances.\n";

    return sampledInstances;
}