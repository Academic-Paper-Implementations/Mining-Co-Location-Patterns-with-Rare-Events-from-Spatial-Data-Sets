/**
 * @file main.cpp
 * @brief Entry point for the joinless colocation pattern mining application
 */

#include "config.h"
#include "data_loader.h"
#include "spatial_index.h"
#include "neighborhood_mgr.h"
#include "miner.h"
#include "utils.h"
#include <iostream>
#include <fstream>
#include <chrono>
#include <windows.h>
#include <iomanip>
#include <psapi.h>

int main(int argc, char* argv[]) {
    // Open file to write results
    std::ofstream out("D:/tai_lieu_hoc_AI/spatial_data_mining/A-Joinless-Approach-for-Mining-Spatial-Colocation-Patterns/mining_report.txt");
    std::streambuf* coutbuf = std::cout.rdbuf();
    std::cout.rdbuf(out.rdbuf());
    auto programStart = std::chrono::high_resolution_clock::now();

    // ========================================================================
    // Step 1: Load Configuration
    // ========================================================================
    std::cout << "[DEBUG] Step 1: Loading configuration...\n";
    std::string config_path = (argc > 1) ? argv[1] : "src/c++/config.txt";
    AppConfig config = ConfigLoader::load(config_path);
    std::cout << "[DEBUG] Step 1: Configuration loaded successfully.\n";
    std::cout << "Running Joinless with d=" << config.neighborDistance << "...\n";

    // ========================================================================
    // Step 2: Load Data
    // ========================================================================
    std::cout << "[DEBUG] Step 2: Loading data from " << config.datasetPath << "...\n";
    auto t_load_start = std::chrono::high_resolution_clock::now();
    auto instances = DataLoader::load_csv(config.datasetPath);
    auto t_load_end = std::chrono::high_resolution_clock::now();
    std::cout << "[DEBUG] Step 2: Loaded " << instances.size() << " instances.\n";

    // ========================================================================
    // Step 3: Build Spatial Index
    // ========================================================================
    // Pass distance parameter d from config to spatial index
    std::cout << "[DEBUG] Step 3: Building spatial index with d=" << config.neighborDistance << "...\n";
    auto t_idx_start = std::chrono::high_resolution_clock::now();
    SpatialIndex spatial_idx(config.neighborDistance);
    auto neighborPairs = spatial_idx.findNeighborPair(instances);
    auto t_idx_end = std::chrono::high_resolution_clock::now();
    std::cout << "[DEBUG] Step 3: Found " << neighborPairs.size() << " neighbor pairs.\n";

    // ========================================================================
    // Step 4: Materialize Neighborhoods
    // ========================================================================
    std::cout << "[DEBUG] Step 4: Materializing neighborhoods...\n";
    auto t_mat_start = std::chrono::high_resolution_clock::now();
    NeighborhoodMgr neighbor_mgr;
    neighbor_mgr.buildFromPairs(neighborPairs);
    auto t_mat_end = std::chrono::high_resolution_clock::now();
    std::cout << "[DEBUG] Step 4: Neighborhoods materialized.\n";

    // ========================================================================
    // Step 5: Mine Colocation Patterns
    // ========================================================================
    std::cout << "[DEBUG] Step 5: Mining colocation patterns with minPrev=" << config.minPrev << "...\n";
    JoinlessMiner miner;
    
    // Define progress callback lambda function
    // This callback reports mining progress to the console
    auto progressCallback = [](int currentStep, int totalSteps, const std::string& message, double percentage) {
        std::cout << "\r[PROGRESS] " << std::fixed << std::setprecision(1) << percentage 
                  << "% (" << currentStep << "/" << totalSteps << ") - " << message;
        std::cout.flush();
        
        // Print newline when completed
        if (percentage >= 100.0) {
            std::cout << std::endl;
        }
    };
    
    auto colocations = miner.mineColocations(config.minPrev, &neighbor_mgr, instances, progressCallback);
    std::cout << "[DEBUG] Step 5: Mining completed.\n";
    
    // ========================================================================
    // Print Results
    // ========================================================================
    auto programEnd = std::chrono::high_resolution_clock::now();
    double totalTimeSec = std::chrono::duration<double>(programEnd - programStart).count();
    double maxMemory = getMemoryUsageMB();

    std::cout << "\n==========================================================\n";
    std::cout << "                   FINAL SUMMARY                          \n";
    std::cout << "==========================================================\n";
    std::cout << "Total Prevalent Colocations Found: " << colocations.size() << "\n";
    std::cout << "Total Execution Time: " << std::fixed << std::setprecision(4) << totalTimeSec << " seconds\n";
    std::cout << "Peak Memory Usage (Approx): " << std::fixed << std::setprecision(2) << maxMemory << " MB\n";

    std::cout << "\n[DETAILED RESULTS]\n";
    for (const auto& col : colocations) {
        for (size_t i = 0; i < col.size(); ++i) {
            if (i > 0) std::cout << " - ";
            std::cout << col[i];
        }
        std::cout << "\n";
    }

    std::cout.rdbuf(coutbuf);
    std::cout << "Mining completed. Report saved to 'mining_report.txt'.\n";

    return 0;
}