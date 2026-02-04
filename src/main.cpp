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
#include <direct.h>
#include <iomanip>
#include <psapi.h>

int main(int argc, char* argv[]) {
    // Open file to write results
    std::ofstream out("mining_report.txt");
    std::streambuf* coutbuf = std::cout.rdbuf();
    std::cout.rdbuf(out.rdbuf());
    auto programStart = std::chrono::high_resolution_clock::now();

    // Print to terminal (not redirected)
    std::cerr << "=== MaxPR Colocation Mining Started ===\n";

    // ========================================================================
    // Step 1: Load Configuration
    // ========================================================================
    std::cerr << "[Step 1] Loading configuration...\n";
    std::cout << "[DEBUG] Step 1: Loading configuration...\n";
    std::cout << "[DEBUG] Current working directory: ";
    char cwd[1024];
    if (_getcwd(cwd, sizeof(cwd))) {
        std::cout << cwd << "\\n";
    }
    std::string config_path = (argc > 1) ? argv[1] : "config/config.txt";
    std::cout << "[DEBUG] Looking for config at: " << config_path << "\\n";
    AppConfig config = ConfigLoader::load(config_path);
    std::cout << "[DEBUG] Step 1: Configuration loaded successfully.\n";
    std::cout << "Running MaxPR with d=" << config.neighborDistance << "...\n";
    std::cerr << "[Step 1] Config loaded: d=" << config.neighborDistance << ", minMaxPR=" << config.minPrev << "\n";

    // ========================================================================
    // Step 2: Load Data
    // ========================================================================
    std::cerr << "[Step 2] Loading data...\n";
    std::cout << "[DEBUG] Step 2: Loading data from " << config.datasetPath << "...\n";
    auto t_load_start = std::chrono::high_resolution_clock::now();
    auto instances = DataLoader::load_csv(config.datasetPath);
    auto t_load_end = std::chrono::high_resolution_clock::now();
    std::cout << "[DEBUG] Step 2: Loaded " << instances.size() << " instances.\n";
    std::cerr << "[Step 2] Loaded " << instances.size() << " instances.\n";

    // ========================================================================
    // Step 3: Build Spatial Index
    // ========================================================================
    // Pass distance parameter d from config to spatial index
    std::cerr << "[Step 3] Building spatial index...\n";
    std::cout << "[DEBUG] Step 3: Building spatial index with d=" << config.neighborDistance << "...\n";
    auto t_idx_start = std::chrono::high_resolution_clock::now();
    SpatialIndex spatial_idx(config.neighborDistance);
    auto neighborPairs = spatial_idx.findNeighborPair(instances);
    auto t_idx_end = std::chrono::high_resolution_clock::now();
    std::cout << "[DEBUG] Step 3: Found " << neighborPairs.size() << " neighbor pairs.\n";
    std::cerr << "[Step 3] Found " << neighborPairs.size() << " neighbor pairs.\n";

    // ========================================================================
    // Step 4: Materialize Neighborhoods
    // ========================================================================
    std::cerr << "[Step 4] Building neighborhoods...\n";
    std::cout << "[DEBUG] Step 4: Materializing neighborhoods...\n";
    auto t_mat_start = std::chrono::high_resolution_clock::now();
    NeighborhoodMgr neighbor_mgr;
    neighbor_mgr.buildFromPairs(neighborPairs);
    auto t_mat_end = std::chrono::high_resolution_clock::now();
    std::cout << "[DEBUG] Step 4: Neighborhoods materialized.\n";
    std::cerr << "[Step 4] Neighborhoods built.\n";

    // ========================================================================
    // Step 5: Mine Colocation Patterns
    // ========================================================================
    std::cerr << "[Step 5] Mining colocations (this may take a while)...\n";
    std::cout << "[DEBUG] Step 5: Mining colocation patterns with minMaxPR=" << config.minPrev << "...\n";
    MaxPRMiner miner;
    
    // Define progress callback lambda function
    // This callback reports mining progress to the TERMINAL (stderr)
    auto progressCallback = [](int currentStep, int totalSteps, const std::string& message, double percentage) {
        std::cerr << "\r[MINING] " << std::fixed << std::setprecision(1) << percentage 
                  << "% (k=" << currentStep << ") - " << message << "          ";
        std::cerr.flush();
        
        // Print newline when completed
        if (percentage >= 100.0) {
            std::cerr << std::endl;
        }
    };
    
    auto colocations = miner.mineColocations(config.minPrev, &neighbor_mgr, instances, progressCallback);
    std::cout << "[DEBUG] Step 5: Mining completed.\n";
    std::cerr << "[Step 5] Mining completed! Found " << colocations.size() << " patterns.\n";
    
    // ========================================================================
    // Print Results
    // ========================================================================
    auto programEnd = std::chrono::high_resolution_clock::now();
    double totalTimeSec = std::chrono::duration<double>(programEnd - programStart).count();
    double maxMemory = getMemoryUsageMB();

    std::cout << "\n==========================================================\n";
    std::cout << "                   FINAL SUMMARY                          \n";
    std::cout << "==========================================================\n";
    std::cout << "Total MaxPR Colocations Found: " << colocations.size() << "\n";
    std::cout << "Total Execution Time: " << std::fixed << std::setprecision(4) << totalTimeSec << " seconds\n";
    std::cout << "Peak Memory Usage (Approx): " << std::fixed << std::setprecision(2) << maxMemory << " MB\n";

    std::cout << "\n[DETAILED RESULTS]\n";
    for (const auto& colPair : colocations) {
        const auto& col = colPair.first;
        double maxPR = colPair.second;
        for (size_t i = 0; i < col.size(); ++i) {
            if (i > 0) std::cout << " - ";
            std::cout << col[i];
        }
        std::cout << " (maxPR: " << std::fixed << std::setprecision(4) << maxPR << ")\n";
    }

    std::cout.rdbuf(coutbuf);
    std::cout << "Mining completed. Report saved to 'mining_report.txt'.\n";

    return 0;
}