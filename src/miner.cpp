/**
 * @file miner.cpp
 * @brief Implementation of the maxPR colocation pattern mining algorithm
 * 
 * This file implements the maxPR algorithm that discovers spatial colocation
 * patterns with rare events using maximal participation ratio measure.
 * Based on: "Mining Co-Location Patterns with Rare Events" (Huang et al., 2006)
 */

#include "miner.h"
#include "utils.h"
#include "neighborhood_mgr.h"
#include "types.h"
#include <algorithm>
#include <set>
#include <string>
#include <iostream>
#include <omp.h> 
#include <iomanip>
#include <chrono>

std::vector<std::pair<Colocation, double>> MaxPRMiner::mineColocations(
    double minMaxPR, 
    NeighborhoodMgr* neighborhoodMgr, 
    const std::vector<SpatialInstance>& instances,
    ProgressCallback progressCb
) {
	// Start timer
    auto minerStart = std::chrono::high_resolution_clock::now();

    // Assign parameters to member variables for use in other methods
    this->progressCallback = progressCb;
    
    // Initialize mining variables
    int k = 2;  // Start with size-2 patterns
    std::vector<FeatureType> types = getAllObjectTypes(instances);
    std::map<FeatureType, int> featureCount = countInstancesByFeature(instances);
    std::vector<Colocation> prevColocations;
    std::vector<ColocationInstance> cliqueInstances;
    std::vector<ColocationInstance> prevCliqueInstances;
    std::vector<std::pair<Colocation, double>> allMaxPRColocations;  // Store pattern with maxPR value
    std::map<Colocation, double> patternMaxPR;  // Track maxPR values for candidate generation

    // Estimate total iterations (max pattern size is number of types)
    int maxK = static_cast<int>(types.size());
    int currentIteration = 0;
    int totalIterations = 0; // Will be updated as we go

    if (progressCallback) {
        progressCallback(0, maxK, "Initializing maxPR mining process...", 0.0);
    }

    // Initialize with size-1 patterns (individual feature types)
    for (auto t : types) {
        prevColocations.push_back({t});
        patternMaxPR[{t}] = 1.0;  // Size-1 patterns have maxPR = 1
    }

    while (!prevColocations.empty()) {
        currentIteration++;
        totalIterations = currentIteration;
        
        // Calculate progress: use a conservative estimate that we might go up to maxK
        // But cap at 95% until we're actually done
        double progressPercent = std::min(95.0, (static_cast<double>(currentIteration) / maxK) * 95.0);
        
        if (progressCallback) {
            progressCallback(currentIteration, maxK, 
                "Processing k=" + std::to_string(k) + " patterns...", 
                progressPercent);
        }
        std::vector<ColocationInstance> starInstances;

		// 1. Generate candidate patterns of size k
        auto t1_start = std::chrono::high_resolution_clock::now();
        std::vector<Colocation> candidates;
        
        if (k == 2) {
            // Step 1 of maxPrune: generate C2 by geometric methods (all pairs)
            for (size_t i = 0; i < types.size(); i++) {
                for (size_t j = i + 1; j < types.size(); j++) {
                    candidates.push_back({types[i], types[j]});
                }
            }
        } else {
            // Step 3 of maxPrune: generate Ck+1 using weak monotonicity
            candidates = generateCandidatesMaxPR(prevColocations, patternMaxPR, minMaxPR);
        }
        auto t1_end = std::chrono::high_resolution_clock::now();
        printDuration("generateCandidates (k=" + std::to_string(k) + ")", t1_start, t1_end);

        if (candidates.empty()) {
            if (progressCallback) {
                progressCallback(currentIteration, maxK, 
                    "No more candidates found. Mining completed.", 
                    100.0);
            }
            break;
        }
        
        if (progressCallback) {
            double progressPercent = std::min(95.0, (static_cast<double>(currentIteration) / maxK) * 95.0);
            progressCallback(currentIteration, maxK, 
                "Filtering star instances for " + std::to_string(candidates.size()) + " candidates...", 
                progressPercent);
        }
        
		// 2. Filter star instances for each candidate
        auto t2_start = std::chrono::high_resolution_clock::now();
        for (auto t : types) {
            for (const auto& starNeigh : neighborhoodMgr->getAllStarNeighborhoods()) {
                if (starNeigh.first == t) {
                    std::vector<ColocationInstance> found = filterStarInstances(candidates, starNeigh);
                    starInstances.insert(starInstances.end(), found.begin(), found.end());
                }
            }
        }
        auto t2_end = std::chrono::high_resolution_clock::now();
        printDuration("filterStarInstances (Total) (k=" + std::to_string(k) + ")", t2_start, t2_end);

        if (k==2){
            cliqueInstances = starInstances;
            if (progressCallback) {
                double progressPercent = std::min(95.0, (static_cast<double>(currentIteration) / maxK) * 95.0);
                progressCallback(currentIteration, maxK, 
                    "Found " + std::to_string(starInstances.size()) + " star instances (k=2)...", 
                    progressPercent);
            }
        }else{
            if (progressCallback) {
                double progressPercent = std::min(95.0, (static_cast<double>(currentIteration) / maxK) * 95.0);
                progressCallback(currentIteration, maxK, 
                    "Selecting maxPR colocations (coarse filter)...", 
                    progressPercent);
            }

			// 3. Select colocations using maxPR coarse filter
            auto t3_start = std::chrono::high_resolution_clock::now();
            std::map<Colocation, double> tempMaxPR;
            candidates = selectMaxPRColocations(candidates, starInstances, minMaxPR, featureCount, tempMaxPR);
            auto t3_end = std::chrono::high_resolution_clock::now();
            printDuration("selectMaxPRColocations (Coarse) (k=" + std::to_string(k) + ")", t3_start, t3_end);
            
            if (progressCallback) {
                double progressPercent = std::min(95.0, (static_cast<double>(currentIteration) / maxK) * 95.0);
                progressCallback(currentIteration, maxK, 
                    "Filtering clique instances...", 
                    progressPercent);
            }

			// 4. Filter clique instances for candidates
            auto t4_start = std::chrono::high_resolution_clock::now();
            cliqueInstances = filterCliqueInstances(
                candidates,
                starInstances,
                prevCliqueInstances
            );
            auto t4_end = std::chrono::high_resolution_clock::now();
            printDuration("filterCliqueInstances (k=" + std::to_string(k) + ")", t4_start, t4_end);
        }
        
        if (progressCallback) {
            double progressPercent = std::min(95.0, (static_cast<double>(currentIteration) / maxK) * 95.0);
            progressCallback(currentIteration, maxK, 
                "Selecting final maxPR colocations...", 
                progressPercent);
        }

        auto t5_start = std::chrono::high_resolution_clock::now();
        std::map<Colocation, double> currentMaxPR;
        prevColocations = selectMaxPRColocations(
            candidates,
            cliqueInstances,
            minMaxPR,
            featureCount,
            currentMaxPR
        );
        // Update global maxPR map
        for (const auto& item : currentMaxPR) {
            patternMaxPR[item.first] = item.second;
        }
        auto t5_end = std::chrono::high_resolution_clock::now();
        printDuration("selectMaxPRColocations (Final) (k=" + std::to_string(k) + ")", t5_start, t5_end);
        

        if (!prevColocations.empty()) {
             // Add patterns with their maxPR values
             for (const auto& col : prevColocations) {
                 double maxPRVal = currentMaxPR.count(col) ? currentMaxPR[col] : 0.0;
                 allMaxPRColocations.push_back({col, maxPRVal});
             }

             std::cout << "[DEBUG] Found " << prevColocations.size() << " maxPR patterns. Details:\n";
             
             if (progressCallback) {
                 double progressPercent = std::min(95.0, (static_cast<double>(currentIteration) / maxK) * 95.0);
                 progressCallback(currentIteration, maxK, 
                     "Found " + std::to_string(prevColocations.size()) + " maxPR k=" + std::to_string(k) + " colocations", 
                     progressPercent);
             }
        } else {
            if (progressCallback) {
                double progressPercent = std::min(95.0, (static_cast<double>(currentIteration) / maxK) * 95.0);
                progressCallback(currentIteration, maxK, 
                    "No maxPR k=" + std::to_string(k) + " colocations found", 
                    progressPercent);
            }
        }
        
        prevCliqueInstances = std::move(cliqueInstances);

        k++;
    }
    
    if (progressCallback) {
        progressCallback(maxK, maxK, 
            "Mining completed! Total maxPR colocations: " + std::to_string(allMaxPRColocations.size()), 
            100.0);
    }

    auto minerEnd = std::chrono::high_resolution_clock::now();
    printDuration("TOTAL MINING TIME", minerStart, minerEnd);
    
    return allMaxPRColocations;
}


std::vector<Colocation> MaxPRMiner::generateCandidatesMaxPR(
    const std::vector<Colocation>& prevMaxPR,
    const std::map<Colocation, double>& patternMaxPR,
    double minMaxPR) 
{
    std::vector<Colocation> candidates;
    
    if (prevMaxPR.empty()) {
        return candidates;
    }
    
    size_t patternSize = prevMaxPR[0].size();
    
    // Create set for quick lookup of patterns above threshold
    std::set<Colocation> prevSet(prevMaxPR.begin(), prevMaxPR.end());
    
    // ========================================================================
    // WEAK MONOTONICITY CANDIDATE GENERATION (Algorithm maxPrune - Step 3)
    // Based on Example 8 and Lemma 3 from paper:
    // 
    // Lemma 3: For a k-pattern P above threshold (k >= 3), there exist two
    // (k-1) patterns P1, P2 such that:
    // (1) P1 ⊂ P, P2 ⊂ P
    // (2) P1 and P2 share their first k-2 features  
    // (3) P1 and P2 share either the kth or (k-1)th or (k-2)th feature in P
    //     but not any two of them
    // (4) Both P1 and P2 are above threshold
    //
    // Therefore: Join P1 and P2 if they differ by one feature in last two positions
    // ========================================================================
    
    for (size_t i = 0; i < prevMaxPR.size(); i++) {
        for (size_t j = i + 1; j < prevMaxPR.size(); j++) {
            const auto& p1 = prevMaxPR[i];
            const auto& p2 = prevMaxPR[j];
            
            // Check if they share first (k-2) features
            bool sharePrefix = true;
            for (size_t idx = 0; idx + 2 < patternSize; idx++) {
                if (p1[idx] != p2[idx]) {
                    sharePrefix = false;
                    break;
                }
            }
            
            if (!sharePrefix) continue;
            
            // Generate candidate by union of p1 and p2
            std::set<FeatureType> candidateSet(p1.begin(), p1.end());
            candidateSet.insert(p2.begin(), p2.end());
            
            // Must produce exactly (k+1) features
            if (candidateSet.size() != patternSize + 1) {
                continue;
            }
            
            Colocation candidate(candidateSet.begin(), candidateSet.end());
            
            // WEAK MONOTONICITY PRUNING (Lemma 2):
            // At most ONE (k-1)-subpattern can be below threshold
            int belowThresholdCount = 0;
            
            for (size_t idx = 0; idx < candidate.size(); idx++) {
                Colocation subset = candidate; 
                subset.erase(subset.begin() + idx);
                
                if (prevSet.find(subset) == prevSet.end()) {
                    belowThresholdCount++;
                    if (belowThresholdCount > 1) break;  // Early exit
                }
            }
            
            if (belowThresholdCount <= 1) {
                candidates.push_back(candidate);
            }
        }
    }
    
    // Remove duplicates
    std::sort(candidates.begin(), candidates.end());
    candidates.erase(std::unique(candidates.begin(), candidates.end()), 
                     candidates.end());

    return candidates;
}


std::vector<ColocationInstance> MaxPRMiner::filterStarInstances(
    const std::vector<Colocation>& candidates, 
    const std::pair<FeatureType, std::vector<StarNeighborhood>>& starNeigh) 
{
    std::vector<ColocationInstance> filteredInstances;
    FeatureType centerType = starNeigh.first;
    
    // Filter candidates to only those with this center type as first element
    std::vector<const Colocation*> relevantCandidates;
    for (const auto& cand : candidates) {
        if (!cand.empty() && cand[0] == centerType) {
            relevantCandidates.push_back(&cand);
        }
    }

    if (relevantCandidates.empty()) return filteredInstances;

    // Iterate through each star neighborhood
    for (const auto& star : starNeigh.second) {
        
        // Build a map of neighbors by feature type for fast lookup
        // Using const pointer since star is const
        std::unordered_map<FeatureType, std::vector<const SpatialInstance*>> neighborMap;
        
        for (auto neighbor : star.neighbors) {
            neighborMap[neighbor->type].push_back(neighbor);
        }

        // Check each relevant candidate pattern
        for (const auto* candPtr : relevantCandidates) {
            const auto& candidate = *candPtr;
            
            std::vector<const SpatialInstance*> currentInstance;
            currentInstance.reserve(candidate.size());
            
            // Add center instance as first element
            currentInstance.push_back(star.center);

			// Recursive function to find combinations
            findCombinations(candidate, 1, currentInstance, neighborMap, filteredInstances);
        }
    }

    return filteredInstances;
}


std::vector<ColocationInstance> MaxPRMiner::filterCliqueInstances(
    const std::vector<Colocation>& candidates,
    const std::vector<ColocationInstance>& instances,
    const std::vector<ColocationInstance>& prevInstances
) {
    // ========================================================================
	// STEP 1: PREPARE LOOKUP STRUCTURES
    // ========================================================================

	// 1.1. Create a set of valid candidate patterns for quick lookup
    std::set<Colocation> validCandidatePatterns(candidates.begin(), candidates.end());

	// 1.2. Create a set of previous instances for quick lookup
    std::set<std::vector<std::string>> validPrevIds;
    for (const auto& prevInst : prevInstances) {
        std::vector<std::string> ids;
        ids.reserve(prevInst.size());
        for (const auto* ptr : prevInst) {
            ids.push_back(ptr->id);
        }
        validPrevIds.insert(ids);
    }

    // ========================================================================
	// STEP 2: PREPARE THREAD BUFFERS
    // ========================================================================
    int num_threads = omp_get_max_threads();
    std::vector<std::vector<ColocationInstance>> thread_buffers(num_threads);

    // ========================================================================
	// STEP 3: PARALLEL FILTERING
    // ========================================================================
    #pragma omp parallel for
    for (size_t i = 0; i < instances.size(); ++i) {
        const auto& instance = instances[i];
        int thread_id = omp_get_thread_num();

        // Safety check
        if (instance.size() < 2) continue;

        Colocation currentPattern;
        currentPattern.reserve(instance.size());
        for (const auto* ptr : instance) {
            currentPattern.push_back(ptr->type);
        }

		// If current pattern is not a valid candidate, skip
        if (validCandidatePatterns.find(currentPattern) == validCandidatePatterns.end()) {
            continue;
        }

        std::vector<std::string> subInstanceIds;
        subInstanceIds.reserve(instance.size() - 1);

		// Generate (k-1)-subset by removing the first instance
        for (size_t j = 1; j < instance.size(); ++j) {
            subInstanceIds.push_back(instance[j]->id);
        }

		// Check if the (k-1)-subset exists in previous instances
        if (validPrevIds.find(subInstanceIds) != validPrevIds.end()) {
            thread_buffers[thread_id].push_back(instance);
        }
    }

    // ========================================================================
	// STEP 4: COMBINE RESULTS
    // ========================================================================
    std::vector<ColocationInstance> filteredInstances;
    for (const auto& buffer : thread_buffers) {
        filteredInstances.insert(filteredInstances.end(), buffer.begin(), buffer.end());
    }

    return filteredInstances;
}


std::vector<Colocation> MaxPRMiner::selectMaxPRColocations(
    const std::vector<Colocation>& candidates, 
    const std::vector<ColocationInstance>& instances, 
    double minMaxPR, 
    const std::map<FeatureType, int>& featureCount,
    std::map<Colocation, double>& patternMaxPR) 
{
    std::vector<Colocation> maxPRColocations;

    // ========================================================================
    // STEP 1: Data structure for aggregation
    // ========================================================================
    // Key: Candidate (Pattern)
    // Value: Map<FeatureType, Set<InstanceID>> - count unique instances for each feature
    std::map<Colocation, std::map<FeatureType, std::set<std::string>>> candidateStats;

    // Initialize stats map for all candidates (ensures every candidate has an entry even without instances)
    // This step costs O(C), very fast compared to O(C*I)
    for (const auto& cand : candidates) {
        candidateStats[cand]; // Create empty entry
    }

    // ========================================================================
    // STEP 2: Single pass through instances
    // ========================================================================
    // Time complexity: O(I * K * log(K)) where K is pattern length (typically very small)
    for (const ColocationInstance& instance : instances) {
        // 2a. Extract pattern from instance
        // Example: Instance has A1, B1 -> Pattern is {A, B}
        Colocation patternKey;
        for (const auto& instPtr : instance) {
            patternKey.push_back(instPtr->type);
        }
        // Ensure patternKey is sorted to match key in map (if candidate is already sorted)
        // std::sort(patternKey.begin(), patternKey.end()); 

        // 2b. Check if this pattern is in the candidates of interest
        auto it = candidateStats.find(patternKey);
        if (it != candidateStats.end()) {
            // 2c. Update participating instances statistics
            // it->second is map<Feature, Set<ID>>
            for (const auto& instPtr : instance) {
                it->second[instPtr->type].insert(instPtr->id);
            }
        }
    }

    // ========================================================================
    // STEP 3: Calculate maxPR and filter (iterate through candidates)
    // ========================================================================
    // Time complexity: O(C * K)
    for (const auto& item : candidateStats) {
        const Colocation& candidate = item.first;
        const auto& participatingMap = item.second; // Map<Feature, Set<ID>>

        double max_participation_ratio = 0.0;  // Changed from min to max
        bool possible = true;
        
        // If no instances participate -> participatingMap may be missing features or have empty sets
        
        for (const auto& feature : candidate) {
            // Get total global instance count for this feature
            auto totalIt = featureCount.find(feature);
            if (totalIt == featureCount.end() || totalIt->second == 0) {
                possible = false;
                break;
            }
            double totalFeatureCount = (double)(totalIt->second);

            // Get count of instances participating in pattern
            int participatedCount = 0;
            auto partIt = participatingMap.find(feature);
            if (partIt != participatingMap.end()) {
                participatedCount = static_cast<int>(partIt->second.size());
            }

            double ratio = (double)participatedCount / totalFeatureCount;
            // maxPR: take maximum instead of minimum
            if (ratio > max_participation_ratio) {
                max_participation_ratio = ratio;
            }
        }

        // Store maxPR value for this pattern
        patternMaxPR[candidate] = max_participation_ratio;

        if (possible && max_participation_ratio >= minMaxPR) {
            maxPRColocations.push_back(candidate);
        }
    }

    return maxPRColocations;
}