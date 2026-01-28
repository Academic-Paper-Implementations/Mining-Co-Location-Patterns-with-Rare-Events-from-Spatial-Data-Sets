/**
 * @file miner.cpp
 * @brief Implementation of the joinless colocation pattern mining algorithm
 * 
 * This file implements the core mining algorithm that discovers prevalent spatial
 * colocation patterns using star neighborhoods to avoid expensive join operations.
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

std::vector<Colocation> JoinlessMiner::mineColocations(
    double minPrev, 
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
    std::vector<Colocation> allPrevalentColocations;

    // Estimate total iterations (max pattern size is number of types)
    int maxK = types.size();
    int currentIteration = 0;
    int totalIterations = 0; // Will be updated as we go

    if (progressCallback) {
        progressCallback(0, maxK, "Initializing mining process...", 0.0);
    }

    // Initialize with size-1 patterns (individual feature types)
    for (auto t : types) prevColocations.push_back({t});

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
        std::vector<Colocation> candidates = generateCandidates(prevColocations);
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
                    "Selecting prevalent colocations (coarse filter)...", 
                    progressPercent);
            }

			// 3. Select prevalent colocations using coarse filter
            auto t3_start = std::chrono::high_resolution_clock::now();
            candidates = selectPrevColocations(candidates, starInstances, minPrev, featureCount);
            auto t3_end = std::chrono::high_resolution_clock::now();
            printDuration("selectPrevColocations (Coarse) (k=" + std::to_string(k) + ")", t3_start, t3_end);
            
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
                "Selecting final prevalent colocations...", 
                progressPercent);
        }

        auto t5_start = std::chrono::high_resolution_clock::now();
        prevColocations = selectPrevColocations(
            candidates,
            cliqueInstances,
            minPrev,
            featureCount
        );
        auto t5_end = std::chrono::high_resolution_clock::now();
        printDuration("selectPrevColocations (Final) (k=" + std::to_string(k) + ")", t5_start, t5_end);
        

        if (!prevColocations.empty()) {
             allPrevalentColocations.insert(
                 allPrevalentColocations.end(), 
                 prevColocations.begin(), 
                 prevColocations.end()
             );

             std::cout << "[DEBUG] Found " << prevColocations.size() << " prevalent patterns. Details:\n";
             
             if (progressCallback) {
                 double progressPercent = std::min(95.0, (static_cast<double>(currentIteration) / maxK) * 95.0);
                 progressCallback(currentIteration, maxK, 
                     "Found " + std::to_string(prevColocations.size()) + " prevalent k=" + std::to_string(k) + " colocations", 
                     progressPercent);
             }
        } else {
            if (progressCallback) {
                double progressPercent = std::min(95.0, (static_cast<double>(currentIteration) / maxK) * 95.0);
                progressCallback(currentIteration, maxK, 
                    "No prevalent k=" + std::to_string(k) + " colocations found", 
                    progressPercent);
            }
        }
        
        prevCliqueInstances = std::move(cliqueInstances);

        k++;
    }
    
    if (progressCallback) {
        progressCallback(maxK, maxK, 
            "Mining completed! Total prevalent colocations: " + std::to_string(allPrevalentColocations.size()), 
            100.0);
    }

    auto minerEnd = std::chrono::high_resolution_clock::now();
    printDuration("TOTAL MINING TIME", minerStart, minerEnd);
    
    return allPrevalentColocations;
}


std::vector<Colocation> JoinlessMiner::generateCandidates(
    const std::vector<Colocation>& prevPrevalent) 
{
    std::vector<Colocation> candidates;
    
    if (prevPrevalent.empty()) {
        return candidates;
    }
    
    size_t patternSize = prevPrevalent[0].size();
    
    std::set<Colocation> prevSet(prevPrevalent.begin(), prevPrevalent.end());
    
    // Generate candidate
    for (size_t i = 0; i < prevPrevalent.size(); i++) {
        for (size_t j = i + 1; j < prevPrevalent.size(); j++) {
            // Take prefix of k-1 first element
            Colocation prefix1(prevPrevalent[i].begin(), 
                             prevPrevalent[i].end() - 1);
            Colocation prefix2(prevPrevalent[j].begin(), 
                             prevPrevalent[j].end() - 1);
            
            // Just join when the prefix is equal
            if (prefix1 != prefix2) {
                continue;
            }
            
            // Generate new candidate
            std::set<FeatureType> candidateSet(prevPrevalent[i].begin(), 
                                              prevPrevalent[i].end());
            candidateSet.insert(prevPrevalent[j].back());
            
            if (candidateSet.size() != patternSize + 1) {
                continue;
            }
            
            Colocation candidate(candidateSet.begin(), candidateSet.end());
            
            // APRIORI PRUNING
            bool allSubsetsValid = true;
            std::vector<FeatureType> candFeatures = candidate;
            
            for (size_t idx = 0; idx < candFeatures.size(); idx++) {
                Colocation subset = candFeatures; 
                subset.erase(subset.begin() + idx);
                
                if (prevSet.find(subset) == prevSet.end()) {
                    allSubsetsValid = false;
                    break;
                }
            }
            
            if (allSubsetsValid) {
                candidates.push_back(candidate);
            }
        }
    }
    
    // Remove duplicate
    std::sort(candidates.begin(), candidates.end());
    candidates.erase(std::unique(candidates.begin(), candidates.end()), 
                     candidates.end());

    return candidates;
}


std::vector<ColocationInstance> JoinlessMiner::filterStarInstances(
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


std::vector<ColocationInstance> JoinlessMiner::filterCliqueInstances(
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


std::vector<Colocation> JoinlessMiner::selectPrevColocations(
    const std::vector<Colocation>& candidates, 
    const std::vector<ColocationInstance>& instances, 
    double minPrev, 
    const std::map<FeatureType, int>& featureCount) 
{
    std::vector<Colocation> coarsePrevalent;

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
    // STEP 3: Calculate ratios and filter (iterate through candidates)
    // ========================================================================
    // Time complexity: O(C * K)
    for (const auto& item : candidateStats) {
        const Colocation& candidate = item.first;
        const auto& participatingMap = item.second; // Map<Feature, Set<ID>>

        double min_participation_ratio = 1.0;
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
                participatedCount = partIt->second.size();
            }

            double ratio = (double)participatedCount / totalFeatureCount;
            if (ratio < min_participation_ratio) {
                min_participation_ratio = ratio;
            }
        }

        if (possible && min_participation_ratio >= minPrev) {
            coarsePrevalent.push_back(candidate);
        }
    }

    return coarsePrevalent;
}