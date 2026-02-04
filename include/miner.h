/**
 * @file miner.h
 * @brief MaxPR colocation pattern mining algorithm implementation
 * 
 * This file contains the core mining algorithm that discovers spatial
 * colocation patterns with rare events using maximal participation ratio.
 * Based on: "Mining Co-Location Patterns with Rare Events" (Huang et al., 2006)
 */

#pragma once
#include "types.h"
#include "neighborhood_mgr.h"
#include <vector>
#include <map>
#include <functional>

/**
 * @brief Progress callback function type
 * 
 * Callback signature: void(currentStep, totalSteps, message, percentage)
 * Used to report mining progress to the caller.
 */
using ProgressCallback = std::function<void(int, int, const std::string&, double)>;

/**
 * @brief MaxPRMiner class implementing the maxPR colocation mining algorithm
 * 
 * This class implements the maxPR approach for mining spatial colocation patterns
 * with rare events. The algorithm uses maximal participation ratio (maxPR) measure
 * and exploits weak monotonicity property for efficient pruning.
 */
class MaxPRMiner {
private:
    double minMaxPR;                         ///< Minimum maxPR threshold
    NeighborhoodMgr* neighborhoodMgr;        ///< Pointer to neighborhood manager
    ProgressCallback progressCallback;        ///< Progress reporting callback

    /**
     * @brief Filter star instances that match candidate patterns
     * 
     * Examines star neighborhoods to find instances that match the candidate
     * colocation patterns. This is the first filtering step.
     * 
     * @param candidates Vector of candidate colocation patterns to check
     * @param starNeigh Pair of feature type and its star neighborhoods
     * @return std::vector<ColocationInstance> Filtered instances matching candidates
     */
    std::vector<ColocationInstance> filterStarInstances(
        const std::vector<Colocation>& candidates,
        const std::pair<FeatureType, std::vector<StarNeighborhood>>& starNeigh
    );

    /**
     * @brief Filter clique instances from star instances
     * 
     * Performs clique filtering to ensure all (k-1) subsets of a k-size pattern
     * exist in the previous level. Uses parallel processing for performance.
     * 
     * @param candidates Vector of candidate patterns
     * @param instances Current level star instances
     * @param prevInstances Previous level clique instances
     * @return std::vector<ColocationInstance> Filtered clique instances
     */
    std::vector<ColocationInstance> filterCliqueInstances(
        const std::vector<Colocation>& candidates,
        const std::vector<ColocationInstance>& instances,
        const std::vector<ColocationInstance>& prevInstances
    );

    /**
     * @brief Select colocations based on maximal participation ratio
     * 
     * Calculates the maxPR for each candidate and selects
     * those that meet the minimum maxPR threshold.
     * 
     * @param candidates Vector of candidate patterns
     * @param instances Colocation instances to evaluate
     * @param minMaxPR Minimum maxPR threshold
     * @param featureCount Map of feature type to total instance count
     * @param[out] patternMaxPR Map to store maxPR values for each pattern
     * @return std::vector<Colocation> Colocations meeting maxPR threshold
     */
    std::vector<Colocation> selectMaxPRColocations(
        const std::vector<Colocation>& candidates,
        const std::vector<ColocationInstance>& instances,
        double minMaxPR,
        const std::map<FeatureType, int>& featureCount,
        std::map<Colocation, double>& patternMaxPR
    );

public:
    /**
     * @brief Mine colocation patterns using maxPR algorithm
     * 
     * Main entry point for the mining algorithm. Discovers all colocation
     * patterns with rare events that meet the minimum maxPR threshold.
     * 
     * @param minMaxPR Minimum maximal participation ratio threshold (0.0 to 1.0)
     * @param nbrMgr Pointer to neighborhood manager containing star neighborhoods
     * @param instances Vector of all spatial instances
     * @param progressCb Optional callback for progress reporting
     * @return std::vector<std::pair<Colocation, double>> Patterns with their maxPR values
     */
    std::vector<std::pair<Colocation, double>> mineColocations(
        double minMaxPR, 
        NeighborhoodMgr* nbrMgr, 
        const std::vector<SpatialInstance>& instances,
        ProgressCallback progressCb = nullptr
    );
    
    /**
     * @brief Generate (k+1)-size candidate patterns using weak monotonicity
     * 
     * Uses weak monotonicity property of maxPR: a k-pattern can have at most
     * one (k-1)-subpattern below the threshold. Generates candidates by joining
     * patterns that differ in one feature in their last two features.
     * 
     * @param prevMaxPR Vector of k-size patterns meeting maxPR threshold
     * @param patternMaxPR Map of patterns to their maxPR values
     * @param minMaxPR Minimum maxPR threshold
     * @return std::vector<Colocation> Generated (k+1)-size candidate patterns
     */
    std::vector<Colocation> generateCandidatesMaxPR(
        const std::vector<Colocation>& prevMaxPR,
        const std::map<Colocation, double>& patternMaxPR,
        double minMaxPR
    );
};